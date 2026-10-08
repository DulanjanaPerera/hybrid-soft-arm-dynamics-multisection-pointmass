function results = analyzeCommonPressureTrimValidation(matPath, outputDir)
%ANALYZECOMMONPRESSURETRIMVALIDATION Analyze saved 199 s trim test only.
% No Simulink model or NI/NDI device is opened or run.

folder = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(matPath)
    matPath = fullfile(folder, 'trim_test_199s.mat');
end
if nargin < 2 || isempty(outputDir)
    [~, runName] = fileparts(matPath);
    outputDir = fullfile(folder, 'results', runName);
end
assert(isfile(matPath), 'Missing MAT-file: %s', matPath);
if ~isfolder(outputDir), mkdir(outputDir); end
loaded = load(matPath, 'out');
assert(isfield(loaded,'out') && isa(loaded.out,'Simulink.SimulationOutput'), ...
    'Expected Simulink.SimulationOutput named out.');
ds = loaded.out.logsout;
assert(isa(ds,'Simulink.SimulationData.Dataset'),'Expected out.logsout.');

[t, length_m] = getSignal(ds,'meas_length');
assert(size(length_m,2)==3,'Expected three measured lengths.');
d.t_s = t;
d.length_m = length_m;
d.command_bar = alignSignal(ds,'des_pressure',t,'linear');
d.regulatorPath_bar = [ ...
    alignSignal(ds,'P1 Regulator',t,'linear'), ...
    alignSignal(ds,'P2 Regulator',t,'linear'), ...
    alignSignal(ds,'P3 Regulator',t,'linear')];
d.phaseID = round(alignSignal(ds,'phaseID',t,'previous'));
d.phi_rad = alignSignal(ds,'phi',t,'linear');
d.theta_rad = alignSignal(ds,'theta',t,'linear');
d.xyz_m = alignSignal(ds,'position',t,'linear');
d.orientationValid = alignSignal(ds,'orientationValid',t,'previous','/Sensor')>0.5;
d.sensorValid = all(alignSignal(ds,'sensorValid',t,'previous')>0.5,2);
d.configurationValid = alignSignal(ds,'configurationValid',t,'previous')>0.5;
d.pollStart_s = alignSignal(ds,'pollStart_s',t,'linear');
[~,kill] = getSignal(ds,'kitt_switch');
d.loggedKillState = kill(1);

radius_m = 0.013;
d.U_rad = (length_m(:,2)+length_m(:,3))/radius_m;
d.V_rad = (length_m(:,3)-length_m(:,2))/(sqrt(3)*radius_m);
d.valid = d.orientationValid & d.sensorValid & d.configurationValid ...
    & all(isfinite(length_m),2) & isfinite(d.phi_rad) ...
    & all(isfinite(d.xyz_m),2);

holdIDs = [12;22;32;42;52;62;72;82;92;102;112];
holds = summarizeHolds(d,holdIDs);
pick = @(id) find(holds.phaseID==id,1);
iA = pick(62); iB = pick(102);
iIntermediate = pick(52);
i0 = pick(12); iMid = pick(82); iEnd = pick(112);
assert(all([iA iB iIntermediate i0 iMid iEnd]>0), ...
    'Missing comparison phase.');

routeARef = mean([holds.U_rad([i0 iMid]),holds.V_rad([i0 iMid])],1);
routeBRef = mean([holds.U_rad([iMid iEnd]),holds.V_rad([iMid iEnd])],1);
candidateA = [holds.U_rad(iA),holds.V_rad(iA)];
candidateB = [holds.U_rad(iB),holds.V_rad(iB)];
results.sourceMat = char(matPath);
results.outputDir = char(outputDir);
results.data = d;
results.holds = holds;
results.routeAReferenceUV_rad = routeARef;
results.routeBReferenceUV_rad = routeBRef;
results.candidateAUV_rad = candidateA;
results.candidateBUV_rad = candidateB;
results.routeAResidualUV_rad = candidateA-routeARef;
results.routeBResidualUV_rad = candidateB-routeBRef;
results.candidateDifferenceUV_rad = candidateB-candidateA;
results.intermediateResidualUV_rad = ...
    [holds.U_rad(iIntermediate),holds.V_rad(iIntermediate)]-routeARef;
results.zeroReferenceDifferenceUV_rad = ...
    [holds.U_rad(iEnd)-holds.U_rad(i0), ...
     holds.V_rad(iEnd)-holds.V_rad(i0)];
results.candidateLengthDifference_mm = [ ...
    holds.l1_mm(iB)-holds.l1_mm(iA), ...
    holds.l2_mm(iB)-holds.l2_mm(iA), ...
    holds.l3_mm(iB)-holds.l3_mm(iA)];
xyz = [holds.x_mm,holds.y_mm,holds.z_mm];
results.routeAReferenceXYZ_mm = mean(xyz([i0 iMid],:),1);
results.routeBReferenceXYZ_mm = mean(xyz([iMid iEnd],:),1);
results.routeAResidualXYZ_mm = xyz(iA,:)-results.routeAReferenceXYZ_mm;
results.routeBResidualXYZ_mm = xyz(iB,:)-results.routeBReferenceXYZ_mm;
results.intermediateResidualXYZ_mm = ...
    xyz(iIntermediate,:)-results.routeAReferenceXYZ_mm;

writetable(holds,fullfile(outputDir,'hold_summary.csv'));
writeReport(results,fullfile(outputDir,'report.txt'));
makePlots(results,outputDir);
save(fullfile(outputDir,'analysis_results.mat'),'results','-v7');
end

function holds = summarizeHolds(d,ids)
n=numel(ids);
phaseID=ids(:); start_s=NaN(n,1); end_s=NaN(n,1);
validFraction=NaN(n,1); nTail=zeros(n,1);
p1_bar=NaN(n,1); p2_bar=NaN(n,1); p3_bar=NaN(n,1);
regulatorError_bar=NaN(n,1);
U_rad=NaN(n,1); V_rad=NaN(n,1); phi_rad=NaN(n,1);
URange_rad=NaN(n,1); VRange_rad=NaN(n,1);
USlope_radPerS=NaN(n,1); VSlope_radPerS=NaN(n,1);
l1_mm=NaN(n,1); l2_mm=NaN(n,1); l3_mm=NaN(n,1);
x_mm=NaN(n,1); y_mm=NaN(n,1); z_mm=NaN(n,1);
for j=1:n
    inPhase=d.phaseID==ids(j);
    assert(any(inPhase),'Missing phase %d.',ids(j));
    start_s(j)=min(d.t_s(inPhase));
    end_s(j)=max(d.t_s(inPhase));
    validFraction(j)=mean(d.valid(inPhase));
    tail=inPhase & d.t_s>=end_s(j)-5.0;
    good=tail & d.valid;
    nTail(j)=nnz(good);
    assert(nTail(j)>=20,'Too few valid samples in phase %d.',ids(j));
    p=median(d.command_bar(tail,:),1,'omitnan');
    p1_bar(j)=p(1); p2_bar(j)=p(2); p3_bar(j)=p(3);
    regulatorError_bar(j)=max(abs(d.command_bar(tail,:) ...
        -d.regulatorPath_bar(tail,:)),[],'all');
    U_rad(j)=mean(d.U_rad(good),'omitnan');
    V_rad(j)=mean(d.V_rad(good),'omitnan');
    phi_rad(j)=mean(d.phi_rad(good),'omitnan');
    URange_rad(j)=range(d.U_rad(good));
    VRange_rad(j)=range(d.V_rad(good));
    tt=d.t_s(good); tt=tt-mean(tt);
    slope=polyfit(tt,d.U_rad(good),1); USlope_radPerS(j)=slope(1);
    slope=polyfit(tt,d.V_rad(good),1); VSlope_radPerS(j)=slope(1);
    l=mean(d.length_m(good,:),1,'omitnan')*1e3;
    l1_mm(j)=l(1); l2_mm(j)=l(2); l3_mm(j)=l(3);
    xyz=mean(d.xyz_m(good,:),1,'omitnan')*1e3;
    x_mm(j)=xyz(1); y_mm(j)=xyz(2); z_mm(j)=xyz(3);
end
holds=table(phaseID,start_s,end_s,validFraction,nTail, ...
    p1_bar,p2_bar,p3_bar,regulatorError_bar, ...
    U_rad,V_rad,phi_rad,URange_rad,VRange_rad, ...
    USlope_radPerS,VSlope_radPerS, ...
    l1_mm,l2_mm,l3_mm,x_mm,y_mm,z_mm);
end

function writeReport(r,path)
d=r.data; h=r.holds;
file=fopen(path,'w'); assert(file>=0,'Could not write report.');
cleanup=onCleanup(@()fclose(file)); %#ok<NASGU>
fprintf(file,'Common-pressure trim validation: %s\n',r.sourceMat);
fprintf(file,'Simulation %.2f..%.2f s; NDI %d samples, %d valid (%.1f%%)\n', ...
    d.t_s(1),d.t_s(end),numel(d.t_s),nnz(d.valid),100*mean(d.valid));
poll=diff(d.pollStart_s); poll=poll(isfinite(poll)&poll>0);
fprintf(file,'Physical NDI poll median interval %.4f s; span %.2f s\n', ...
    median(poll),d.pollStart_s(end)-d.pollStart_s(1));
fprintf(file,'Raw logged Kill dashboard state %.0f\n',d.loggedKillState);
fprintf(file,'Maximum command vs regulator-path logged difference %.5f bar\n', ...
    max(abs(d.command_bar-d.regulatorPath_bar),[],'all'));
fprintf(file,'The pressure signals are software commands, not chamber measurements.\n');
fprintf(file,'Each hold row averages its final 5 simulation seconds.\n\n');
fprintf(file,'Phase  Pressure P1,P2,P3(bar)    U(rad)    V(rad)   phi(rad)   U range   V range   U slope   V slope\n');
for j=1:height(h)
    fprintf(file,'%3d    %5.3f,%5.3f,%5.3f    %+8.5f %+8.5f  %8.5f   %.5f   %.5f   %+.5f  %+.5f\n', ...
        h.phaseID(j),h.p1_bar(j),h.p2_bar(j),h.p3_bar(j), ...
        h.U_rad(j),h.V_rad(j),h.phi_rad(j), ...
        h.URange_rad(j),h.VRange_rad(j), ...
        h.USlope_radPerS(j),h.VSlope_radPerS(j));
end
fprintf(file,'\nRoute A zero-reference U,V: %.5f, %.5f rad\n',r.routeAReferenceUV_rad);
fprintf(file,'Route B zero-reference U,V: %.5f, %.5f rad\n',r.routeBReferenceUV_rad);
fprintf(file,'Candidate A residual U,V: %+8.5f,%+8.5f rad\n',r.routeAResidualUV_rad);
fprintf(file,'Candidate B residual U,V: %+8.5f,%+8.5f rad\n',r.routeBResidualUV_rad);
fprintf(file,'Candidate B-A difference U,V: %+8.5f,%+8.5f rad\n',r.candidateDifferenceUV_rad);
fprintf(file,'Intermediate [1.7,1.7,1.4] residual U,V: %+8.5f,%+8.5f rad\n', ...
    r.intermediateResidualUV_rad);
fprintf(file,'Candidate B-A length differences l1,l2,l3: %+7.3f,%+7.3f,%+7.3f mm\n', ...
    r.candidateLengthDifference_mm);
fprintf(file,'Tip XYZ residual norm: intermediate %.2f mm; candidate A %.2f mm; candidate B %.2f mm\n', ...
    norm(r.intermediateResidualXYZ_mm), ...
    norm(r.routeAResidualXYZ_mm),norm(r.routeBResidualXYZ_mm));
fprintf(file,'Initial-to-final zero-reference difference U,V: %+8.5f,%+8.5f rad\n', ...
    r.zeroReferenceDifferenceUV_rad);
fprintf(file,'Candidate comparison isolates repeatability only to the extent\n');
fprintf(file,'that the two final-five-second holds are settled and reference\n');
fprintf(file,'conditions match. Pressure history, valve behavior, and model\n');
fprintf(file,'error are not separately identified without chamber sensors.\n');
end

function makePlots(r,outputDir)
d=r.data; h=r.holds;
f=figure('Visible','off','Position',[80 80 1300 900]);
tt=tiledlayout(f,3,1,'TileSpacing','compact');
ax=nexttile(tt);
plot(ax,d.t_s,d.command_bar,'LineWidth',1.0); grid(ax,'on');
ylabel(ax,'Command (bar)'); ylim(ax,[0 3]);
legend(ax,{'P1','P2','P3'},'Location','northwest');
ax=nexttile(tt);
plot(ax,d.t_s,d.U_rad,d.t_s,d.V_rad,'LineWidth',1.0); grid(ax,'on');
ylabel(ax,'Signed bend (rad)');
lg=legend(ax,{'U','V'},'Location','northwest');
lg.AutoUpdate='off';
hold(ax,'on');
for id=[62 102]
    i=find(h.phaseID==id,1);
    xline(ax,h.start_s(i),'k:');
    xline(ax,h.end_s(i),'k:');
end
ax=nexttile(tt);
plot(ax,d.t_s,d.phi_rad,'LineWidth',1.0); grid(ax,'on');
ylabel(ax,'phi (rad)'); xlabel(ax,'Simulation time (s)');
title(tt,'Trim validation: two approaches to the same pressure triplet');
exportgraphics(f,fullfile(outputDir,'validation_trajectories.png'),'Resolution',160);
close(f);

f=figure('Visible','off','Position',[80 80 1100 700]);
plot(d.U_rad(d.valid),d.V_rad(d.valid),'Color',[0.72 0.72 0.72]);
hold on; grid on; axis equal;
plot(r.routeAReferenceUV_rad(1),r.routeAReferenceUV_rad(2), ...
    'kx','MarkerSize',11,'LineWidth',2);
plot(r.candidateAUV_rad(1),r.candidateAUV_rad(2), ...
    'bo','MarkerSize',9,'MarkerFaceColor','b');
plot(r.candidateBUV_rad(1),r.candidateBUV_rad(2), ...
    'ro','MarkerSize',9,'MarkerFaceColor','r');
xlabel('U (rad)'); ylabel('V (rad)');
legend({'Trajectory','Zero reference','Candidate A','Candidate B'}, ...
    'Location','best');
title('Candidate repeatability in signed bending coordinates');
exportgraphics(f,fullfile(outputDir,'candidate_repeatability.png'),'Resolution',160);
close(f);
end

function [t,x]=getSignal(ds,name,blockSuffix)
if nargin<3, blockSuffix=''; end
names=ds.getElementNames;
candidates=find(strcmp(names,name));
matches=zeros(0,1);
for i=candidates(:)'
    e=ds.get(i);
    path=char(e.BlockPath.getBlock(1));
    if isempty(blockSuffix) || endsWith(path,blockSuffix)
        matches(end+1,1)=i; %#ok<AGROW>
    end
end
assert(numel(matches)==1, ...
    'Expected one %s matching "%s"; found %d.',name,blockSuffix,numel(matches));
v=ds.get(matches).Values;
assert(isa(v,'timeseries'),'Expected timeseries for %s.',name);
t=double(v.Time(:));
raw=double(v.Data);
if size(raw,1)==numel(t) && ndims(raw)<=2
    x=reshape(raw,numel(t),[]);
else
    x=reshape(raw,[],numel(t)).';
end
assert(size(x,1)==numel(t),'Time/data mismatch for %s.',name);
end

function x=alignSignal(ds,name,t,method,blockSuffix)
if nargin<5, blockSuffix=''; end
[sourceTime,source]=getSignal(ds,name,blockSuffix);
[sourceTime,keep]=unique(sourceTime,'last');
source=source(keep,:);
if numel(sourceTime)==1
    x=repmat(source,numel(t),1);
else
    x=interp1(sourceTime,source,t,method,NaN);
end
end
