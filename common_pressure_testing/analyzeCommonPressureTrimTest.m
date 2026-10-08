function results = analyzeCommonPressureTrimTest(matPath, outputDir)
%ANALYZECOMMONPRESSURETRIMTEST Analyze saved open-loop neutral-trim data.
% Reads logged SimulationOutput only. Does not run Simulink or hardware.

folder = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(matPath)
    matPath = fullfile(folder, 'trim_test_387s.mat');
end
if nargin < 2 || isempty(outputDir)
    [~, runName] = fileparts(matPath);
    outputDir = fullfile(folder, 'results', runName);
end
assert(isfile(matPath), 'Missing MAT-file: %s', matPath);
if ~isfolder(outputDir)
    mkdir(outputDir);
end
loaded = load(matPath, 'out');
assert(isfield(loaded, 'out') && isa(loaded.out, 'Simulink.SimulationOutput'), ...
    'Expected a Simulink.SimulationOutput named out.');
ds = loaded.out.logsout;
assert(isa(ds, 'Simulink.SimulationData.Dataset'), 'Expected out.logsout.');

[t, lengths] = getSignal(ds, 'meas_length');
assert(size(lengths,2) == 3, 'Expected three measured lengths.');
data.t_s = t;
data.length_m = lengths;
data.command_bar = alignSignal(ds, 'des_pressure', t, 'linear');
data.regulatorPath_bar = [ ...
    alignSignal(ds, 'P1 Regulator', t, 'linear'), ...
    alignSignal(ds, 'P2 Regulator', t, 'linear'), ...
    alignSignal(ds, 'P3 Regulator', t, 'linear')];
data.phaseID = round(alignSignal(ds, 'phaseID', t, 'previous'));
data.phi_rad = alignSignal(ds, 'phi', t, 'linear');
data.theta_rad = alignSignal(ds, 'theta', t, 'linear');
data.xyz_m = alignSignal(ds, 'position', t, 'linear');
data.orientationValid = alignSignal(ds, 'orientationValid', t, ...
    'previous', '/Sensor') > 0.5;
data.sensorValid = all(alignSignal(ds, 'sensorValid', t, 'previous') > 0.5, 2);
data.configurationValid = alignSignal(ds, 'configurationValid', ...
    t, 'previous') > 0.5;
data.pollStart_s = alignSignal(ds, 'pollStart_s', t, 'linear');
[~, kill] = getSignal(ds, 'kitt_switch');
data.loggedKillState = kill(1); % Raw logged dashboard state only.

radius_m = 0.013; % PMA radial offset from backbone.
data.U_rad = (lengths(:,2) + lengths(:,3))/radius_m;
data.V_rad = (lengths(:,3) - lengths(:,2))/(sqrt(3)*radius_m);
data.valid = data.orientationValid & data.sensorValid ...
    & data.configurationValid & all(isfinite(lengths),2) ...
    & isfinite(data.phi_rad) & all(isfinite(data.xyz_m),2);

holdIDs = [0; 102;112;114;122;124;132;134;142;144; ...
    202;212;214;222;224;232;234;242;244; ...
    302;312;314;322;324;332;334;342;344; 901];
phaseSummary = summarizeHolds(data, holdIDs, 5.0);
zeroRows = phaseSummary.phaseID == 0 | phaseSummary.phaseID == 901;
targetUV_rad = [mean(phaseSummary.U_rad(zeroRows)); ...
    mean(phaseSummary.V_rad(zeroRows))];
trimSummary = estimateTrim(phaseSummary, targetUV_rad);

results.sourceMat = char(matPath);
results.outputDir = char(outputDir);
results.data = data;
results.phaseSummary = phaseSummary;
results.trimSummary = trimSummary;
results.targetUV_rad = targetUV_rad;
results.assumptions.radius_m = radius_m;
results.assumptions.tailWindow_s = 5.0;
results.assumptions.commandFloorForFeasibility_bar = 1.0;
results.assumptions.commandCeiling_bar = 3.0;
writetable(phaseSummary, fullfile(outputDir, 'phase_summary.csv'));
writetable(trimSummary, fullfile(outputDir, 'trim_estimates.csv'));
writeReport(results, fullfile(outputDir, 'report.txt'));
makePlots(results, outputDir);
save(fullfile(outputDir, 'analysis_results.mat'), 'results', '-v7');
end

function summary = summarizeHolds(data, ids, tailSeconds)
n = numel(ids);
phaseID = ids(:);
level_bar = NaN(n,1); start_s = NaN(n,1); end_s = NaN(n,1);
validFraction = NaN(n,1); nTail = zeros(n,1);
p1_bar = NaN(n,1); p2_bar = NaN(n,1); p3_bar = NaN(n,1);
regulatorError_bar = NaN(n,1);
U_rad = NaN(n,1); V_rad = NaN(n,1); phi_rad = NaN(n,1);
URange_rad = NaN(n,1); VRange_rad = NaN(n,1);
l1_mm = NaN(n,1); l2_mm = NaN(n,1); l3_mm = NaN(n,1);
x_mm = NaN(n,1); y_mm = NaN(n,1); z_mm = NaN(n,1);
for j = 1:n
    inPhase = data.phaseID == ids(j);
    assert(any(inPhase), 'Missing logged phase %d.', ids(j));
    start_s(j) = min(data.t_s(inPhase));
    end_s(j) = max(data.t_s(inPhase));
    validFraction(j) = mean(data.valid(inPhase));
    tail = inPhase & data.t_s >= end_s(j) - tailSeconds;
    good = tail & data.valid;
    nTail(j) = nnz(good);
    assert(nTail(j) >= 20, 'Too few valid tail samples in phase %d.', ids(j));
    p = median(data.command_bar(tail,:),1,'omitnan');
    p1_bar(j) = p(1); p2_bar(j) = p(2); p3_bar(j) = p(3);
    level_bar(j) = mean(p);
    regulatorError_bar(j) = max(abs(data.regulatorPath_bar(tail,:) ...
        - data.command_bar(tail,:)),[],'all');
    U_rad(j) = mean(data.U_rad(good),'omitnan');
    V_rad(j) = mean(data.V_rad(good),'omitnan');
    phi_rad(j) = mean(data.phi_rad(good),'omitnan');
    URange_rad(j) = range(data.U_rad(good));
    VRange_rad(j) = range(data.V_rad(good));
    l = mean(data.length_m(good,:),1,'omitnan')*1e3;
    l1_mm(j)=l(1); l2_mm(j)=l(2); l3_mm(j)=l(3);
    xyz = mean(data.xyz_m(good,:),1,'omitnan')*1e3;
    x_mm(j)=xyz(1); y_mm(j)=xyz(2); z_mm(j)=xyz(3);
end
summary = table(phaseID,level_bar,start_s,end_s,validFraction,nTail, ...
    p1_bar,p2_bar,p3_bar,regulatorError_bar,U_rad,V_rad,phi_rad, ...
    URange_rad,VRange_rad,l1_mm,l2_mm,l3_mm,x_mm,y_mm,z_mm);
end

function summary = estimateTrim(phases, targetUV_rad)
common_bar = [1.2;1.6;2.0];
n = numel(common_bar);
baselineU_rad = NaN(n,1); baselineV_rad = NaN(n,1);
targetU_rad = repmat(targetUV_rad(1),n,1);
targetV_rad = repmat(targetUV_rad(2),n,1);
baselineUDrift_rad = NaN(n,1); baselineVDrift_rad = NaN(n,1);
baselineUSpan_rad = NaN(n,1); baselineVSpan_rad = NaN(n,1);
Js_U_radPerBar = NaN(n,1); Js_V_radPerBar = NaN(n,1);
Ja_U_radPerBar = NaN(n,1); Ja_V_radPerBar = NaN(n,1);
symDirectionDisagreement_radPerBar = NaN(n,1);
antiDirectionDisagreement_radPerBar = NaN(n,1);
JacobianCondition = NaN(n,1);
candidateS_bar = NaN(n,1); candidateA_bar = NaN(n,1);
candidateP1_bar = NaN(n,1); candidateP2_bar = NaN(n,1);
candidateP3_bar = NaN(n,1);
withinTestedPerturbations = false(n,1);
withinNominalBounds = false(n,1);
for j = 1:n
    base = 100*j;
    h = @(suffix) holdUV(phases, base+suffix);
    B0=h(2); B1=h(14); B2=h(24); B3=h(34); B4=h(44);
    baselines = [B0 B1 B2 B3 B4];
    y0 = mean(baselines,2);
    baselineU_rad(j)=y0(1); baselineV_rad(j)=y0(2);
    baselineUSpan_rad(j)=range(baselines(1,:));
    baselineVSpan_rad(j)=range(baselines(2,:));
    drift=B4-B0;
    baselineUDrift_rad(j)=drift(1);
    baselineVDrift_rad(j)=drift(2);
    symPlus=(h(12)-(B0+B1)/2)/0.2;
    symMinus=(h(22)-(B1+B2)/2)/(-0.2);
    antiPlus=(h(32)-(B2+B3)/2)/0.1;
    antiMinus=(h(42)-(B3+B4)/2)/(-0.1);
    J=[(symPlus+symMinus)/2, (antiPlus+antiMinus)/2];
    Js_U_radPerBar(j)=J(1,1); Js_V_radPerBar(j)=J(2,1);
    Ja_U_radPerBar(j)=J(1,2); Ja_V_radPerBar(j)=J(2,2);
    symDirectionDisagreement_radPerBar(j)=norm(symPlus-symMinus);
    antiDirectionDisagreement_radPerBar(j)=norm(antiPlus-antiMinus);
    JacobianCondition(j)=cond(J);
    if all(isfinite(J),'all') && rcond(J)>1e-3
        correction=-J\(y0-targetUV_rad);
        s=correction(1); a=correction(2);
        p=[common_bar(j)-2*s/3; common_bar(j)+s/3-a; ...
            common_bar(j)+s/3+a];
        candidateS_bar(j)=s; candidateA_bar(j)=a;
        candidateP1_bar(j)=p(1); candidateP2_bar(j)=p(2);
        candidateP3_bar(j)=p(3);
        withinTestedPerturbations(j)=abs(s)<=0.2 && abs(a)<=0.1;
        withinNominalBounds(j)=all(p>=1.0) && all(p<=3.0);
    end
end
summary=table(common_bar,targetU_rad,targetV_rad, ...
    baselineU_rad,baselineV_rad, ...
    baselineUDrift_rad,baselineVDrift_rad, ...
    baselineUSpan_rad,baselineVSpan_rad, ...
    Js_U_radPerBar,Js_V_radPerBar, ...
    Ja_U_radPerBar,Ja_V_radPerBar, ...
    symDirectionDisagreement_radPerBar, ...
    antiDirectionDisagreement_radPerBar,JacobianCondition, ...
    candidateS_bar,candidateA_bar,candidateP1_bar, ...
    candidateP2_bar,candidateP3_bar, ...
    withinTestedPerturbations,withinNominalBounds);
end

function y = holdUV(phases, id)
ix=find(phases.phaseID==id);
assert(numel(ix)==1, 'Expected one summary row for phase %d.', id);
y=[phases.U_rad(ix); phases.V_rad(ix)];
end

function writeReport(results, path)
d=results.data; p=results.phaseSummary; s=results.trimSummary;
file=fopen(path,'w'); assert(file>=0,'Could not write report.');
cleanup=onCleanup(@()fclose(file)); %#ok<NASGU>
fprintf(file,'Neutral-trim test: %s\n',results.sourceMat);
fprintf(file,'Simulation %.1f..%.1f s; NDI samples %d; valid %d (%.1f%%)\n', ...
    d.t_s(1),d.t_s(end),numel(d.t_s),nnz(d.valid),100*mean(d.valid));
poll=diff(d.pollStart_s); poll=poll(isfinite(poll)&poll>0);
fprintf(file,'Physical NDI poll interval median %.4f s; span %.2f s\n', ...
    median(poll),d.pollStart_s(end)-d.pollStart_s(1));
fprintf(file,'Logged raw Kill dashboard state %.0f\n',d.loggedKillState);
fprintf(file,'Maximum command versus regulator-path logged difference %.5f bar\n', ...
    max(abs(d.command_bar-d.regulatorPath_bar),[],'all'));
fprintf(file,'These are software commands, not measured chamber pressures.\n');
fprintf(file,'Hold summaries average each final 5 simulation seconds.\n');
fprintf(file,'Trim target is the mean of initial/final zero-pressure U,V: %.5f, %.5f rad.\n', ...
    results.targetUV_rad(1),results.targetUV_rad(2));
fprintf(file,'Maximum hold U range %.5f rad; V range %.5f rad.\n\n', ...
    max(p.URange_rad),max(p.VRange_rad));
fprintf(file,'Common  baseline U,V (rad)  first-last U,V drift  baseline U,V span\n');
for j=1:height(s)
    fprintf(file,'%4.1f    %+7.4f,%+7.4f   %+7.4f,%+7.4f   %7.4f,%7.4f\n', ...
        s.common_bar(j),s.baselineU_rad(j),s.baselineV_rad(j), ...
        s.baselineUDrift_rad(j),s.baselineVDrift_rad(j), ...
        s.baselineUSpan_rad(j),s.baselineVSpan_rad(j));
end
fprintf(file,'\nCommon     Js=[dU,dV]/ds       Ja=[dU,dV]/da\n');
for j=1:height(s)
    fprintf(file,'%4.1f    %+7.4f,%+7.4f     %+7.4f,%+7.4f\n', ...
        s.common_bar(j), ...
        s.Js_U_radPerBar(j),s.Js_V_radPerBar(j), ...
        s.Ja_U_radPerBar(j),s.Ja_V_radPerBar(j));
end
fprintf(file,'\nLinearized neutral-trim candidates (NOT commands to replay without validation):\n');
fprintf(file,'Common   s(bar)  a(bar)    P1,P2,P3(bar)         cond(J)  within tested range  1..3 bar\n');
for j=1:height(s)
    fprintf(file,'%4.1f   %+7.3f %+7.3f   %6.3f,%6.3f,%6.3f   %8.2f       %d                %d\n', ...
        s.common_bar(j),s.candidateS_bar(j),s.candidateA_bar(j), ...
        s.candidateP1_bar(j),s.candidateP2_bar(j),s.candidateP3_bar(j), ...
        s.JacobianCondition(j),s.withinTestedPerturbations(j), ...
        s.withinNominalBounds(j));
end
fprintf(file,'\nJ is the local command-pressure-to-orientation response, not a force map.\n');
fprintf(file,'Its estimate compares +/- pulses with adjacent baseline holds.\n');
fprintf(file,'Baseline drift, loading/unloading differences, and extrapolation\n');
fprintf(file,'limit use of these candidates. Verify any feasible trim in a new hold.\n');
end

function makePlots(results, outputDir)
d=results.data; s=results.trimSummary;
f=figure('Visible','off','Position',[80 80 1300 850]);
tt=tiledlayout(f,3,1,'TileSpacing','compact');
ax=nexttile(tt); plot(ax,d.t_s,d.command_bar,'LineWidth',1.0);
ylabel(ax,'Command (bar)'); grid(ax,'on'); ylim(ax,[0 3]);
legend(ax,{'P1','P2','P3'},'Location','northwest');
ax=nexttile(tt); plot(ax,d.t_s,d.U_rad,d.t_s,d.V_rad,'LineWidth',1.0);
ylabel(ax,'Signed bend (rad)'); grid(ax,'on');
legend(ax,{'U','V'},'Location','northwest');
ax=nexttile(tt); plot(ax,d.t_s,d.phi_rad,'LineWidth',1.0);
ylabel(ax,'phi (rad)'); xlabel(ax,'Simulation time (s)'); grid(ax,'on');
title(tt,'Open-loop pressure and measured bending');
exportgraphics(f,fullfile(outputDir,'trim_trajectories.png'),'Resolution',160);
close(f);

f=figure('Visible','off','Position',[80 80 1100 700]);
hold on; grid on; axis equal;
plot(d.U_rad(d.valid),d.V_rad(d.valid),'Color',[0.7 0.7 0.7]);
plot(s.baselineU_rad,s.baselineV_rad,'ko','MarkerFaceColor','k');
text(s.baselineU_rad,s.baselineV_rad, ...
    compose('  %.1f bar',s.common_bar));
plot(results.targetUV_rad(1),results.targetUV_rad(2), ...
    'rx','MarkerSize',12,'LineWidth',2);
xlabel('U (rad)'); ylabel('V (rad)');
title('Measured bend trajectory, zero-pressure target, and baseline means');
exportgraphics(f,fullfile(outputDir,'bend_plane.png'),'Resolution',160);
close(f);
end

function [t, x] = getSignal(ds, name, blockSuffix)
if nargin<3, blockSuffix=''; end
names=ds.getElementNames;
candidates=find(strcmp(names,name));
matches=zeros(0,1);
for i=candidates(:)'
    element=ds.get(i);
    path=char(element.BlockPath.getBlock(1));
    if isempty(blockSuffix) || endsWith(path,blockSuffix)
        matches(end+1,1)=i; %#ok<AGROW>
    end
end
assert(numel(matches)==1, ...
    'Expected one logged %s matching "%s"; found %d.', ...
    name,blockSuffix,numel(matches));
values=ds.get(matches).Values;
assert(isa(values,'timeseries'),'Expected timeseries for %s.',name);
t=double(values.Time(:));
raw=double(values.Data);
if size(raw,1)==numel(t) && ndims(raw)<=2
    x=reshape(raw,numel(t),[]);
else
    x=reshape(raw,[],numel(t)).';
end
assert(size(x,1)==numel(t),'Time/data mismatch for %s.',name);
end

function x = alignSignal(ds,name,t,method,blockSuffix)
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
