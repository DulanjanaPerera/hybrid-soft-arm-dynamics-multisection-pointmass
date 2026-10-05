function results = analyzeInverseDynamics(inputFile)
% Offline analysis only. Does not load Simulink models or access hardware.
% Run from any folder: addpath('.../inverDynamic_test'); analyzeInverseDynamics
here = fileparts(mfilename('fullpath'));
if nargin < 1, inputFile = fullfile(here,'simulationresulta_20261001_1936.mat'); end
root = fileparts(here); addpath(root);
outputDir = fullfile(here,'results'); if ~isfolder(outputDir), mkdir(outputDir); end
cfg.MinHold_s=5; cfg.Tail_s=3; cfg.TargetTolerance_m=1e-7;
cfg.SettlingTolerance_m=0.001; cfg.SettlingDwell_s=2;
cfg.StationaryVelocity_m_s=0.0002; cfg.StationarySpread_m=0.0005;
cfg.ThetaMinPhi_rad=0.5*pi/180;
s=load(inputFile,'out'); ds=s.out.logsout;
[t,qdes]=signal(ds,'sampled des_length');
q=aligned(ds,'lengthChnage_m',t,'linear'); q=q(:,2:3);
dq=aligned(ds,'filtered_d_length',t,'linear');
p=aligned(ds,'pressureCommand_bar',t,'previous');
tau=aligned(ds,'tauDesired_N',t,'previous');
available=aligned(ds,'tauAvailable_N',t,'previous');
limited=logical(aligned(ds,'pressureLimited',t,'previous'));
valid=logical(aligned(ds,'orientationValid',t,'previous')) & logical(aligned(ds,'inverseModelValid',t,'previous'));
physical=aligned(ds,'pollStart',t,'linear'); physical=physical-physical(1);
assert(all(diff(physical)>0),'Physical poll times must increase strictly.');
phi=aligned(ds,'phiOrientation',t,'linear'); theta=aligned(ds,'thetaOrientation',t,'linear');
observable=logical(aligned(ds,'thetaObservable',t,'previous'));
desPhi=aligned(ds,'des_phi',t,'previous'); desTheta=aligned(ds,'des_theta',t,'previous');
desVelocity=aligned(ds,'dl',t,'previous'); desAcceleration=aligned(ds,'ddl',t,'previous');
assert(all(abs(desVelocity(:))<1e-12) && all(abs(desAcceleration(:))<1e-12), ...
    'This static analysis requires logged desired dl and ddl to be zero.');
% Target holds are defined from the actual sampled length reference.
starts=[1;find(any(abs(diff(qdes,1,1))>cfg.TargetTolerance_m,2))+1];
ends=[starts(2:end)-1;numel(t)];
params.N=1; params.L=.17411; params.r=.013; params.mi=.1;
params.g=[0;0;-9.81]; params.K=1350*[2 1;1 2];
params.mu=2000; params.lKbounds=[-.03;.03;1e6];
pressureRadius_m=.0065; B=pi*pressureRadius_m^2*1e5*[-1 1 0;-1 0 1];
rows=struct([]);
for k=1:numel(starts)
 a=starts(k); b=ends(k); duration=physical(b)-physical(a);
 if duration<cfg.MinHold_s, continue; end
 ix=(a:b)'; tail=ix(physical(ix)>=physical(b)-cfg.Tail_s); good=tail(valid(tail)&all(isfinite(q(tail,:)),2)&all(isfinite(dq(tail,:)),2));
 if isempty(good), continue; end
 qm=mean(q(good,:),1); target=qdes(a,:); err=target-qm;
 spread=max(q(good,:),[],1)-min(q(good,:),[],1);
 speed=sqrt(sum(dq(good,:).^2,2));
 stationary=max(speed)<=cfg.StationaryVelocity_m_s && max(spread)<=cfg.StationarySpread_m;
 [elastic,gravity]=staticForce(qm,params); required=elastic+gravity;
 [ed,gd]=staticForce(target,params); predicted=ed+gd;
 effective=max(p(good,:)-.8,0); mapped=effective*B';
 % Settling requires both coordinates in tolerance continuously to hold end.
 inside=valid(ix)&all(abs(q(ix,:)-target)<=cfg.SettlingTolerance_m,2);
 lastOutside=find(~inside,1,'last'); if isempty(lastOutside), firstInside=1; else, firstInside=lastOutside+1; end
 settle=NaN;
 if firstInside<=numel(ix) && physical(b)-physical(ix(firstInside))>=cfg.SettlingDwell_s
  settle=physical(ix(firstInside))-physical(a);
 end
 thetaGood=good(observable(good)&phi(good)>cfg.ThetaMinPhi_rad&desPhi(good)>cfg.ThetaMinPhi_rad);
 thetaError=NaN; if ~isempty(thetaGood), thetaError=atan2(mean(sin(desTheta(thetaGood)-theta(thetaGood))),mean(cos(desTheta(thetaGood)-theta(thetaGood)))); end
 r.Hold=k; r.SimStart_s=t(a); r.PhysicalStart_s=physical(a); r.Duration_s=duration;
 r.DesTheta_rad=desTheta(a); r.DesPhi_rad=desPhi(a);
 r.MeasPhi_rad=mean(phi(good)); r.ThetaError_rad=thetaError;
 r.DesL2_mm=1000*target(1); r.DesL3_mm=1000*target(2);
 r.MeasL2_mm=1000*qm(1); r.MeasL3_mm=1000*qm(2);
 r.ErrorL2_mm=1000*err(1); r.ErrorL3_mm=1000*err(2);
 r.TailSpread_mm=1000*max(spread); r.TailMaxSpeed_mm_s=1000*max(speed);
 r.StationaryTail=stationary; r.ValidFraction=mean(valid(ix));
 r.LimitedFraction=mean(limited(ix)); r.SettlingTime_s=settle;
 r.MaxTrackingError_mm=1000*max(vecnorm(q(ix(valid(ix)),:)-target,2,2));
 pm=mean(p(good,:),1); r.P1_bar=pm(1); r.P2_bar=pm(2); r.P3_bar=pm(3);
 av=mean(available(good,:),1); residual=required-av;
 r.Elastic2_N=elastic(1); r.Elastic3_N=elastic(2); r.Gravity2_N=gravity(1); r.Gravity3_N=gravity(2);
 r.Required2_N=required(1); r.Required3_N=required(2);
 r.Available2_N=av(1); r.Available3_N=av(2);
 r.Residual2_N=residual(1); r.Residual3_N=residual(2);
 r.DesiredForceCheck_N=norm(mean(tau(good,:),1)-predicted);
 r.PressureMappingCheck_N=max(vecnorm(mapped-available(good,:),2,2));
 if isempty(rows), rows=r; else, rows(end+1)=r; end %#ok<AGROW>
end
if isempty(rows), error('No valid holds meet minimum duration.'); end
holds=struct2table(rows);
results.InputFile=inputFile; results.Config=cfg; results.ModelParameters=params;
results.PressureRadius_m=pressureRadius_m; results.Holds=holds;
results.SampleCount=numel(t); results.ValidFraction=mean(valid);
results.PhysicalDuration_s=physical(end); results.SimulationDuration_s=t(end)-t(1);
results.PollMedian_s=median(diff(physical)); results.PollMax_s=max(diff(physical));
results.LimitedFraction=mean(limited);
results.TimeSimulation_s=t; results.TimePhysical_s=physical;
results.DesiredLengths_m=qdes; results.MeasuredLengths_m=q; results.PressureCommand_bar=p;
writetable(holds,fullfile(outputDir,'hold_summary.csv'));
save(fullfile(outputDir,'analysis_results.mat'),'results');
f=figure('Visible','off','Position',[50 50 1100 850]);
tiledlayout(4,1);
nexttile; plot(physical,[desPhi phi]); ylabel('phi (rad)'); legend('desired','NDI'); grid on;
nexttile; plot(physical,[desTheta theta]); ylabel('theta (rad)'); legend('desired','NDI'); grid on;
nexttile; plot(physical,1000*[qdes q]); ylabel('length change (mm)'); legend('desired l2','desired l3','NDI l2','NDI l3'); grid on;
nexttile; plot(physical,p); hold on; plot(physical,3*double(limited),'k:'); ylabel('command (bar)'); xlabel('physical poll elapsed time (s)'); legend('P1','P2','P3','limited x3'); grid on;
exportgraphics(f,fullfile(outputDir,'tracking.png'),'Resolution',150); close(f);
f=figure('Visible','off','Position',[50 50 1100 700]); tiledlayout(2,1);
nexttile; plot(physical,1000*(qdes-q)); ylabel('desired - measured (mm)'); legend('l2','l3'); grid on;
nexttile; plot(physical,[tau available]); ylabel('generalized force (N)'); xlabel('physical poll elapsed time (s)'); legend('desired 2','desired 3','available 2','available 3'); grid on;
exportgraphics(f,fullfile(outputDir,'errors_and_forces.png'),'Resolution',150); close(f);
fid=fopen(fullfile(outputDir,'report.txt'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Input: %s\nSamples: %d; physical duration %.3f s; simulation duration %.3f s\n',inputFile,numel(t),physical(end),t(end)-t(1));
fprintf(fid,'Valid %.2f%%; limited %.2f%%; median/max poll interval %.5f/%.5f s\n',100*mean(valid),100*mean(limited),results.PollMedian_s,results.PollMax_s);
fprintf(fid,'Candidate holds: %d; stationary tails: %d; targets reached and maintained: %d\n',height(holds),sum(holds.StationaryTail),sum(isfinite(holds.SettlingTime_s)));
fprintf(fid,'Static residual = nominal required force at mean measured pose minus logged available force. Interpret only stationary tails.\n');
fprintf(fid,'Available force is a model estimate, not measured chamber force. No post-Kill command or measured chamber pressure is assumed.\n');
fprintf(fid,'Transient peak tracking error is not signed overshoot; changing directions and multiple coordinates make a single overshoot percentage ambiguous.\n');
fprintf(fid,'Near-straight theta errors are excluded. Dynamic identification is not supported by static inputs.\n');
fprintf(fid,'Tail rules: last %.1f s, max speed <= %.2f mm/s, range <= %.2f mm; settling each coordinate <= %.2f mm, maintained to hold end for >= %.1f s.\n',cfg.Tail_s,1000*cfg.StationaryVelocity_m_s,1000*cfg.StationarySpread_m,1000*cfg.SettlingTolerance_m,cfg.SettlingDwell_s);
disp(holds(:,{'Hold','DesTheta_rad','DesPhi_rad','Duration_s','MeasPhi_rad','ErrorL2_mm','ErrorL3_mm','StationaryTail','LimitedFraction'}));
fprintf('Analysis saved to %s\n',outputDir);
end

function [t,a]=signal(ds,name)
names=ds.getElementNames; ix=find(strcmp(names,name));
assert(numel(ix)==1,'Expected one signal named %s, found %d.',name,numel(ix));
v=ds.getElement(ix).Values; t=double(v.Time(:)); a=double(v.Data);
if ~v.IsTimeFirst
 a=reshape(a,[],numel(t))';
else
 a=reshape(a,numel(t),[]);
end
assert(size(a,1)==numel(t));
end
function a=aligned(ds,name,t,method)
[tt,x]=signal(ds,name); [tt,ix]=unique(tt,'last'); x=x(ix,:);
assert(t(1)>=tt(1) && (strcmp(method,'previous') || t(end)<=tt(end)+.1),'Signal %s does not cover the analysis interval.',name);
% Constant parameter logs stop at the last edit; hold that value afterwards.
a=interp1(tt,x,t,method,'extrap');
end
function [elastic,gravity]=staticForce(q,p)
q=q(:); [~,~,G]=armS_single_core(q,zeros(2,1),p); K=p.K;
for i=1:2
 K(i,i)=K(i,i)+.5*p.lKbounds(3)*(2+tanh(p.mu*(q(i)-p.lKbounds(2)))-tanh(p.mu*(q(i)-p.lKbounds(1))));
end
elastic=(K*q)'; gravity=G';
end
