function results = estimatePlanarStiffness(inputFile)
% Offline P1 static stiffness grid search; no model/hardware execution.
here=fileparts(mfilename('fullpath')); addpath(fileparts(here));
if nargin<1, error('Pass the full path of the saved planar test MAT-file.'); end
cfg.Tail_s=3; cfg.MinHold_s=5; cfg.MinValidFraction=.9;
cfg.Kgrid_N_m=100:10:10000; cfg.Deadzone_bar=0; cfg.MinLengthNorm_m=.0005;
cfg.MaxTailRange_m=.0005;
s=load(inputFile,'out'); ds=s.out.logsout;
[t,q3]=getSignal(ds,'lengthChnage_m'); q=q3(:,2:3);
names=ds.getElementNames;
if any(strcmp(names,'des_pressure'))
 pressure=align(ds,'des_pressure',t,'previous'); pressureSource='des_pressure';
else
 pressure=[align(ds,'P1',t,'previous'),align(ds,'P2',t,'previous'),align(ds,'P3',t,'previous')];
 pressureSource='P1/P2/P3 input commands';
end
% Commanded inputs are not measured chamber pressures.
valid=logical(align(ds,'orientationValid',t,'previous'));
physical=align(ds,'pollStart',t,'linear'); physical=physical-physical(1);
phi=align(ds,'phiOrientation',t,'linear'); theta=align(ds,'thetaOrientation',t,'linear');
assert(all(diff(physical)>0),'Poll timestamps must increase.');
assert(all(abs(pressure(:,2:3))<1e-6,'all'),'This analysis requires P2=P3=0.');
p.N=1; p.L=.17411; p.r=.013; p.mi=.1; p.g=[0;0;-9.81];
B=pi*.0065^2*1e5*[-1 1 0;-1 0 1]; A=[2 1;1 2];
starts=[1;find(any(abs(diff(pressure))>1e-6,2))+1]; ends=[starts(2:end)-1;numel(t)];
rows=struct([]); objective=[];
for j=1:numel(starts)
 a=starts(j); b=ends(j); duration=physical(b)-physical(a);
 if duration<cfg.MinHold_s,continue;end
 tail=find((1:numel(t))'>=a & (1:numel(t))'<=b & physical>=physical(b)-cfg.Tail_s);
 good=tail(valid(tail)&all(isfinite(q(tail,:)),2));
 if isempty(good),continue;end
 qm=mean(q(good,:),1)'; pm=mean(pressure(good,:),1)';
 [~,~,G]=armS_single_core(qm,zeros(2,1),p);
 % Keep known bound penalty distinct from searched nominal stiffness.
 bound=.5e6*(2+tanh(2000*(qm-.035))-tanh(2000*(qm+.035))).*qm;
 force=B*max(pm-cfg.Deadzone_bar,0); z=A*qm;
 residuals=force-G-bound-z*cfg.Kgrid_N_m;
 costs=sum(residuals.^2,1); [~,best]=min(costs);
 identifiable=norm(qm)>=cfg.MinLengthNorm_m && pm(1)>0;
 range=max(q(good,:))-min(q(good,:));
 r=struct(); r.Hold=j; r.SimStart_s=t(a); r.PhysicalStart_s=physical(a);r.Duration_s=duration;
 r.P1_bar=round(pm(1)*1e6)/1e6; r.Phi_rad=mean(phi(good));
 r.Theta_rad=atan2(mean(sin(theta(good))),mean(cos(theta(good))));
 r.L1_mm=-1000*sum(qm); r.L2_mm=1000*qm(1); r.L3_mm=1000*qm(2);
 r.ValidFraction=numel(good)/numel(tail); r.TailRange_mm=1000*max(range);
 r.StationaryTail=max(range)<=cfg.MaxTailRange_m;
 r.Identifiable=identifiable; r.Stiffness_N_m=NaN; r.ForceResidual_N=NaN;r.GridEdge=false;
 if identifiable
  r.Stiffness_N_m=cfg.Kgrid_N_m(best);r.ForceResidual_N=sqrt(costs(best));r.GridEdge=best==1||best==numel(cfg.Kgrid_N_m);
 end
 r.UseForMean=identifiable&&r.StationaryTail&&r.ValidFraction>=cfg.MinValidFraction;
 if isempty(rows),rows=r;else,rows(end+1)=r;end %#ok<AGROW>
 objective(end+1,:)=costs; %#ok<AGROW>
end
assert(~isempty(rows),'No qualifying pressure holds were found.');
holds=struct2table(rows); levels=unique(holds.P1_bar(holds.UseForMean)); summary=struct([]);
for j=1:numel(levels)
 ix=holds.UseForMean&holds.P1_bar==levels(j);
 r=struct(); r.P1_bar=levels(j);r.Holds=sum(ix);r.Phi_rad=mean(holds.Phi_rad(ix));
 r.L1_mm=mean(holds.L1_mm(ix));r.L2_mm=mean(holds.L2_mm(ix));r.L3_mm=mean(holds.L3_mm(ix));
 r.StiffnessMean_N_m=mean(holds.Stiffness_N_m(ix));r.StiffnessStd_N_m=std(holds.Stiffness_N_m(ix));
 r.MeanForceResidual_N=mean(holds.ForceResidual_N(ix));r.AnyGridEdge=any(holds.GridEdge(ix));
 if isempty(summary),summary=r;else,summary(end+1)=r;end %#ok<AGROW>
end
variation=struct2table(summary);
results.InputFile=inputFile;results.Config=cfg;results.Parameters=p;results.Holds=holds;
results.PressureSource=pressureSource;
results.Variation=variation;results.ForceResidualSquared_N2=objective;
outdir=fullfile(here,'results');if ~isfolder(outdir),mkdir(outdir);end
writetable(holds,fullfile(outdir,'hold_estimates.csv'));writetable(variation,fullfile(outdir,'stiffness_variation.csv'));
save(fullfile(outdir,'stiffness_variation.mat'),'results');
if ~isempty(summary)
 f=figure('Visible','off'); errorbar(variation.Phi_rad,variation.StiffnessMean_N_m,variation.StiffnessStd_N_m,'o');
 xlabel('mean measured phi (rad)');ylabel('effective stiffness (N/m)');grid on;
 exportgraphics(f,fullfile(outdir,'stiffness_vs_phi_expanded.png'),'Resolution',150);close(f);
end
disp(variation);
end
function [t,a]=getSignal(ds,name)
ix=find(strcmp(ds.getElementNames,name));
if numel(ix)>1
 raw=[];
 for j=ix(:)'
  e=ds.getElement(j);
  if isa(e,'Simulink.SimulationData.Signal') && e.BlockPath.getLength>0
   if endsWith(e.BlockPath.getBlock(1),'/MATLAB System'),raw(end+1)=j;end %#ok<AGROW>
  end
 end
 if numel(raw)==1,ix=raw;end
end
assert(numel(ix)==1,'Need one unambiguous logged signal named %s.',name);
element=ds.getElement(ix); if isa(element,'timeseries'),v=element;else,v=element.Values;end
t=double(v.Time(:));a=double(v.Data);
if v.IsTimeFirst,a=reshape(a,numel(t),[]);else,a=reshape(a,[],numel(t))';end
[t,ix]=unique(t,'last');a=a(ix,:);
end
function a=align(ds,name,t,method)
[tt,x]=getSignal(ds,name);assert(t(1)>=tt(1)&&t(end)<=tt(end)+.001,'Incomplete %s coverage.',name);
a=interp1(tt,x,t,method,'extrap');
end
