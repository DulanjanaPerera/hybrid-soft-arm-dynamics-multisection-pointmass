function validatePacTuneEvaluator()
%VALIDATEPACTUNEEVALUATOR Offline regression and failure/masking checks.
here=fileparts(mfilename('fullpath'));addpath(here,fileparts(here));
data=fullfile(fileparts(here),'experimentations','adaptCommonP');
% Alignment never crosses a gap, extrapolates, or fills an invalid native value.
s=struct('t',[0;.05;.5],'x',[1;NaN;3]);
x=pacTuneAlignSignal(s,[0;.05;.1;.2;.5;.55],.075);
assert(isequaln(x,[1;NaN;NaN;NaN;3;NaN]));
x=pacTuneAlignSignal(struct('t',0,'x',2),[0;1],.075);assert(isequaln(x,[2;NaN]));
% Synthetic Dataset tests: explicit duplicate resolution and time-last matrices.
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
ds=Simulink.SimulationData.Dataset;
for k=1:2
    e=Simulink.SimulationData.Signal;e.Name='parameterEstimate';
    e.BlockPath=Simulink.SimulationData.BlockPath(sprintf('fixture/source%d',k));
    e.Values=timeseries(reshape(1:12,[2 2 3]),[0 .05 .1],'IsTimeFirst',false);
    ds=ds.addElement(e);
end
out=struct('logsout',ds);file=fullfile(folder,'fixture.mat');save(file,'out');
cfg=pacTuneCampaignConfig();failed=false;
try,pacTuneReadTrial(file,cfg);catch err,failed=strcmp(err.identifier,'pacTune:AmbiguousSignal');end
assert(failed,'Duplicate names must require a selector.');
cfg.selectors.parameterEstimate='fixture/source2';trial=pacTuneReadTrial(file,cfg);
assert(isequal(trial.signals.parameterEstimate.x,reshape(1:12,4,3)'));
assert(trial.signals.parameterEstimate.index==2&&strlength(trial.sha256)==64);
missing=pacTuneEvaluateTrial(trial,cfg);assert(strcmp(missing.classification,'invalid_data'));
% Real v8: compare each lap against saved audit at exactly matched boundaries.
cfg=pacTuneCampaignConfig();trial=pacTuneReadTrial(fullfile(data,'1rad_10s_gainTesting_v8.mat'),cfg);
r=pacTuneEvaluateTrial(trial,cfg);old=readmatrix(fullfile(data,'1rad_10s_gainTesting_v8_analysis','lap_metrics.csv'));
assert(abs(r.period_s-10)<1e-9&&height(r.tuningEvents)==2);
assert(max(abs(r.laps.thetaMae_rad-old(:,2)))<1e-8);
assert(max(abs(r.laps.q2Rmse_mm-old(:,4)))<1e-8);
assert(~any(r.tuningEvents.eventSignal=="supervisorOptions"),'Constant Inf must not create changes.');
assert(~r.success&&strcmp(r.classification,'protocol_deviation'));
% Remove live edits in memory only: unset targets cannot claim success.
fixed=trial;fixed.signals.adaptationGamma.x(:,:)=repmat(trial.signals.adaptationGamma.x(1,:), ...
    size(trial.signals.adaptationGamma.x,1),1);
r=pacTuneEvaluateTrial(fixed,cfg);assert(strcmp(r.classification,'unconfigured_targets')&&~r.success);
for f=fieldnames(cfg.targets)',cfg.targets.(f{1})=1e6;end
r=pacTuneEvaluateTrial(fixed,cfg);assert(r.success);
cfg.targets.velocityRmse_mPerS=0;
r=pacTuneEvaluateTrial(fixed,cfg);assert(strcmp(r.classification,'valid_infeasible')&&~r.success);
% Missing/censored measurements cannot improve acceptance by masking errors.
broken=fixed;broken.signals.meas_length.x(1000:3000,:)=NaN;
r=pacTuneEvaluateTrial(broken,cfg);assert(strcmp(r.classification,'invalid_data')&&~r.success);
% 2 s run: exact historical percentages and nominal allocator replay.
cfg=pacTuneCampaignConfig();cfg.trackingWindow_s=[20 80];cfg.scoringWindow_s=[74 80];
trial=pacTuneReadTrial(fullfile(data,'1rad_10s_gainTesting_v9.mat'),cfg);r=pacTuneEvaluateTrial(trial,cfg);
assert(abs(r.period_s-2)<1e-9);
assert(abs(r.tracking.activeFraction-65/1200)<1e-10);
assert(abs(r.tracking.limitedFraction-1071/1200)<1e-10);
assert(r.allocationReplay.slewTargetMismatch==1071&&r.allocationReplay.allocationLimited==0);
assert(r.allocationReplay.loggedFlagMismatch==0&&r.gateActiveMismatch==0);
% 4 s run: fault is separate from increasing pre-fault roughness.
cfg.trackingWindow_s=[20 140];cfg.scoringWindow_s=[36 44];
trial=pacTuneReadTrial(fullfile(data,'1rad_10s_gainTesting_v10.mat'),cfg);early=pacTuneEvaluateTrial(trial,cfg);
cfg.scoringWindow_s=[52 60];late=pacTuneEvaluateTrial(trial,cfg);
assert(abs(late.period_s-4)<1e-9&&abs(late.firstFault_s-104.85)<1e-8);
assert(strcmp(late.classification,'aborted_or_faulted')&&~late.success);
assert(abs(early.scoring.thetaMae_rad-.44945)<1e-5);
assert(abs(late.scoring.thetaMae_rad-.34942)<1e-5);
assert(abs(early.scoring.thetaRateStd_radPerS-.87627)<1e-5);
assert(abs(late.scoring.thetaRateStd_radPerS-1.16627)<1e-5);
assert(late.scoring.thetaMae_rad<early.scoring.thetaMae_rad&& ...
    late.scoring.thetaRateStd_radPerS>early.scoring.thetaRateStd_radPerS);
fprintf('PASS: parser, provenance, gaps, targets, masking, v8 laps, v9 slew/gates, v10 roughness/fault.\n');
end
