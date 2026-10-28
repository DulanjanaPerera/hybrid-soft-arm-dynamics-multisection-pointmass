function baseline=pacTuneBaselineFromTrial(trial,start_s)
%PACTUNEBASELINEFROMTRIAL Extract a recorded snapshot, not a recommendation.
if nargin<2,start_s=20;end
names={'adaptationGamma','parameterScale','parameterMaxRate','parameterLower', ...
 'parameterUpper','parameterInitial','adaptationOptions','hystoryOptions', ...
 'pressureOptions','supervisorOptions','lambda_gains','damping_gains', ...
 'enableAdaptation','enableCommonPressureAdaptation','freezeHysteresisWeights'};
baseline=struct;
for f=names
    assert(isfield(trial.signals,f{1}),'pacTune:Baseline','Missing baseline field %s.',f{1});
    signal=trial.signals.(f{1});j=find(signal.t<=start_s+1e-9,1,'last');
    assert(~isempty(j),'pacTune:Baseline','No baseline value at requested time.');
    baseline.(f{1})=signal.x(j,:)';
end
end
