function C=pacTuneSmoothCampaignDraft()
%PACTUNESMOOTHCAMPAIGNDRAFT New version only; cannot alter the active campaign.
% Candidate bounds are proposals, not physically validated operating limits.
C=pacTuneManualConfig();C.id='smoothness_v2_draft';
C.gammaBounds=[.10 1.0;.0005 .008;1.0 10.0];
C.evaluator.trackingWindow_s=[20 140];C.evaluator.scoringWindow_s=[110 140];
C.evaluator.expectedPeriod_s=10;
C.evaluator.objectiveMetric='velocityErrorRateRmse_mPerS2';
C.evaluator.smoothnessNormalizer_mPerS2=.025; % Scale only, not new acceptance.
C.evaluator.motionMetricSource='fixedRawLength';
C.evaluator.analysisDerivativeTau_s=[.20 .10]; % Measured/reference; fixed across runs.
C.evaluator.targets.q2Rmse_mm=2;C.evaluator.targets.q3Rmse_mm=2; % User approved.
C.evaluator.targets.velocityRmse_mPerS=1e-2;C.evaluator.targets.limitedFraction=0;
C.evaluator.targets.velocityErrorRateRmse_mPerS2=.025; % Provisional carry-over; recheck baseline.
C.filterDefinition=struct('name','pacMeasuredDerivativesOriginal','tauV_s',.20, ...
    'tauA_s',.30,'maxGap_s',Inf,'velocityOrder',1,'discretization','backward Euler', ...
    'physicalValidation','Original logged filter; user restoration recommended after matched-gain comparison.');
% All limits/objective use the fixed offline analysis filter, not live outputs.
% Purpose remains draft; no optimizer campaign or hardware run is authorized here.
end
