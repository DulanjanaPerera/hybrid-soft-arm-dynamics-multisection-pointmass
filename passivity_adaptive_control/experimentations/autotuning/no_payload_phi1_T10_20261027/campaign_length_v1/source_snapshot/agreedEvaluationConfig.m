function cfg=agreedEvaluationConfig()
%AGREEDEVALUATIONCONFIG User-agreed length targets for the unloaded 10 s trial.
% No hardware execution. User agreed zero pressure limiting for this campaign.
cfg=pacTuneCampaignConfig();
cfg.trackingWindow_s=[20 140];cfg.scoringWindow_s=[110 140];
cfg.expectedPeriod_s=10;
cfg.targets.q2Rmse_mm=1;
cfg.targets.q3Rmse_mm=1;
cfg.targets.velocityRmse_mPerS=1e-2;
cfg.targets.velocityErrorRateRmse_mPerS2=.025;
% Angular and vector-RMSE limits remain optional/disabled.
cfg.targets.limitedFraction=0;
end
