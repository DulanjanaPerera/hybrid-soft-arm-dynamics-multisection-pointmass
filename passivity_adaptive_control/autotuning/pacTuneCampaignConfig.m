function cfg=pacTuneCampaignConfig()
%PACTUNECAMPAIGNCONFIG Offline evaluator defaults; targets intentionally unset.
cfg.schemaVersion='pac-tune-2';
cfg.trackingWindow_s=[20 320]; % Explicit protocol, not inferred from filename.
cfg.scoringWindow_s=[290 320];
cfg.maxAlignmentAge_s=.075;
cfg.lengthNormalizer_mm=1; % Fixed across a campaign.
cfg.minimumCoverage=.95;
cfg.expectedPeriod_s=NaN;
cfg.expectedTuning=struct;
cfg.selectors=struct; % Field = exact Dataset index or exact block path.
% Only these declared constant inputs may extend beyond their last log time.
% Reader verifies every native value is identical; varying signals never extend.
cfg.constantInputFields={'enableAdaptation','enableCommonPressureAdaptation', ...
    'freezeHysteresisWeights'};
% Required: per-coordinate tracking, length-velocity error/roughness, limiting.
cfg.targets=struct('q2Rmse_mm',NaN,'q3Rmse_mm',NaN, ...
    'velocityRmse_mPerS',NaN,'velocityErrorRateRmse_mPerS2',NaN, ...
    'limitedFraction',NaN, ...
    'qVectorRmse_mm',NaN,'thetaMae_rad',NaN,'phiMae_rad',NaN, ...
    'thetaRateStd_radPerS',NaN); % Optional: NaN disables acceptance limit.
end
