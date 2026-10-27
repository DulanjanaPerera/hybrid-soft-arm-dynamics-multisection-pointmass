function cfg=pacTuneCampaignConfig()
%PACTUNECAMPAIGNCONFIG Offline evaluator defaults; targets intentionally unset.
cfg.schemaVersion='pac-tune-1';
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
cfg.targets=struct('qVectorRmse_mm',NaN,'thetaMae_rad',NaN, ...
    'phiMae_rad',NaN,'velocityRmse_mPerS',NaN,'limitedFraction',NaN);
end
