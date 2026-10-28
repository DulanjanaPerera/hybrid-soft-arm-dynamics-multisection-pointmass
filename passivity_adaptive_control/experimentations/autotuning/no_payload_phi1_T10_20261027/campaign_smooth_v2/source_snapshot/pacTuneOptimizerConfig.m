function cfg=pacTuneOptimizerConfig()
%PACTUNEOPTIMIZERCONFIG Bounded, proposal-only surrogate optimizer settings.
cfg.schemaVersion='pac-gp-1';
cfg.gammaBounds=NaN(3,2); % Bound order gammaK, gammaD, gammaH.
cfg.baselineGamma=NaN(3,1);
cfg.seed=23;
cfg.initialDesignSize=6; % Includes the mandatory baseline.
cfg.poolSize=1500;
cfg.lengthScale=.35; % Fixed RBF length in normalized log coordinates.
cfg.noiseFloor=.05; % Standardized objective standard deviation per trial.
cfg.jitter=1e-9;
cfg.improvementOffset=.01; % Standardized expected-improvement offset.
cfg.feasibilityLengthScale=.2;
cfg.minimumDistance=.01; % Exclude seen/pending points in normalized log space.
cfg.maxObservations=30;
cfg.continueSearchAfterFeasible=false; % New campaigns may optimize beyond target satisfaction.
end
