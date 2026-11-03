function cfg=ftConfig(model)
%FTCONFIG Declared before fitting; architecture/physics/loss/normalizer locked.
cfg=model.cfg;cfg.schema='force-finetune-1';cfg.learningRate=1e-4;
cfg.maxUpdates=1200;cfg.validationEvery=150;cfg.seed=113;cfg.batchSize=32;
cfg.oldRegressionFactor=1.10;cfg.cohortSelectionWeights=[.5 .5];
cfg.newCohort='day4_20261103';cfg.optimizerRestart=true;
end
