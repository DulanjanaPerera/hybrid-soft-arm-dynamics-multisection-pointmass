function cfg=sgConfig()
%SGCONFIG Frozen first streaming comparison protocol, old development only.
cfg=forceResidualConfig();cfg.schema='streaming-force-2';
cfg.widths=32;cfg.horizon_s=.5;cfg.burnIn_s=.5;cfg.commonBurnIn_s=1;
cfg.chunk_s=.5;cfg.maxHorizon_s=1;cfg.batchSize=32;
cfg.learningRate=.0003;cfg.gradientClip=1;cfg.maxUpdates=1200;
cfg.validationEvery=150;cfg.minimumUpdates=600;cfg.patience=4;
cfg.selectionAnchorsPerRun=16;cfg.primaryWeight=.5;cfg.rolloutWeight=.5;
cfg.seed=113;cfg.selectionOneStepWeight=.7;cfg.selectionRolloutWeight=.3;
end
