function c=pacTuneOnlineConfig()
%PACTUNEONLINECONFIG Host-side batch supervisor; deliberately unconfigured.
% Not a MATLAB Function block. No device, model or workspace writes.
c.schemaVersion='pac-online-1';
c.optimizer=pacTuneOptimizerConfig();
c.period_s=5;c.circles=10;c.scoreCircles=3;
c.holdTimeout_s=20;c.applyTimeout_s=10;c.optimizerTimeout_s=30;
c.settle_s=5;c.maxGammaFactor=1.5;c.confirmations=2;
c.targets=[2 2 .01 .025]; % l2/l3 mm, vector velocity m/s, roughness m/s^2
c.minimumCoverage=.95;
c.statePolicy='preserve';
c.workerConnectTimeout_s=1800; % Host may wait while operator preflights/compiles.
% Caller must supply and verify bounds/baseline, complete tuning snapshot,
% filter identity and trajectory identity before using a hardware adapter.
end
