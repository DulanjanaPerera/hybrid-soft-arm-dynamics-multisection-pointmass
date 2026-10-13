function cfg=pacAdaptiveConfig()
%PACADAPTIVECONFIG Conservative initial update settings, not hardware-tuned.
%#codegen
% theta=[k1 k2 k3 d2 d3 h1 h2 h3 b2 b3]'; k3 is shared spring.
cfg.initial=[650;650;650;40;40;0;0;0;0;0];
cfg.lower=[50;50;0;1;1;0;0;0;-8;-8];
cfg.upper=[6000;6000;6000;250;250;10;10;10;8;8];
cfg.scale=[1000;1000;1000;100;100;1;1;1;1;1];
cfg.gamma=[.2;.2;.2;.0002;.0002;.5;.5;.5;1;1];
cfg.maxRate=[200;200;200;10;10;.2;.2;.2;.2;.2];
cfg.sampleTime_s=.05; cfg.maxUpdateInterval_s=.075;
cfg.deadband_mPerS=.005; cfg.residualTolerance_N=.02;
cfg.historyLength_m=.005; cfg.historyBeta=1;cfg.historyGamma=.5;
cfg.historyBound=2;
end
