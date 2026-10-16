function [adaptationGamma,parameterScale,parameterMaxRate,parameterLower, ...
 parameterUpper,parameterInitial,adaptationOptions,historyOptions, ...
 pressureOptions,supervisorOptions]=pacTunableDefaults()
%PACTUNABLEDEFAULTS Workspace values for Constant blocks; no device access.
a=pacAdaptiveConfig();p=pacPressureConfig();
adaptationGamma=a.gamma;parameterScale=a.scale;parameterMaxRate=a.maxRate;
parameterLower=a.lower;parameterUpper=a.upper;parameterInitial=a.initial;
adaptationOptions=[a.deadband_mPerS;a.residualTolerance_N;a.maxUpdateInterval_s];
historyOptions=[a.historyLength_m;a.historyBeta;a.historyGamma;a.historyBound];
pressureOptions=[p.common_bar;p.sTrim_bar;p.aTrim_bar;p.minimum_bar;p.maximum_bar; ...
 p.maxDeltaS_bar;p.maxDeltaA_bar;p.slew_barPerS;p.activation_bar;p.forcePerBar_N];
supervisorOptions=[p.pressureStart_s;p.feedbackStart_s;p.blendTime_s; ...
 p.invalidGrace_s;p.maxSampleAge_s;p.maxAnchorSpeed_mPerS;p.stopTime_s];
end
