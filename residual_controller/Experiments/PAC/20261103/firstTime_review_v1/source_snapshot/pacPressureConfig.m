function cfg = pacPressureConfig()
%PACPRESSURECONFIG Explicit initial command-path settings; not hardware tuned.
%#codegen
cfg.common_bar=1.6;
cfg.sTrim_bar=-0.0983981693363845;
cfg.aTrim_bar=-0.0819816714043049;
cfg.minimum_bar=0.2; % Saved base-rotation model value; BELOW activation.
cfg.maximum_bar=3.0;
cfg.maxDeltaS_bar=2.5;
cfg.maxDeltaA_bar=2.5;
cfg.slew_barPerS=2.0;
cfg.activation_bar=0.8; % Nominal analysis map; not measured chamber behavior.
cfg.forcePerBar_N=pi*0.0065^2*1e5;
cfg.pressureStart_s=10;
cfg.feedbackStart_s=15;
cfg.blendTime_s=5;
cfg.stopTime_s=inf; % No implicit experiment timeout. Disable externally.
cfg.invalidGrace_s=0.10;
cfg.maxSampleAge_s=0.075;
cfg.maxAnchorSpeed_mPerS=0.0005;
end
