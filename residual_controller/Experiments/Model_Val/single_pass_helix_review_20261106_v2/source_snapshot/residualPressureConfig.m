function cfg=residualPressureConfig(pressure,supervisor,variableCommon)
%RESIDUALPRESSURECONFIG Same vector contract as PAC tuning, no feedback gains.
%#codegen
assert(all(isfinite(pressure))&&all(isfinite(supervisor(1:6))));
assert(pressure(4)>=0&&pressure(5)>pressure(4)&&pressure(9)>=0&&pressure(9)<pressure(5));
assert(pressure(1)>=max(pressure(4),pressure(9))&&pressure(1)<=pressure(5));
assert(pressure(6)>abs(pressure(2))&&pressure(7)>abs(pressure(3))&&pressure(8)>0&&pressure(10)>0);
assert(all(supervisor(1:6)>=0)&&supervisor(2)>=supervisor(1)&&supervisor(3)>0);
assert(~isnan(supervisor(7))&&supervisor(7)>supervisor(2));
cfg=pacPressureConfig();
cfg.common_bar=pressure(1);cfg.sTrim_bar=pressure(2);cfg.aTrim_bar=pressure(3);
cfg.minimum_bar=pressure(4);cfg.maximum_bar=pressure(5);cfg.maxDeltaS_bar=pressure(6);
cfg.maxDeltaA_bar=pressure(7);cfg.slew_barPerS=pressure(8);cfg.activation_bar=pressure(9);cfg.forcePerBar_N=pressure(10);
cfg.pressureStart_s=supervisor(1);cfg.feedbackStart_s=supervisor(2);cfg.blendTime_s=supervisor(3);
cfg.invalidGrace_s=supervisor(4);cfg.maxSampleAge_s=supervisor(5);cfg.maxAnchorSpeed_mPerS=supervisor(6);cfg.stopTime_s=supervisor(7);
cfg.enableCommonAdaptation=logical(variableCommon);
trim=[cfg.common_bar-2*cfg.sTrim_bar/3;cfg.common_bar+cfg.sTrim_bar/3-cfg.aTrim_bar;cfg.common_bar+cfg.sTrim_bar/3+cfg.aTrim_bar];
assert(all(trim>=cfg.minimum_bar&trim<=cfg.maximum_bar));
end
