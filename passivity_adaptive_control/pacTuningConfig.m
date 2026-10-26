function [cfg,pc]=pacTuningConfig(gamma,scale,maxRate,lower,upper,initial,options,history,pressure,supervisor)
%PACTUNINGCONFIG Map fixed-size external tuning vectors. See TUNABLE_BLOCK.md.
%#codegen
assert(all(isfinite([gamma;scale;maxRate;lower;upper;initial;options;history;pressure;supervisor(1:6)])));
assert(all(gamma>=0) && all(scale>0) && all(maxRate>=0));
assert(all(lower<upper) && all(initial>=lower & initial<=upper));
% Experimental signed k2 compensation; box bounds no longer guarantee K SPD.
assert(lower(1)>0 && lower(2)>=-500 && lower(3)>=0 && ...
 all(lower(4:5)>0) && all(lower(6:8)>=0));
assert(all(options>=0) && options(3)>0);
assert(history(1)>0 && history(2)>=abs(history(3)) && history(4)>0);
assert(pressure(4)>=0 && pressure(5)>pressure(4) && pressure(9)>=0 && pressure(9)<pressure(5));
assert(pressure(1)>=max(pressure(4),pressure(9)) && pressure(1)<=pressure(5));
assert(pressure(6)>abs(pressure(2)) && pressure(7)>abs(pressure(3)) && all(pressure(8:10)>0));
assert(all(supervisor(1:6)>=0) && supervisor(2)>=supervisor(1) && supervisor(3)>0);
assert(~isnan(supervisor(7)) && supervisor(7)>supervisor(2));

trim=[pressure(1)-2*pressure(2)/3;pressure(1)+pressure(2)/3-pressure(3);pressure(1)+pressure(2)/3+pressure(3)];
assert(all(trim>=pressure(4) & trim<=pressure(5)));
cfg=pacAdaptiveConfig();pc=pacPressureConfig();
cfg.gamma=gamma;cfg.scale=scale;cfg.maxRate=maxRate;
cfg.lower=lower;cfg.upper=upper;cfg.initial=initial;
cfg.deadband_mPerS=options(1);cfg.residualTolerance_N=options(2);cfg.maxUpdateInterval_s=options(3);
cfg.historyLength_m=history(1);cfg.historyBeta=history(2);
cfg.historyGamma=history(3);cfg.historyBound=history(4);
pc.common_bar=pressure(1);pc.sTrim_bar=pressure(2);pc.aTrim_bar=pressure(3);
pc.minimum_bar=pressure(4);pc.maximum_bar=pressure(5);
pc.maxDeltaS_bar=pressure(6);pc.maxDeltaA_bar=pressure(7);
pc.slew_barPerS=pressure(8);pc.activation_bar=pressure(9);pc.forcePerBar_N=pressure(10);
pc.pressureStart_s=supervisor(1);pc.feedbackStart_s=supervisor(2);pc.blendTime_s=supervisor(3);
pc.invalidGrace_s=supervisor(4);pc.maxSampleAge_s=supervisor(5);
pc.maxAnchorSpeed_mPerS=supervisor(6);pc.stopTime_s=supervisor(7);
end
