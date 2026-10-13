function [target, equivalentForce, residual, limited] = ...
    pacAllocatePressure(requested, anchorForce, forceBias, cfg)
%PACALLOCATEPRESSURE Fixed-mean local allocation with calibrated force offset.
% forceBias=estimated trim pressure force - nominal static anchor force.
% Force-equivalent delivered estimate is pressure force minus this bias.
% Uses a nominal max(p-activation,0) map; this is an estimate, not a sensor.
%#codegen
delta=requested-anchorForce;
dsWanted=(delta(1)+delta(2))/(2*cfg.forcePerBar_N);
daWanted=(delta(2)-delta(1))/(2*cfg.forcePerBar_N);
ds=min(max(dsWanted,-cfg.maxDeltaS_bar),cfg.maxDeltaS_bar);
da=min(max(daWanted,-cfg.maxDeltaA_bar),cfg.maxDeltaA_bar);
s=cfg.sTrim_bar+ds; a=cfg.aTrim_bar+da;
raw=[cfg.common_bar-2*s/3;cfg.common_bar+s/3-a;cfg.common_bar+s/3+a];
target=min(max(raw,cfg.minimum_bar),cfg.maximum_bar);
limited=any(abs(raw-target)>1e-10) || abs(ds-dsWanted)>1e-10 || abs(da-daWanted)>1e-10;
effective=max(target-cfg.activation_bar,0);
pressureForce=cfg.forcePerBar_N*[effective(2)-effective(1);effective(3)-effective(1)];
equivalentForce=pressureForce-forceBias;
residual=equivalentForce-requested;
end
