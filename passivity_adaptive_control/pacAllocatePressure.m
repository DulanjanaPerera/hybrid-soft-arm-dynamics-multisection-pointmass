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
% Invert the same activation map used below. Effective pressures must satisfy
% e2-e1 and e3-e1 = (requested+forceBias)/forcePerBar. Retain fixed mean;
% inactive muscles can take up the remaining pressure below activation.
offset=[0;(requested+forceBias)/cfg.forcePerBar_N];
xMean=cfg.common_bar-cfg.activation_bar-sum(offset)/3;
x=max(xMean,-min(offset));
ep=max(x+offset,0);
candidate=cfg.activation_bar+ep;
inactive=ep<=1e-10;
nInactive=sum(inactive);
excess=sum(candidate)-3*cfg.common_bar;
feasible=true;
if excess>1e-9
    if nInactive>0
        candidate(inactive)=candidate(inactive)-excess/nInactive;
    else
        feasible=false;
    end
end
candidateS=(candidate(2)+candidate(3))/2-candidate(1);
candidateA=(candidate(3)-candidate(2))/2;
feasible=feasible && all(candidate>=cfg.minimum_bar-1e-10) && ...
    all(candidate<=cfg.maximum_bar+1e-10) && ...
    abs(candidateS-cfg.sTrim_bar)<=cfg.maxDeltaS_bar+1e-10 && ...
    abs(candidateA-cfg.aTrim_bar)<=cfg.maxDeltaA_bar+1e-10;
if feasible
    target=min(max(candidate,cfg.minimum_bar),cfg.maximum_bar);
    limited=false;
end
effective=max(target-cfg.activation_bar,0);
pressureForce=cfg.forcePerBar_N*[effective(2)-effective(1);effective(3)-effective(1)];
equivalentForce=pressureForce-forceBias;
residual=equivalentForce-requested;
% A rail/correction/activation-infeasible request is a constrained allocation,
% even when the legacy bounded fallback happens not to touch a pressure rail.
limited=limited || norm(residual)>1e-8;
end
