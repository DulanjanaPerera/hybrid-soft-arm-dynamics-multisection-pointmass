function [l, dl, ddl] = pacMeasuredDerivativesOriginal(lengths3_m, orientationValid, pollStart_s)
%#codegen
% NDI lengths are [l1;l2;l3]; model coordinates are [l2;l3].
% Invalid or repeated BX frames do not update the derivative states.
persistent lastL lastT velocity acceleration
if isempty(lastT)
    lastL = nan(2,1);
    lastT = -1;
    velocity = zeros(2,1);
    acceleration = zeros(2,1);
end
l = lastL;
dl = velocity;
ddl = acceleration;
if ~orientationValid || ~all(isfinite(lengths3_m)) || ~isfinite(pollStart_s)
    return;
end
q = [lengths3_m(2); lengths3_m(3)];
if lastT < 0
    lastL = q;
    lastT = pollStart_s;
    l = q;
    return;
end
dt = pollStart_s-lastT;
if dt > 1e-6
    rawVelocity = (q-lastL)/dt;
    oldVelocity = velocity;
    alphaV = dt/(0.20+dt);
    velocity = velocity + alphaV*(rawVelocity-velocity);
    rawAcceleration = (velocity-oldVelocity)/dt;
    alphaA = dt/(0.30+dt);
    acceleration = acceleration + alphaA*(rawAcceleration-acceleration);
    lastL = q;
    lastT = pollStart_s;
end
l = lastL;
dl = velocity;
ddl = acceleration;
end
