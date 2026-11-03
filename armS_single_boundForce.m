function [force_N,coefficient,actuatorForce_N] = ...
    armS_single_boundForce(q,mu,lengthBounds,boundMode)
%ARMS_SINGLE_BOUNDFORCE Soft actuator-length restoring force in q=[l2;l3].
% Default mode1: all l=[-l2-l3;l2;l3] use the same bounds, in metres.
% Mode0 reproduces historical coordinate-only l2/l3 bounds for old analyses.
% coefficient*q equals force_N; coefficient is NOT the tangent stiffness.
%#codegen
if nargin<4,boundMode=1;end
assert(boundMode==0||boundMode==1);
q=q(:);l=[-q(1)-q(2);q(1);q(2)];
c=.5*lengthBounds(3)*(2+tanh(mu*(l-lengthBounds(2)))- ...
    tanh(mu*(l-lengthBounds(1))));
if boundMode==0,c(1)=0;end
actuatorForce_N=c.*l;
% J=dl/dq=[-1 -1;1 0;0 1]; generalized force is J'*actuatorForce.
force_N=[actuatorForce_N(2)-actuatorForce_N(1); ...
         actuatorForce_N(3)-actuatorForce_N(1)];
coefficient=[c(1)+c(2),c(1);c(1),c(1)+c(3)];
end
