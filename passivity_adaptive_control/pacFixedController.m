function [tau, s, referenceRate, referenceAcceleration, inertiaForce, ...
    coriolisForce, dampingForce, staticForce, feedbackForce, valid] = ...
    pacFixedController(measured, desired, params, lambda, dampingGain, inputValid)
%PACFIXEDCONTROLLER Nominal reduced-length passivity controller, force only.
% measured: [l2;l3;dl2;dl3], desired: [l2;l3;dl2;dl3;ddl2;ddl3].
% tau and force outputs: N; s and referenceRate: m/s; acceleration: m/s^2.
% lambda and dampingGain are diagonal entries, respectively 1/s and N*s/m.
% Invalid input produces zero outputs and valid=false. A downstream supervisor
% must hold/vent pressure appropriately; zero force is NOT a vent command.
%#codegen
tau=zeros(2,1); s=zeros(2,1); referenceRate=zeros(2,1);
referenceAcceleration=zeros(2,1); inertiaForce=zeros(2,1);
coriolisForce=zeros(2,1); dampingForce=zeros(2,1);
staticForce=zeros(2,1); feedbackForce=zeros(2,1); valid=false;
if ~inputValid || ~all(isfinite(measured(:))) || ...
        ~all(isfinite(desired(:))) || ~all(isfinite(lambda(:))) || ...
        ~all(isfinite(dampingGain(:))) || any(lambda(:)<=0) || ...
        any(dampingGain(:)<=0) || ~all(isfinite(params.g(:))) || ...
        ~all(isfinite(params.K(:))) || ~all(isfinite(params.D(:)))
    return;
end
q=measured(1:2); q=q(:); dq=measured(3:4); dq=dq(:);
qd=desired(1:2); qd=qd(:); dqd=desired(3:4); dqd=dqd(:);
ddqd=desired(5:6); ddqd=ddqd(:);
[M,C,G]=armS_single_core(q,dq,params);
Keff=params.K;
for i=1:2
    Keff(i,i)=Keff(i,i)+0.5*params.lKbounds(3)*(2+ ...
        tanh(params.mu*(q(i)-params.lKbounds(2)))- ...
        tanh(params.mu*(q(i)-params.lKbounds(1))));
end
% Reject an unusable inertia matrix before constructing controller outputs.
if ~all(isfinite(M(:))) || ~all(isfinite(C(:))) || ...
        ~all(isfinite(G(:))) || ~all(isfinite(Keff(:))) || ...
        M(1,1)<=0 || det(M)<=0 || norm(M-M.','fro')>1e-8*max(1,norm(M,'fro'))
    return;
end
referenceRate=dqd-lambda(:).*(q-qd);
referenceAcceleration=ddqd-lambda(:).*(dq-dqd);
s=dq-referenceRate;
inertiaForce=M*referenceAcceleration;
coriolisForce=C*referenceRate;
dampingForce=params.D*dq;
staticForce=G+Keff*q;
feedbackForce=-dampingGain(:).*s;
tau=inertiaForce+coriolisForce+dampingForce+staticForce+feedbackForce;
valid=all(isfinite(tau));
if ~valid
    tau(:)=0; s(:)=0; referenceRate(:)=0; referenceAcceleration(:)=0;
    inertiaForce(:)=0; coriolisForce(:)=0; dampingForce(:)=0;
    staticForce(:)=0; feedbackForce(:)=0;
end
end
