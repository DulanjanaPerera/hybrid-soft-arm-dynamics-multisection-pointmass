function dX = armS_single_dynamics(t,X,params)
% Same mass, stiffness, damping, gravity, and actuation laws as standard arm.
q=X(1:2); dq=X(3:4);
[M,C,G]=armS_single_core(q,dq,params);
K=zeros(2);
for i=1:2
    K(i,i)=params.K(i,i)+0.5*params.lKbounds(3)*(2 ...
        +tanh(params.mu*(q(i)-params.lKbounds(2))) ...
        -tanh(params.mu*(q(i)-params.lKbounds(1))));
end
dX=[dq;M\(params.tau-(C+params.D)*dq-G-K*q)];
end
