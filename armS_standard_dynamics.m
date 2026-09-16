function dX = armS_standard_dynamics(t,X,params)
% Same state ordering, damping, input and stiffness law as point-mass model.
q=X(1:6); dq=X(7:12);
[M,C,G]=armS_standard_core(q,dq,params);
K=zeros(6);
for i=1:6
    K(i,i)=params.K(i,i)+0.5*params.lKbounds(3)*(2 ...
        +tanh(params.mu*(q(i)-params.lKbounds(2))) ...
        -tanh(params.mu*(q(i)-params.lKbounds(1))));
end
dX=[dq; M\(params.tau-(C+params.D)*dq-G-K*q)];
end
