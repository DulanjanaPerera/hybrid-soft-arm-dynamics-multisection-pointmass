function dX = armS_single_dynamics(t,X,params)
% One-section dynamics with full reduced-coordinate linear stiffness.
% Existing nonlinear bounds remain diagonal in q2/q3.
q=X(1:2); dq=X(3:4);
[M,C,G]=armS_single_core(q,dq,params);
K=params.K; % Preserve physical coupling; add bound penalties on the diagonal.
for i=1:2
    K(i,i)=params.K(i,i)+0.5*params.lKbounds(3)*(2 ...
        +tanh(params.mu*(q(i)-params.lKbounds(2))) ...
        -tanh(params.mu*(q(i)-params.lKbounds(1))));
end
dX=[dq;M\(params.tau-(C+params.D)*dq-G-K*q)];
end
