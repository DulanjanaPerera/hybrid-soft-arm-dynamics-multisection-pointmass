function dX = armS_hybrid_dynamics(t,X,params) %#ok<INUSD>
% Same six coordinates, input law, damping and stiffness as standard model.
q=X(1:6); dq=X(7:12);
[M,C,G]=armS_hybrid_core(q,dq,params);
K=zeros(6);
for i=1:6
    K(i,i)=params.K(i,i)+0.5*params.lKbounds(3)*(2 ...
        +tanh(params.mu*(q(i)-params.lKbounds(2))) ...
        -tanh(params.mu*(q(i)-params.lKbounds(1))));
end
dX=[dq; M\(params.tau-(C+params.D)*dq-G-K*q)];
end
