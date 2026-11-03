function dX = armS_single_dynamics(t,X,params)
% One-section dynamics with full reduced-coordinate linear stiffness.
% Equal nonlinear actuator bounds include l1=-l2-l3 through virtual work.
% Optional params.boundMode=0 reproduces the legacy coordinate-only model.
q=X(1:2); dq=X(3:4);
[M,C,G]=armS_single_core(q,dq,params);
boundMode=1;if isfield(params,'boundMode'),boundMode=params.boundMode;end
bound=armS_single_boundForce(q,params.mu,params.lKbounds,boundMode);
dX=[dq;M\(params.tau-(C+params.D)*dq-G-params.K*q-bound)];
end
