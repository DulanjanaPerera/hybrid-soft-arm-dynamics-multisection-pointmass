function dX = armS_dynamics_recursive(t, X, params)
% Recursive dynamics for multisection continuum arm (2 DoF per section).
% Uses HTM_nume and LocalJacob_nume for section-local kinematics and Jacobians.
%
% State:
%   X = [q; qd], where q is stacked as [l2_1; l3_1; l2_2; l3_2; ...]
%
% Dynamics:
%   M(q) qdd + C(q,qd) qd + G(q) + D qd + K q = tau

N   = params.N;
nd  = 2*N;

q  = X(1:nd);
qd = X(nd+1:end);

% Compute M(q) and G(q) using recursion (analytic Jacobians)
[M, G] = compute_M_and_G_recursive(q, params);

% Compute C(q,qd) using finite-difference dM/dq (robust)
C = compute_C_christoffel_fd(q, qd, params);

% Damping + stiffness + actuation
D   = params.D;
K   = params.K;
tau = params.tau;

% Solve for acceleration
rhs = tau - C*qd - G - D*qd - K*q;

% Numerical safety: symmetric M and regularization if needed
M = 0.5*(M + M.');
if rcond(M) < 1e-12
    M = M + 1e-8*eye(size(M));
end

qdd = M \ rhs;

dX = [qd; qdd];
end
