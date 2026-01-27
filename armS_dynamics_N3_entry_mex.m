function dX = armS_dynamics_N3_entry_mex(t, X, L, r, cog_xi, mi, g, Kmin, D, tau, mu, lKbounds)
% The function to compute the EoM matrices.
% 
% The model of ******** 3-section *********
% 
% continuum arm is computed. The joint-space varibles (length changes) are
% considered for the generalized coordinates.
% 
% Inputs:
%   t : time
%   X : length and length velocity [4Nx1]
%       [l2, l3, ..., dl2, dl3, ...]
%   params : parameters for matrices (structure)


N = length(cog_xi);

K = zeros(size(Kmin));

lmin = lKbounds(1);
lmax = lKbounds(2);
Kmax = lKbounds(3);

l = reshape(X(1:2*N,1)',[2,N])';
dl = reshape(X(2*N+1:end,1)',[2,N])';

flat_l = X(1:2*N,1);
for i=1:2*N
    K(i,i) = Kmin(i,i) + 0.5 * Kmax * (2 +tanh(mu * (flat_l(i) - lmax)) ...
        - tanh(mu * (flat_l(i)-lmin)));
end

[M, C, G] = armS_core_N3_mex(t, l, dl, L, r, cog_xi, mi, g, K);

% flat_dl = reshape(dl(1:N,:)',[2*N,1]);
Ge = K * flat_l;
flat_dl = X(2*N+1:end,1);
eps_reg = 0;
ddl = (M + eps_reg*eye(size(M))) \ (tau - (C + D)*flat_dl - (G+Ge));


% persistent minRc
% if isempty(minRc), minRc = inf; end
% rc = rcond(M);
% minRc = min(minRc, rc);
% if rc < 1e-10
%     fprintf("t=%g  rcond(M)=%e  min=%e\n", t, rc, minRc);
% end


dX = [flat_dl; ddl];

end