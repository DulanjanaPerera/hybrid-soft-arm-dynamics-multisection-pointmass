function C = compute_C_christoffel_fd(q, qd, params)
nd = length(q);
C  = zeros(nd, nd);

% step size (tune if needed)
if isfield(params, 'fd_eps') && ~isempty(params.fd_eps)
    eps = params.fd_eps;
else
    eps = 1e-6;
end

% Precompute M at +/- perturbations to form dM/dq_h
dM = zeros(nd, nd, nd);

for h = 1:nd
    dq = zeros(nd,1);
    dq(h) = eps;

    [M_plus, ~]  = compute_M_and_G_recursive(q + dq, params);
    [M_minus, ~] = compute_M_and_G_recursive(q - dq, params);

    dM(:,:,h) = (M_plus - M_minus) / (2*eps);
end

% Christoffel construction:
% C_{jk} = 1/2 * sum_h ( dM(k,j,h) + dM(k,h,j) - dM(h,j,k) ) * qd_h
for j = 1:nd
    for k = 1:nd
        s = 0;
        for h = 1:nd
            s = s + 0.5 * ( dM(k,j,h) + dM(k,h,j) - dM(h,j,k) ) * qd(h);
        end
        C(j,k) = s;
    end
end
end
