function dX = armS_standard_entry(t,X,L,r,mi,g,K,D,tau,mu,lKbounds)
%#codegen
% Fixed three-section interface. All physical parameters are runtime inputs.
params.N = 3;
params.L = L;
params.r = r;
params.mi = mi;
params.g = g;
params.K = K;
params.D = D;
params.tau = tau;
params.mu = mu;
params.lKbounds = lKbounds;
dX = armS_standard_dynamics(t,X,params);
end
