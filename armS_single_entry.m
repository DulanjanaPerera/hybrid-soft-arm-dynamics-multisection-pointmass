function dX = armS_single_entry(t,X,L,r,mi,g,K,D,tau,mu,lKbounds)
%#codegen
% One-section, two-coordinate distributed-mass RHS.
params.N=1;
params.L=L;
params.r=r;
params.mi=mi;
params.g=g;
params.K=K;
params.D=D;
params.tau=tau;
params.mu=mu;
params.lKbounds=lKbounds;
dX=armS_single_dynamics(t,X,params);
end
