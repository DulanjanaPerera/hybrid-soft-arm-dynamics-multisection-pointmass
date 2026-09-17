function dX = armS_hybrid_entry(t,X,L,r,mi,g,K,D,tau,mu,lKbounds,addedMass,addedXi)
%#codegen
% Runtime interface: distributed mi, then additional physical point masses.
params.N=3; params.L=L; params.r=r; params.mi=mi; params.g=g;
params.K=K; params.D=D; params.tau=tau; params.mu=mu;
params.lKbounds=lKbounds; params.addedMass=addedMass;
params.addedXi=addedXi;
dX=armS_hybrid_dynamics(t,X,params);
end
