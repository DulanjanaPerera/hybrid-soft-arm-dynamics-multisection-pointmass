function [nextQ,nextV]=forceResidualTrainStep(q,v,p,r,dt,cfg)
%FORCERESIDUALTRAINSTEP Exact numeric flow with explicit local solver Jacobian.
% Central differences differentiate the implemented nonlinear solver; the
% affine tangent reconnects q/v/force to the neural automatic-diff graph.
% Forward values are EXACTLY forceResidualStep. First-order gradients only.
x=[q;v;r];data=extractdata(x);N=size(data,2);
epsilon=[1e-6;1e-6;1e-5;1e-5;1e-4;1e-4];
perturbed=repmat(data,1,13);
for j=1:6
    plus=j*N+(1:N);minus=(j+6)*N+(1:N);
    perturbed(j,plus)=perturbed(j,plus)+epsilon(j);
    perturbed(j,minus)=perturbed(j,minus)-epsilon(j);
end
[a,b]=forceResidualStep(perturbed(1:2,:),perturbed(3:4,:),repmat(p,1,13), ...
    perturbed(5:6,:),repmat(reshape(dt,1,[]),1,13),cfg);
all=[a;b];base=all(:,1:N);out=base;
for j=1:6
    plus=j*N+(1:N);minus=(j+6)*N+(1:N);
    jacobian=(all(:,plus)-all(:,minus))/(2*epsilon(j));
    out=out+jacobian.*(x(j,:)-data(j,:));
end
nextQ=out(1:2,:);nextV=out(3:4,:);
end
