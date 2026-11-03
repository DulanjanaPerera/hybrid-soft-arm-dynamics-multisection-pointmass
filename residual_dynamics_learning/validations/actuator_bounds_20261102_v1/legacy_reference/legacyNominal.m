function [ddq,pInverse,dInverse,limited]=legacyNominal(q,dq,acc,p)
%RESIDUALNOMINAL Fixed collection-model equations; no actuator-memory fit.
params=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81], ...
    'K',650*[2 1;1 2],'D',40*eye(2),'mu',2000,'lKbounds',[-.035;.035;1e6]);
[M,C,G]=armS_single_core(q,dq,params);
assert(all(isfinite(M(:)))&&all(isfinite(C(:)))&&all(isfinite(G)));
K=params.K;
for j=1:2
    K(j,j)=K(j,j)+.5e6*(2+tanh(2000*(q(j)-.035))-tanh(2000*(q(j)+.035)));
end
alpha=pi*.0065^2*1e5;A=[-1 1 0;-1 0 1];
h=(C+params.D)*dq+G+K*q;
ddq=M\(alpha*A*min(max(p,0),3)-h);
dInverse=(M*acc+h)/alpha;
p1=max([0;-dInverse]);pInverse=[p1;p1+dInverse];
% Preserve UNCLIPPED inverse for label definition; expose feasibility flag.
limited=any(pInverse>3+1e-10)|any(pInverse< -1e-10);
end
