function [ddq,pInverse,dInverse,limited]=residualNominal(q,dq,acc,p,boundMode)
%RESIDUALNOMINAL Fixed collection equations; optional corrected bound model.
% Default0 preserves two_day_v2 labels and historical trained models.
% Explicit mode1 uses equal three-actuator bounds for a NEW baseline version.
if nargin<5,boundMode=0;end
params=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81], ...
    'K',650*[2 1;1 2],'D',40*eye(2),'mu',2000,'lKbounds',[-.035;.035;1e6]);
[M,C,G]=armS_single_core(q,dq,params);
assert(all(isfinite(M(:)))&&all(isfinite(C(:)))&&all(isfinite(G)));
bound=armS_single_boundForce(q,params.mu,params.lKbounds,boundMode);
alpha=pi*.0065^2*1e5;A=[-1 1 0;-1 0 1];
h=(C+params.D)*dq+G+params.K*q+bound;
ddq=M\(alpha*A*min(max(p,0),3)-h);
dInverse=(M*acc+h)/alpha;
p1=max([0;-dInverse]);pInverse=[p1;p1+dInverse];
% Preserve UNCLIPPED inverse for label definition; expose feasibility flag.
limited=any(pInverse>3+1e-10)|any(pInverse< -1e-10);
end
