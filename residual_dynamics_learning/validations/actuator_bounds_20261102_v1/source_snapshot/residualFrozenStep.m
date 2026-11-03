function next=residualFrozenStep(q,v,p,dt,boundMode)
%RESIDUALFROZENSTEP Exact local affine flow with M/C/K/G frozen at current q/v.
% Offline nominal baseline, ZOH command; not exact nonlinear integration.
% Default0 preserves frozen historical screening weights/context. Mode1 is
% corrected all-actuator bounds and requires a separately versioned model fit.
if nargin<5,boundMode=0;end
params=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81]);
[M,C,G]=armS_single_core(q,v,params);K=650*[2 1;1 2];
[~,boundCoefficient]=armS_single_boundForce(q,2000,[-.035;.035;1e6],boundMode);
K=K+boundCoefficient;
A=[zeros(2) eye(2);-M\K -M\(C+40*eye(2))];
b=[zeros(2,1);M\(pi*.0065^2*1e5*[-1 1 0;-1 0 1]*min(max(p,0),3)-G)];
z=expm([A b;zeros(1,5)]*dt)*[q;v;1];next=z(1:4);
end
