function next=residualFrozenStep(q,v,p,dt)
%RESIDUALFROZENSTEP Exact local affine flow with M/C/K/G frozen at current q/v.
% Offline nominal baseline, ZOH command; not exact nonlinear integration.
params=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81]);
[M,C,G]=armS_single_core(q,v,params);K=650*[2 1;1 2];
for i=1:2,K(i,i)=K(i,i)+.5e6*(2+tanh(2000*(q(i)-.035))-tanh(2000*(q(i)+.035)));end
A=[zeros(2) eye(2);-M\K -M\(C+40*eye(2))];
b=[zeros(2,1);M\(pi*.0065^2*1e5*[-1 1 0;-1 0 1]*min(max(p,0),3)-G)];
z=expm([A b;zeros(1,5)]*dt)*[q;v;1];next=z(1:4);
end
