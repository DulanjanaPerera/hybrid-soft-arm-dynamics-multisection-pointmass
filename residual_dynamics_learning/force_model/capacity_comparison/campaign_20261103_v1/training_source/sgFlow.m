function [q,v]=sgFlow(q,v,p,r,dt,cfg)
if ~isa(r,'dlarray'),[q,v]=forceResidualStep(q,v,p,r,dt,cfg);return;end
if isa(q,'dlarray')||isa(v,'dlarray')
    [q,v]=forceResidualTrainStep(q,v,p,r,dt,cfg);return;
end
% Measured-state branch: only force is learnable, so only its Jacobian is needed.
data=extractdata(r);N=size(data,2);epsilon=1e-4;rr=repmat(data,1,5);
for j=1:2
    rr(j,j*N+(1:N))=rr(j,j*N+(1:N))+epsilon;
    rr(j,(j+2)*N+(1:N))=rr(j,(j+2)*N+(1:N))-epsilon;
end
[a,b]=forceResidualStep(repmat(q,1,5),repmat(v,1,5),repmat(p,1,5),rr,repmat(dt,1,5),cfg);
out=[a(:,1:N);b(:,1:N)];all=[a;b];
for j=1:2
    jac=(all(:,j*N+(1:N))-all(:,(j+2)*N+(1:N)))/(2*epsilon);
    out=out+jac.*(r(j,:)-data(j,:));
end
q=out(1:2,:);v=out(3:4,:);
end
