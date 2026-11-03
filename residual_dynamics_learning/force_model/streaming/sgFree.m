function [Q,R]=sgFree(params,B,state,cfg,n)
%SGFREE Auxiliary state copy; never advances the committed measured-state path.
q=B.q0;v=B.v0;Q=cell(1,size(B.auxP,3));R=Q;
for j=1:numel(Q)
    x=[q;v;B.auxP(:,:,j);B.auxPrevDt(:,:,j)];
    [state,r]=sgCell(params,(x-n.mean)./n.std,state,cfg);
    [q,v]=sgFlow(q,v,B.auxP(:,:,j),r,B.auxDt(:,:,j),cfg);
    Q{j}=q;R{j}=r;
end
Q=cat(3,Q{:});R=cat(3,R{:});
end
