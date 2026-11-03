function [Q,R]=forceResidualRollout(params,B,cfg,normalizer)
%FORCERESIDUALROLLOUT Actual integrated q/v states; no measured future resets.
count=size(B.history,2);hidden=zeros(cfg.hiddenUnits,count);
for j=1:cfg.historyFrames
    x=(B.history(:,:,j)-normalizer.mean)./normalizer.std;
    [hidden,r]=forceResidualCell(params,x,hidden,cfg);
end
q=B.q0;v=B.v0;Q=cell(1,size(B.pressure,3));R=Q;
for j=1:numel(Q)
    if j>1
        values=[q;v;B.pressure(:,:,j);B.previousDt(:,:,j)];
        [hidden,r]=forceResidualCell(params,(values-normalizer.mean)./normalizer.std,hidden,cfg);
    end
    if isa(r,'dlarray')
        [q,v]=forceResidualTrainStep(q,v,B.pressure(:,:,j),r,B.dt(:,:,j),cfg);
    else
        [q,v]=forceResidualStep(q,v,B.pressure(:,:,j),r,B.dt(:,:,j),cfg);
    end
    Q{j}=q;R{j}=r;
end
Q=cat(3,Q{:});R=cat(3,R{:});
end
