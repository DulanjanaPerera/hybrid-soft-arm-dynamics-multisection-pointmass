function [one,anchors,Q,R,quality]=frozenForceMeasuredReplay(S,model)
%FROZENFORCEMEASUREDREPLAY Identical streaming contract to Nov2 frozen evaluation.
% No fitting, interpolation or future measured inputs. Future1s is an eligibility
% guard only, matching the existing Nov2 plots; outputs are one-step predictions.
cfg=model.cfg;anchors=[];R=[];segmentIds=[];
for sid=unique(S.segment(S.mature))'
    rows=find(S.mature&S.segment==sid);if isempty(rows),continue;end
    first=rows(1);last=rows(end);
    start=find(S.physicalTime>=S.physicalTime(first)+cfg.commonBurnIn_s,1);
    if isempty(start)||start>=last,continue;end
    take=start:last-1;take=take(S.physicalTime(last)-S.physicalTime(take)>=1);
    if isempty(take),continue;end
    assert(all(S.valid(take))&&all(S.segment(take)==sid));
    seg=struct('starts',take,'firstMature',first);h=sgWarm(model.params,S,seg,take(1),cfg,model.normalizer);
    rr=zeros(2,numel(take));
    for j=1:numel(take)
        k=take(j);x=[S.q(k,:) S.v(k,:) S.p(k,:) S.dt(k)]';
        [h,r]=sgCell(model.params,(x-model.normalizer.mean)./model.normalizer.std,h,cfg);rr(:,j)=r;
    end
    anchors=[anchors take];R=[R rr];segmentIds=[segmentIds repmat(sid,1,numel(take))]; %#ok<AGROW>
end
N=numel(anchors);Q=zeros(2,N,3);one=NaN(numel(S.t),2,3);
for first=1:1024:N
    take=first:min(N,first+1023);k=anchors(take);q=S.q(k,:)';v=S.v(k,:)';p=S.p(k,:)';dt=S.dt(k+1)';
    Q(:,take,1)=sgFlow(q,v,p,R(:,take),dt,cfg);Q(:,take,2)=sgFlow(q,v,p,zeros(2,numel(k)),dt,cfg);Q(:,take,3)=q+v.*dt;
end
for m=1:3,one(anchors+1,:,m)=Q(:,:,m)';end
labels={'force_gru64','nominal','constant_velocity'};truth=S.q(anchors+1,:)';quality=table;
x=[S.q(anchors,:) S.v(anchors,:) S.p(anchors,:) S.dt(anchors)];z=(x-model.normalizer.mean')./model.normalizer.std';
for m=1:3
    e=1000*(Q(:,:,m)-truth);bad=any(~isfinite(Q(:,:,m))|abs(Q(:,:,m))>cfg.coordinateGuard_m,1);rm=sqrt(mean(e.^2,2));
    quality=[quality;table(string(S.runId),string(labels{m}),N,nnz(bad),sqrt(mean(sum(e.^2,1))),rm(1),rm(2), ...
        max(abs(z),[],'all'),mean(any(abs(z)>5,2)), ...
        'VariableNames',{'runId','predictor','anchors','guardFailures','lengthVectorRMSE_mm','l2RMSE_mm','l3RMSE_mm','maximumTrainingZ','fractionAnyTrainingZAbove5'})]; %#ok<AGROW>
end
assert(all(isfinite(R),'all'));
end
