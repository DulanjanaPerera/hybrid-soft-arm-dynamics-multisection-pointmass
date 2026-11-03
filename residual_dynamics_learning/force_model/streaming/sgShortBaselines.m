function sgShortBaselines(D,detail,cfg,folder)
%SGSHORTBASELINES Matched nominal/CV sanity references, no learned fit.
B=detail.batch;N=size(B.q0,2);K=size(B.auxP,3);Q=zeros(2,N,K,2);q=B.q0;v=B.v0;
for k=1:K
    [q,v]=sgFlow(q,v,B.auxP(:,:,k),zeros(2,N),B.auxDt(:,:,k),cfg);Q(:,:,k,1)=q;
end
Q(:,:,:,2)=B.q0+B.v0.*cumsum(B.auxDt,3);T=table;labels={'nominal','constant_velocity'};
for run=unique(detail.index.run)'
    for h=[.1 .25 .5 1]
        for m=1:2
            errors=[];failed=0;
            for i=find(detail.index.run==run)'
                n=nnz(cumsum(B.auxDt(1,i,:))<=h+1e-10&B.auxMask(1,i,:)>0);
                z=Q(:,i,1:n,m);e=1000*(z-B.auxTruth(:,i,1:n));errors=[errors reshape(e,2,[])]; %#ok<AGROW>
                failed=failed+any(~isfinite(z)|abs(z)>cfg.coordinateGuard_m,'all');
            end
            rm=sqrt(mean(errors.^2,2));vr=sqrt(mean(sum(errors.^2,1)));if failed,vr=Inf;end
            T=[T;table(string(D.sequences{run}.runId),string(labels{m}),h,vr,rm(1),rm(2),failed, ...
                'VariableNames',{'runId','predictor','horizon_s','pathVectorRMSE_mm','q2RMSE_mm','q3RMSE_mm','failures'})]; %#ok<AGROW>
        end
    end
end
writetable(T,fullfile(folder,'short_baselines.csv'));save(fullfile(folder,'short_baseline_predictions.mat'),'Q','B','-v7.3');
end
