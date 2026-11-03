function [score,T,detail]=sgAssess(params,D,cfg)
%SGASSESS Old-only ranking: .7 one-step + .3 half-second run-macro vector RMSE.
ix=D.index;H=sgStates(params,D,ix,cfg);jobs=zeros(height(ix),2);
for i=1:height(ix)
    sid=ix.segmentId(i);jobs(i,:)=[sid find(D.segments(sid).starts==ix.anchor(i),1)];
end
evalcfg=cfg;evalcfg.horizon_s=1;B=sgBatch(D,jobs,evalcfg);
[~,R]=sgCell(params,(B.x(:,:,1)-D.normalizer.mean)./D.normalizer.std,H,cfg);
[one,~]=sgFlow(B.q0,B.v0,B.pressure(:,:,1),R,B.dt(:,:,1),cfg);
[Q,forces]=sgFree(params,B,H,cfg,D.normalizer);
E=1000*(one-B.truth(:,:,1));T=table;horizons=[.1 .25 .5 1];
for run=unique(ix.run)'
    take=find(ix.run==run);S=D.sequences{run};
    for hz=horizons
        errors=[];rates=[];bounds=[];fail=0;
        for i=take'
            count=nnz(cumsum(B.auxDt(1,i,:))<=hz+1e-10 & B.auxMask(1,i,:)>0);
            if count<1,continue;end
            path=Q(:,i,1:count);err=1000*(path-B.auxTruth(:,i,1:count));
            fail=fail+any(~isfinite(path)|abs(path)>cfg.coordinateGuard_m,'all');
            errors=[errors reshape(err,2,[])]; %#ok<AGROW>
            f=reshape(forces(:,i,1:count),2,[]);bounds=[bounds reshape(abs(f)>=.99*cfg.forceBound_N,1,[])]; %#ok<AGROW>
            if count>1,rates=[rates reshape(diff(f,1,2)./reshape(B.auxDt(1,i,1:count-1),1,[]),1,[])];end %#ok<AGROW>
        end
        rm=sqrt(mean(errors.^2,2));vr=sqrt(mean(sum(errors.^2,1)));
        if fail>0,vr=Inf;end
        T=[T;table(string(S.runId),hz,numel(take),sqrt(mean(sum(E(:,take).^2,1))), ...
            vr,rm(1),rm(2),prctile(abs(errors(1,:)),95),prctile(abs(errors(2,:)),95), ...
            fail,prctile(abs(rates),95),mean(bounds), ...
            'VariableNames',{'runId','horizon_s','anchors','oneStepVectorRMSE_mm','pathVectorRMSE_mm', ...
            'q2RMSE_mm','q3RMSE_mm','q2p95_mm','q3p95_mm','failures','forceRateP95_Nps','forceBoundFraction'})]; %#ok<AGROW>
    end
end
half=T(abs(T.horizon_s-.5)<1e-8,:);
score=cfg.selectionOneStepWeight*mean(half.oneStepVectorRMSE_mm)+cfg.selectionRolloutWeight*mean(half.pathVectorRMSE_mm);
detail=struct('index',ix,'oneStepError_mm',E,'Q',Q,'R',forces,'batch',B,'preStates',{H});
end
