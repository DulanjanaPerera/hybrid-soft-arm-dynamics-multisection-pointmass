function sgExternal(root)
base=fullfile(root,'residual_dynamics_learning','force_model');src=fullfile(base,'streaming');addpath(base,src);
campaign=fullfile(src,'campaign_20261103_v2');assert(isfile(fullfile(campaign,'FINALISTS_FROZEN.mat')));
assert(isfile(fullfile(src,'release_20261103_v1','RELEASE_FROZEN.mat')));
F=sgRead(fullfile(campaign,'finalists.csv'),'TextType','string');folder=fullfile(src,'evaluation_20261103_v1');
assert(~isfolder(folder));mkdir(folder);cfg=sgConfig();D=sgExternalData(root,cfg);
protocol=struct('scope','November2 previously inspected evaluation, not blind confirmation', ...
    'usedForFitOrSelection',false,'status','evaluating','commonAnchors',true,'horizons_s',[.1 .25 .5 1], ...
    'speedBins_mps',[.005 .015],'pressure','recorded commands; not measured chamber pressure');
save(fullfile(folder,'protocol.mat'),'protocol');writetable(D.sources,fullfile(folder,'sources.csv'));
allOne=table;allShort=table;allPhase=table;quality=table;
for fit=1:height(F)
    path=fullfile(F.folder(fit),'force_model.mat');assert(strcmp(sgHash(path),F.sha256(fit)));
    a=load(path);m=a.model;cfg=m.cfg;D.normalizer=m.normalizer;
    [~,short,detail]=sgAssess(m.params,D,cfg);
    if fit==1,sgShortBaselines(D,detail,cfg,folder);end
    short.id=repmat(F.id(fit),height(short),1);short.seed=repmat(F.seed(fit),height(short),1);
    allShort=[allShort;short];save(fullfile(folder,sprintf('%s_seed%d_short.mat',F.id(fit),F.seed(fit))),'detail','-v7.3');
    for run=1:numel(D.sequences)
        S=D.sequences{run};anchors=[];R=[];segmentIds=[];elapsed=tic;
        for sid=find([D.segments.run]==run)
            seg=D.segments(sid);h=sgWarm(m.params,S,seg,seg.starts(1),cfg,m.normalizer);
            first=seg.starts(1);last=seg.starts(end);rr=zeros(2,last-first+1);
            for k=first:last
                x=[S.q(k,:) S.v(k,:) S.p(k,:) S.dt(k)]';
                [h,r]=sgCell(m.params,(x-m.normalizer.mean)./m.normalizer.std,h,cfg);rr(:,k-first+1)=r;
            end
            anchors=[anchors first:last];R=[R rr];segmentIds=[segmentIds repmat(sid,1,last-first+1)]; %#ok<AGROW>
        end
        inference=toc(elapsed);N=numel(anchors);Q=zeros(2,N,3);truth=S.q(anchors+1,:)';
        for first=1:1024:N
            take=first:min(N,first+1023);k=anchors(take);q=S.q(k,:)';v=S.v(k,:)';p=S.p(k,:)';dt=S.dt(k+1)';
            Q(:,take,1)=sgFlow(q,v,p,R(:,take),dt,cfg);
            Q(:,take,2)=sgFlow(q,v,p,zeros(2,numel(k)),dt,cfg);Q(:,take,3)=q+v.*dt;
        end
        labels={'streaming_force','nominal','constant_velocity'};speed=vecnorm(S.v(anchors,:),2,2);
        for pred=1:3
            E=1000*(Q(:,:,pred)-truth);failed=any(~isfinite(Q(:,:,pred))|abs(Q(:,:,pred))>cfg.coordinateGuard_m,1);
            row=metric(E,failed,string(S.runId),string(labels{pred}));row.id=F.id(fit);row.seed=F.seed(fit);allOne=[allOne;row]; %#ok<AGROW>
            if pred==1
                for phase=unique(S.phase(anchors))'
                    mask=S.phase(anchors)==phase;a=metric(E(:,mask),failed(mask),string(S.runId),"phase"+phase);a.id=F.id(fit);a.seed=F.seed(fit);allPhase=[allPhase;a]; %#ok<AGROW>
                end
                for bin=1:3
                    mask=(1+(speed>.005)+(speed>.015))==bin;
                    if any(mask),a=metric(E(:,mask),failed(mask),string(S.runId),"speedBin"+bin);a.id=F.id(fit);a.seed=F.seed(fit);allPhase=[allPhase;a];end %#ok<AGROW>
                end
            end
        end
        adjacent=diff(segmentIds)==0;rate=diff(R,1,2)./S.dt(anchors(2:end))';rate=rate(:,adjacent);
        quality=[quality;table(F.id(fit),F.seed(fit),string(S.runId),N,max(abs(R),[],'all'), ...
            mean(abs(R(:))>=.99*cfg.forceBound_N),prctile(abs(rate(:)),95),max(abs(rate),[],'all'), ...
            1e3*inference/N,'VariableNames',{'id','seed','runId','anchors','maxForce_N','boundFraction','rateP95_Nps','maxRate_Nps','meanGRUReplay_ms'})]; %#ok<AGROW>
        save(fullfile(folder,sprintf('%s_seed%d_run%d_stream.mat',F.id(fit),F.seed(fit),run)), ...
            'anchors','Q','R','truth','labels','segmentIds','speed','-v7.3');
    end
    writetable(allOne,fullfile(folder,'one_step.csv'));writetable(allShort,fullfile(folder,'short_paths.csv'));
    writetable(allPhase,fullfile(folder,'phase_speed.csv'));writetable(quality,fullfile(folder,'force_quality.csv'));
    assert(strcmp(sgHash(path),F.sha256(fit)));fprintf('EXTERNAL complete: %s seed%d.\n',F.id(fit),F.seed(fit));
end
one=allOne(allOne.predictor=="streaming_force",:);perSeed=groupsummary(one,{'id','seed'},'mean','vectorRMSE_mm');
summary=groupsummary(perSeed,'id',{'mean','std'},'mean_vectorRMSE_mm');writetable(summary,fullfile(folder,'one_step_seed_summary.csv'));
for j=1:height(D.sources),assert(strcmp(sgHash(D.sources.source(j)),D.sources.sha256(j)));end
protocol.status='complete';save(fullfile(folder,'protocol.mat'),'protocol');
sourceFiles={'sgExternal.m','sgExternalData.m','sgShortBaselines.m','sgAssess.m','sgStates.m','sgWarm.m','sgBatch.m','sgFree.m','sgFlow.m','sgCell.m','sgRead.m'};
code=table;mkdir(fullfile(folder,'evaluation_source'));
for j=1:numel(sourceFiles)
    file=fullfile(src,sourceFiles{j});copyfile(file,fullfile(folder,'evaluation_source',sourceFiles{j}));
    code=[code;table(string(file),string(sgHash(file)),'VariableNames',{'source','sha256'})]; %#ok<AGROW>
end
writetable(code,fullfile(folder,'evaluation_code_hashes.csv'));
disp(summary);fprintf('EXTERNAL COMPLETE; no fitting or post-test selection.\n');
end
function T=metric(E,failed,id,predictor)
rm=sqrt(mean(E.^2,2));vr=sqrt(mean(sum(E.^2,1)));if any(failed),vr=Inf;end
T=table(id,predictor,size(E,2),nnz(failed),vr,rm(1),rm(2),prctile(abs(E(1,:)),95),prctile(abs(E(2,:)),95), ...
    'VariableNames',{'runId','predictor','anchors','failures','vectorRMSE_mm','q2RMSE_mm','q3RMSE_mm','q2p95_mm','q3p95_mm'});
end
