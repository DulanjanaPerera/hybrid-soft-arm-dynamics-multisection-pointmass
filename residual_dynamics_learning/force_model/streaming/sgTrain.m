function model=sgTrain(dataFolder,folder,cfg)
%SGTRAIN Ordered measured-state TBPTT; auxiliary hidden states never committed.
assert(~isfolder(folder),'Retain attempts; use a fresh folder.');mkdir(folder);
D=load(fullfile(dataFolder,'data.mat'));assert(all(ismember(D.manifest.cohort,{'day1_root','day2_20261027'})));
protocol=struct('status','training_old_only','externalLoaded',false,'dataFolder',dataFolder, ...
    'gradient','exact forward flow; central-difference first-order solver Jacobian', ...
    'selection','.7 macro one-step vector RMSE + .3 macro .5s path vector RMSE', ...
    'state','observed streaming TBPTT; auxiliary copy discarded');
save(fullfile(folder,'protocol.mat'),'cfg','protocol');
params=sgInit(cfg);average=[];squared=[];history=table;best=Inf;stale=0;
training=find(arrayfun(@(s)strcmp(s.partition,'train'),D.segments));
rng(cfg.seed+900);queue=training(randperm(numel(training)));cursor=1;
slots=zeros(cfg.batchSize,2);H=sgZero(cfg,cfg.batchSize);total=tic;last=tic;
lossSum=zeros(1,5);selectedUpdate=0;
for update=1:cfg.maxUpdates
    numeric=sgPlain(params);
    for i=1:cfg.batchSize
        if slots(i,1)==0
            if cursor>numel(queue),queue=training(randperm(numel(training)));cursor=1;end
            sid=queue(cursor);cursor=cursor+1;slots(i,:)=[sid 1];seg=D.segments(sid);
            h=sgWarm(numeric,D.sequences{seg.run},seg,seg.starts(1),cfg,D.normalizer);
            for l=1:numel(H),H{l}(:,i)=h{l};end
        end
    end
    B=sgBatch(D,slots,cfg);
    [loss,g,carry,parts]=dlfeval(@sgLoss,params,B,H,cfg,D.normalizer);
    value=extractdata(loss);assert(isfinite(value),'Nonfinite loss; attempt retained.');norm2=0;
    for f=fieldnames(g)'
        v=extractdata(g.(f{1}));assert(all(isfinite(v),'all'),'Nonfinite gradient.');norm2=norm2+sum(v.^2,'all');
    end
    factor=min(1,cfg.gradientClip/max(1e-30,sqrt(norm2)));
    for f=fieldnames(g)',g.(f{1})=g.(f{1})*factor;end
    [params,average,squared]=adamupdate(params,g,average,squared,update,cfg.learningRate);
    for l=1:numel(H),H{l}=extractdata(carry{l});end % detach measured stream only
    for i=1:size(slots,1)
        slots(i,2)=slots(i,2)+1;
        if slots(i,2)>numel(D.segments(slots(i,1)).starts),slots(i,:)=[0 0];end
    end
    lossSum=lossSum+[value extractdata(parts)];
    if mod(update,50)==0,fprintf('Update %d/%d, loss %.5g, total %.1fs.\n',update,cfg.maxUpdates,value,toc(total));end
    if mod(update,cfg.validationEvery)==0||update==cfg.maxUpdates
        numeric=sgPlain(params);[score,selection]=sgAssess(numeric,D,cfg);
        history=[history;table(update,lossSum(1)/cfg.validationEvery,lossSum(2)/cfg.validationEvery, ...
            lossSum(3)/cfg.validationEvery,score,toc(last),toc(total), ...
            'VariableNames',{'updates','loss','oneStepLoss','rolloutLoss','oldScore_mm','interval_s','total_s'})]; %#ok<AGROW>
        writetable(history,fullfile(folder,'training_history.csv'));
        writetable(selection,fullfile(folder,sprintf('selection_%04d.csv',update)));
        save(fullfile(folder,sprintf('checkpoint_%04d.mat',update)),'numeric','cfg','selection','score','-v7.3');
        if score<best-1e-5,best=score;selectedParams=numeric;selectedUpdate=update;stale=0;else,stale=stale+1;end
        fprintf('SELECTION update %d: %.5f mm, best %.5f, elapsed %.1fs.\n',update,score,best,toc(total));
        lossSum=zeros(1,5);last=tic;
        if update>=cfg.minimumUpdates&&stale>=cfg.patience,break;end
    end
end
assert(selectedUpdate>0&&isfinite(best));
model=struct('params',selectedParams,'cfg',cfg,'normalizer',D.normalizer,'selectedUpdate',selectedUpdate, ...
    'oldSelectionScore_mm',best,'trainingDays',{{'20261019','20261027'}},'externalUsedForFitOrSelection',false, ...
    'residualUnits','N','residualPlacement','M*ddq+(C+D)*dq+K*q+bound+G=B*p+r_tau');
save(fullfile(folder,'force_model.mat'),'model','-v7.3');
protocol.status='frozen_before_November2';protocol.modelSHA256=sgHash(fullfile(folder,'force_model.mat'));
protocol.elapsed_s=toc(total);save(fullfile(folder,'protocol.mat'),'cfg','protocol');
fid=fopen(fullfile(folder,'FROZEN.md'),'w');fprintf(fid,'Old-only checkpoint selected at update %d; score %.9g mm.\nSHA256 %s\nNo November2 fit/selection.\n',selectedUpdate,best,protocol.modelSHA256);fclose(fid);
end
