function result=ftTrain(dataFolder,initialModel,folder,cfg,mode,baseline)
%FTTRAIN Weight warm start, fresh Adam, unchanged loss; balanced replay streams.
assert(ismember(mode,{'old_only','mixed_replay'}));assert(~isfolder(folder));mkdir(folder);
D=load(fullfile(dataFolder,'data.mat'));assert(all(D.runGroup==1|D.runGroup==2));
assert(all(strcmp(D.addedManifest.cohort,cfg.newCohort)));
a=load(initialModel);initial=a.model;assert(isequaln(initial.normalizer,D.normalizer));
params=initial.params;
for f=fieldnames(params)',params.(f{1})=dlarray(params.(f{1}));end
assert(isequaln(sgPlain(params),initial.params),'Fine-tune must start at exact frozen weights.');
average=[];squared=[];history=table;rng(cfg.seed+900);total=tic;
train=find(strcmp({D.segments.partition},'train'));
pool=cell(1,2);queue=pool;cursor=[1 1];
for group=1:2
    pool{group}=train(D.runGroup([D.segments(train).run])==group);
    assert(~isempty(pool{group}));queue{group}=pool{group}(randperm(numel(pool{group})));
end
slotGroup=ones(cfg.batchSize,1);if strcmp(mode,'mixed_replay'),slotGroup(cfg.batchSize/2+1:end)=2;end
slots=zeros(cfg.batchSize,2);H=sgZero(cfg,cfg.batchSize);
best=baseline.composite;selectedUpdate=0;selectedParams=initial.params;selectedCohorts=baseline.cohorts;
protocol=struct('status','fitting','mode',mode,'initialModel',initialModel,'initialSHA256',sgHash(initialModel), ...
    'optimizer','fresh Adam; exact weight warm start','normalizerRefit',false,'Nov2Loaded',false, ...
    'selection','equal old/Nov3 cohort means, each .7 one-step + .3 half-second path', ...
    'oldGuard',cfg.oldRegressionFactor*baseline.cohorts(1),'baselineEligible',true);
save(fullfile(folder,'protocol.mat'),'cfg','protocol','baseline');lossSum=zeros(1,5);clock=tic;
for update=1:cfg.maxUpdates
    numeric=sgPlain(params);
    for i=1:cfg.batchSize
        if slots(i,1)~=0,continue;end
        group=slotGroup(i);
        if cursor(group)>numel(queue{group}),queue{group}=pool{group}(randperm(numel(pool{group})));cursor(group)=1;end
        sid=queue{group}(cursor(group));cursor(group)=cursor(group)+1;slots(i,:)=[sid 1];seg=D.segments(sid);
        assert(D.runGroup(seg.run)==group&&strcmp(seg.partition,'train'));
        h=sgWarm(numeric,D.sequences{seg.run},seg,seg.starts(1),cfg,D.normalizer);
        for layer=1:numel(H),H{layer}(:,i)=h{layer};end
    end
    B=sgBatch(D,slots,cfg);[loss,g,carry,parts]=dlfeval(@sgLoss,params,B,H,cfg,D.normalizer);
    value=extractdata(loss);assert(isfinite(value),'Nonfinite fine-tune objective.');norm2=0;
    for f=fieldnames(g)'
        gradient=extractdata(g.(f{1}));assert(all(isfinite(gradient),'all'));norm2=norm2+sum(gradient.^2,'all');
    end
    factor=min(1,cfg.gradientClip/max(1e-30,sqrt(norm2)));
    for f=fieldnames(g)',g.(f{1})=g.(f{1})*factor;end
    [params,average,squared]=adamupdate(params,g,average,squared,update,cfg.learningRate);
    for layer=1:numel(H),H{layer}=extractdata(carry{layer});end
    for i=1:cfg.batchSize
        slots(i,2)=slots(i,2)+1;
        if slots(i,2)>numel(D.segments(slots(i,1)).starts),slots(i,:)=[0 0];end
    end
    lossSum=lossSum+[value extractdata(parts)];
    if mod(update,50)==0,fprintf('%s update %d/%d, loss %.5g, %.1fs.\n',mode,update,cfg.maxUpdates,value,toc(total));end
    if mod(update,cfg.validationEvery)==0
        numeric=sgPlain(params);[score,T,cohorts]=ftAssess(numeric,D,cfg);
        eligible=isfinite(score)&&cohorts(1)<=cfg.oldRegressionFactor*baseline.cohorts(1);
        history=[history;table(update,lossSum(1)/cfg.validationEvery,lossSum(2)/cfg.validationEvery, ...
            lossSum(3)/cfg.validationEvery,score,cohorts(1),cohorts(2),eligible,toc(clock),toc(total), ...
            'VariableNames',{'updates','trainingLoss','oneStepLoss','rolloutLoss','selectionScore_mm','oldScore_mm', ...
            'nov3Score_mm','passesOldGuard','interval_s','total_s'})]; %#ok<AGROW>
        writetable(history,fullfile(folder,'training_history.csv'));writetable(T,fullfile(folder,sprintf('selection_%04d.csv',update)));
        optimizer=struct('average',sgPlain(average),'squared',sgPlain(squared),'update',update,'slots',slots, ...
            'H',{H},'queue',{queue},'cursor',cursor,'rng',rng(),'slotGroup',slotGroup); %#ok<NASGU>
        save(fullfile(folder,sprintf('checkpoint_%04d.mat',update)),'numeric','optimizer','cfg','score','cohorts','eligible','-v7.3');
        if eligible&&score<best-1e-5
            best=score;selectedUpdate=update;selectedParams=numeric;selectedCohorts=cohorts;
        end
        fprintf('VALIDATION %s %d: old %.5f, Nov3 %.5f, composite %.5f mm, eligible%d, best %.5f.\n', ...
            mode,update,cohorts(1),cohorts(2),score,eligible,best);
        lossSum=zeros(1,5);clock=tic;
    end
end
model=initial;model.params=selectedParams;model.cfg=cfg;model.selectedUpdate=selectedUpdate;
model.oldSelectionScore_mm=selectedCohorts(1);
model.fineTune=struct('mode',mode,'initialSHA256',protocol.initialSHA256,'extraUpdates',selectedUpdate, ...
    'budgetUpdates',cfg.maxUpdates,'optimizerRestart',true,'selectionComposite_mm',best, ...
    'oldSelection_mm',selectedCohorts(1),'nov3Selection_mm',selectedCohorts(2),'normalizerRefit',false);
if strcmp(mode,'mixed_replay'),model.trainingDays={'20261019','20261027','20261103'};end
model.selectionDays={'20261019','20261027','20261103'};model.externalUsedForFitOrSelection=false;
save(fullfile(folder,'force_model.mat'),'model','-v7.3');protocol.status='FROZEN_BEFORE_NOV2_EVALUATION';
protocol.selectedUpdate=selectedUpdate;protocol.modelSHA256=sgHash(fullfile(folder,'force_model.mat'));
protocol.elapsed_s=toc(total);save(fullfile(folder,'protocol.mat'),'cfg','protocol','baseline');
result=struct('mode',mode,'folder',folder,'score',best,'cohorts',selectedCohorts,'selectedUpdate',selectedUpdate,'sha256',protocol.modelSHA256);
fid=fopen(fullfile(folder,'FROZEN.md'),'w');fprintf(fid,'Fine-tune arm %s frozen; selection uses old/Nov3 only.\nUpdate %d, old %.9g, Nov3 %.9g, composite %.9g mm.\nSHA256 %s\n',mode,selectedUpdate,selectedCohorts(1),selectedCohorts(2),best,result.sha256);fclose(fid);
end
