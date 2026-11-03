function trainForceResidualModel(root,dataFolder,folder,overrides)
%TRAINFORCERESIDUALMODEL Two older days only; freeze BEFORE external reading.
assert(~isfolder(folder),'Use a fresh training output folder.');mkdir(folder);
D=load(fullfile(dataFolder,'training_data.mat'));cfg=D.cfg;normalizer=D.normalizer;
cfg.selectionPatience=5;
if nargin>=4
    for f=fieldnames(overrides)'
        assert(ismember(f{1},{'learningRate','epochs','gradientClip','selectionPatience'}));
        cfg.(f{1})=overrides.(f{1});
    end
end
assert(all(ismember(D.manifest.cohort,{'day1_root','day2_20261027'})));
protocol=struct('status','training','externalLoaded',false,'trainingData',dataFolder, ...
    'selection','macro run-level length rollout RMSE on older whole-trial holdouts only', ...
    'newDataPolicy','November2 evaluation only, after force_model.mat frozen', ...
    'gradient','GRU autodiff plus central-difference solver Jacobian with exact forward values', ...
    'loss','integrated length + physical residual-force-rate and magnitude penalties');
save(fullfile(folder,'protocol.mat'),'protocol','cfg');
sourceFolder=fileparts(mfilename('fullpath'));sources=dir(fullfile(sourceFolder,'*.m'));
sources=[sources;dir(fullfile(sourceFolder,'*.py'))];mkdir(fullfile(folder,'source_snapshot'));
hashes=table;
for j=1:numel(sources)
    source=fullfile(sources(j).folder,sources(j).name);copyfile(source,fullfile(folder,'source_snapshot',sources(j).name));
    hashes=[hashes;table(string(source),string(filehash(source)),'VariableNames',{'source','sha256'})]; %#ok<AGROW>
end
file=fullfile(dataFolder,'training_data.mat');
hashes=[hashes;table(string(file),string(filehash(file)),'VariableNames',{'source','sha256'})];
writetable(hashes,fullfile(folder,'training_source_hashes.csv'));writetable(D.manifest,fullfile(folder,'raw_sources.csv'));
writetable(D.index,fullfile(folder,'window_manifest.csv'));
params=forceResidualInit(cfg);average=[];squared=[];iteration=0;history=table;
index=D.index(D.index.partition=="train",:);best=Inf;stale=0;selectedEpoch=0;
trainingClock=tic;
for epoch=1:cfg.epochs
    rng(cfg.seed+epoch);order=randperm(height(index));lossSum=0;count=0;partsSum=zeros(1,3);epochClock=tic;
    for first=1:cfg.batchSize:numel(order)
        take=order(first:min(first+cfg.batchSize-1,numel(order)));
        B=forceResidualBatch(D.sequences,index(take,:),cfg,cfg.trainSteps);
        [loss,g,parts]=dlfeval(@forceResidualLoss,params,B,cfg,normalizer);
        value=extractdata(loss);assert(isfinite(value),'Nonfinite training objective; do not hide failed batches.');
        normSquared=0;
        for f=fieldnames(g)'
            v=extractdata(g.(f{1}));
            if ~all(isfinite(v),'all')
                save(fullfile(folder,'failed_batch.mat'),'B','take','epoch','iteration','cfg','params');
                error('Nonfinite gradient in %s at epoch %d update %d; failed batch retained.',f{1},epoch,iteration);
            end
            normSquared=normSquared+sum(v.^2,'all');
        end
        factor=min(1,cfg.gradientClip/max(1e-30,sqrt(normSquared)));
        if factor<1,for f=fieldnames(g)',g.(f{1})=factor*g.(f{1});end,end
        iteration=iteration+1;
        [params,average,squared]=adamupdate(params,g,average,squared,iteration,cfg.learningRate);
        lossSum=lossSum+numel(take)*value;partsSum=partsSum+numel(take)*extractdata(parts);count=count+numel(take);
    end
    numeric=plain(params);selection=select(numeric,D,cfg);score=mean(selection.vectorRMSE_mm);
    history=[history;table(epoch,iteration,lossSum/count,partsSum(1)/count,partsSum(2)/count, ...
        score,toc(epochClock),'VariableNames',{'epoch','updates','trainingObjective','lengthLoss', ...
        'forceRateLoss','oldSelectionMacroRMSE_mm','elapsedEpoch_s'})]; %#ok<AGROW>
    writetable(history,fullfile(folder,'training_history.csv'));
    writetable(selection,fullfile(folder,sprintf('selection_epoch%02d.csv',epoch)));
    if score<best-1e-4
        best=score;selectedEpoch=epoch;selectedParams=numeric;stale=0; %#ok<NASGU>
        save(fullfile(folder,sprintf('checkpoint_epoch%02d.mat',epoch)),'numeric','cfg','normalizer','selection');
    else,stale=stale+1;end
    fprintf('Force GRU epoch %d/%d: training loss %.5g, old-day selection %.4f mm, %.1fs.\n', ...
        epoch,cfg.epochs,lossSum/count,score,toc(epochClock));
    if stale>=cfg.selectionPatience,fprintf('Stopping: old-day selection patience exhausted.\n');break;end
end
assert(selectedEpoch>0&&isfinite(best),'No finite qualified training checkpoint.');
model=struct('params',selectedParams,'cfg',cfg,'normalizer',normalizer, ...
    'selectedEpoch',selectedEpoch,'selectionMacroLengthRMSE_mm',best, ...
    'residualUnits','N','residualPlacement','M*ddq+(C+D)*dq+K*q+bound+G=B*p+r_tau', ...
    'trainingDays',{{'20261019','20261027'}},'externalUsedForFitOrSelection',false);
save(fullfile(folder,'force_model.mat'),'model','-v7.3');
protocol.status='frozen_before_new_data';protocol.elapsedTraining_s=toc(trainingClock);
protocol.selectedEpoch=selectedEpoch;protocol.selectionScore_mm=best;
save(fullfile(folder,'protocol.mat'),'protocol','cfg');
for j=1:height(hashes),assert(strcmp(filehash(hashes.source(j)),hashes.sha256(j)));end
for j=1:height(D.manifest),assert(strcmp(filehash(D.manifest.source(j)),D.manifest.sha256(j)));end
fid=fopen(fullfile(folder,'FROZEN_BEFORE_NEW_DATA.md'),'w');
fprintf(fid,'Force GRU frozen after selection on older days only.\nSelected epoch %d; old-day macro length RMSE %.8g mm.\nNo November2 recording was read by this trainer.\nModel SHA256: %s\n',selectedEpoch,best,filehash(fullfile(folder,'force_model.mat')));fclose(fid);
fprintf('FROZEN: selected epoch %d, %.4f mm old-day selection; new data never loaded.\n',selectedEpoch,best);
end
function selection=select(params,D,cfg)
selection=table;ix=D.index(D.index.partition=="selection",:);
for run=unique(ix.run)'
    take=ix(ix.run==run,:);B=forceResidualBatch(D.sequences,take,cfg,cfg.selectionSteps);
    [Q,R]=forceResidualRollout(params,B,cfg,D.normalizer);
    failed=squeeze(any(any(~isfinite(Q)|abs(Q)>cfg.coordinateGuard_m,1),3));
    e=1000*(Q-B.truth);score=sqrt(mean(sum(e.^2,1),'all'));
    if any(failed),score=Inf;end
    selection=[selection;table(string(D.sequences{run}.runId),height(take),nnz(failed),score, ...
        max(abs(R),[],'all'),'VariableNames',{'runId','windows','failures','vectorRMSE_mm','maximumResidual_N'})]; %#ok<AGROW>
end
end
function p=plain(p)
for f=fieldnames(p)',p.(f{1})=extractdata(p.(f{1}));end
end
function h=filehash(path)
f=fopen(path,'rb');assert(f>=0);c=onCleanup(@()fclose(f)); %#ok<NASGU>
bytes=fread(f,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
