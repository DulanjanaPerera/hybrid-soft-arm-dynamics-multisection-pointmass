function releaseForceResidualModel(dataFolder,candidates,folder)
%RELEASEFORCERESIDUALMODEL Select only by old-trial score, freeze exact bytes.
assert(~isfolder(folder));comparison=table;
for j=1:numel(candidates)
    assert(isfile(fullfile(candidates{j},'FROZEN_BEFORE_NEW_DATA.md')));
    a=load(fullfile(candidates{j},'force_model.mat'));m=a.model;
    assert(~m.externalUsedForFitOrSelection);
    comparison=[comparison;table(string(candidates{j}),m.cfg.learningRate,m.selectedEpoch, ...
        m.selectionMacroLengthRMSE_mm,'VariableNames',{'trainingFolder','learningRate','epoch','oldSelectionMacroRMSE_mm'})]; %#ok<AGROW>
end
[~,j]=min(comparison.oldSelectionMacroRMSE_mm);selected=candidates{j};
for name={'forceResidualCell.m','forceResidualPhysics.m','forceResidualStep.m','forceResidualRollout.m'}
    assert(strcmp(filehash(which(name{1})),filehash(fullfile(selected,'source_snapshot',name{1}))), ...
        'Inference code differs from training snapshot.');
end
mkdir(folder);copyfile(fullfile(selected,'force_model.mat'),fullfile(folder,'force_model.mat'));
writetable(comparison,fullfile(folder,'old_only_model_comparison.csv'));
sourceModel=fullfile(selected,'force_model.mat');selectedHash=filehash(sourceModel);
assert(strcmp(filehash(fullfile(folder,'force_model.mat')),selectedHash));
copyfile(fullfile(selected,'source_snapshot'),fullfile(folder,'source_snapshot'));
copyfile(fullfile(selected,'training_history.csv'),fullfile(folder,'training_history.csv'));
copyfile(fullfile(selected,'raw_sources.csv'),fullfile(folder,'raw_sources.csv'));
copyfile(fullfile(selected,'window_manifest.csv'),fullfile(folder,'window_manifest.csv'));
a=load(fullfile(folder,'force_model.mat'));model=a.model;D=load(fullfile(dataFolder,'training_data.mat'));
cfg=model.cfg;report=table;index=D.index(D.index.partition=="selection",:);
for run=unique(index.run)'
    take=index(index.run==run,:);B=forceResidualBatch(D.sequences,take,cfg,cfg.selectionSteps);
    [learned,~]=forceResidualRollout(model.params,B,cfg,model.normalizer);
    nominal=zeros(size(learned));q=B.q0;v=B.v0;time=cumsum(B.dt,3);
    cv=B.q0+B.v0.*time;
    for step=1:cfg.selectionSteps
        [q,v]=forceResidualStep(q,v,B.pressure(:,:,step),zeros(2,height(take)),B.dt(:,:,step),cfg);
        nominal(:,:,step)=q;
    end
    preds={nominal,cv,learned};names={'nominal','constant_velocity','force_gru'};
    for j=1:3
        q=preds{j};bad=squeeze(any(any(~isfinite(q)|abs(q)>cfg.coordinateGuard_m,1),3));
        error=1000*(q-B.truth);rmse=sqrt(mean(sum(error.^2,1),'all'));
        report=[report;table(string(D.sequences{run}.runId),string(names{j}),height(take),nnz(bad),rmse, ...
            'VariableNames',{'runId','predictor','windows','failures','vectorRMSE_mm'})]; %#ok<AGROW>
    end
end
writetable(report,fullfile(folder,'old_selection_baselines.csv'));
fid=fopen(fullfile(folder,'FROZEN_BEFORE_NEW_DATA.md'),'w');
fprintf(fid,'Selected between two fits by OLD-TRIAL length error only.\nSource: %s\nExact model SHA256: %s\nNo November2 evaluation used in fitting or this selection.\n',selected,selectedHash);fclose(fid);
copyfile(mfilename('fullpath')+".m",fullfile(folder,'source_snapshot','releaseForceResidualModel.m'));
disp(comparison);disp(report);fprintf('RELEASE FROZEN: %s\n',selected);
end
function h=filehash(path)
f=fopen(path,'rb');assert(f>=0);c=onCleanup(@()fclose(f)); %#ok<NASGU>
bytes=fread(f,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
