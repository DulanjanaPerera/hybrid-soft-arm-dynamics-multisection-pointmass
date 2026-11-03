function sgFreezeRelease(root)
%SGFREEZERELEASE Research checkpoint chosen on older data, BEFORE external replay.
src=fullfile(root,'residual_dynamics_learning','force_model','streaming');camp=fullfile(src,'campaign_20261103_v2');
assert(isfile(fullfile(camp,'FINALISTS_FROZEN.mat')));
folder=fullfile(src,'release_20261103_v1');
if isfolder(folder)
    listing=dir(folder);assert(all(ismember({listing.name},{'.','..'})),'Never overwrite a release.');
else,mkdir(folder);end
F=sgRead(fullfile(camp,'finalists.csv'),'TextType','string');G=groupsummary(F,'id','mean','oldScore_mm');
[~,j]=min(G.mean_oldScore_mm);selected=F(F.id==G.id(j)&F.seed==113,:);
path=fullfile(selected.folder,'force_model.mat');assert(strcmp(sgHash(path),selected.sha256));
copyfile(path,fullfile(folder,'force_model.mat'));copyfile(fullfile(camp,'training_source'),fullfile(folder,'training_source'));
copyfile(fullfile(camp,'source_hashes.csv'),fullfile(folder,'training_source_hashes.csv'));
copyfile(fullfile(camp,'finalists.csv'),fullfile(folder,'finalists.csv'));
protocol=struct('researchOnly',true,'selectedId',selected.id,'representativeSeed',113, ...
    'selection','lowest mean old holdout score across three fixed seeds; representative seed113', ...
    'selectedMeanOldScore_mm',G.mean_oldScore_mm(j),'modelSHA256',selected.sha256,'created',char(datetime('now')), ...
    'externalLoaded',false,'previousReleaseUnchanged',true,'hardwareQualified',false);
save(fullfile(folder,'RELEASE_FROZEN.mat'),'protocol','selected','G');
fid=fopen(fullfile(folder,'FROZEN.md'),'w');fprintf(fid,'Research checkpoint selected before November2 replay.\nConfiguration %s; fixed representative seed113.\nMean old holdout score %.9g mm.\nSHA256 %s\nNo hardware/controller/inverse qualification; previous release untouched.\n',selected.id,G.mean_oldScore_mm(j),selected.sha256);fclose(fid);
disp(protocol);
end
