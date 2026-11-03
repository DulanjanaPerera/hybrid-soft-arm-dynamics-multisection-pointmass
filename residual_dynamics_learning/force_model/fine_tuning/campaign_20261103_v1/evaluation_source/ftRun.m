function ftRun(root)
base=fullfile(root,'residual_dynamics_learning','force_model');src=fullfile(base,'fine_tuning');stream=fullfile(base,'streaming');
addpath(base,src,stream);folder=fullfile(src,'campaign_20261103_v1');assert(~isfolder(folder));mkdir(folder);
initialFile=fullfile(stream,'release_20261103_v1','force_model.mat');initialHash=sgHash(initialFile);
assert(strcmp(initialHash,'8150352dd47230a2a4094b240e1e886b641b8c22ebad0cbfe300f63e9f526c6a'));
a=load(initialFile);initial=a.model;cfg=ftConfig(initial);
protocol=struct('status','predeclared_before_fit','initialModel',initialFile,'initialHash',initialHash, ...
    'arms',{{'frozen','old_only','mixed_replay'}},'newDay','20261103','heldOutExternalDay','20261102', ...
    'newSplit','middle listed trial per channel, same deterministic rule as original older split', ...
    'replay','16 old +16 Nov3 streams per mixed update; old-only32 old streams', ...
    'selection','50/50 cohort means; per cohort .7 one-step + .3 half-second path, macro equal recordings', ...
    'guard','old holdout score <=110% of initial old score; update0 eligible', ...
    'budget','1200 extra updates per arm, checkpoint150, no budget extension from external scores', ...
    'normalizer','unchanged older TRAIN-only','optimizer','restart Adam from exact selected weights', ...
    'noImputation',true,'externalUsedForFitOrSelection',false,'seed',cfg.seed);
save(fullfile(folder,'PREDECLARED_PROTOCOL.mat'),'protocol','cfg');
snapshot=fullfile(folder,'training_source');mkdir(snapshot);hashes=table;
files=dir(fullfile(src,'*.m'));for j=1:numel(files),p=fullfile(files(j).folder,files(j).name);copyfile(p,fullfile(snapshot,files(j).name));hashes=[hashes;hashrow(p)];end %#ok<AGROW>
dependencies={'sgAssess','sgStates','sgWarm','sgBatch','sgLoss','sgFree','sgFlow','sgZero','sgPlain','sgHash','sgRead', ...
    'sgCell','forceResidualStep','forceResidualTrainStep','forceResidualPhysics'};
for j=1:numel(dependencies)
    p=which(dependencies{j});assert(~isempty(p));[~,name,ext]=fileparts(p);copyfile(p,fullfile(snapshot,[name ext]));hashes=[hashes;hashrow(p)]; %#ok<AGROW>
end
hashes=[hashes;hashrow(initialFile)];
data=fullfile(folder,'data');ftPrepare(root,data,cfg);hashes=[hashes;hashrow(fullfile(data,'data.mat'))];
writetable(hashes,fullfile(folder,'training_source_hashes.csv'));D=load(fullfile(data,'data.mat'));
ftCheck(D,initial,cfg,folder);
[composite,T,cohorts]=ftAssess(initial.params,D,cfg);baseline=struct('composite',composite,'cohorts',cohorts);
assert(abs(cohorts(1)-initial.oldSelectionScore_mm)<1e-8,'Old baseline/anchors changed.');
writetable(T,fullfile(folder,'baseline_internal.csv'));save(fullfile(folder,'baseline.mat'),'baseline','cfg');
fprintf('INITIAL old %.6f, Nov3 %.6f, composite %.6f mm.\n',cohorts(1),cohorts(2),composite);
fid=fopen(fullfile(folder,'PROTOCOL.md'),'w');
fprintf(fid,'# Predeclared fine-tune experiment\n\nExact GRU64 weight warm start: %s\nOld October19/27 splits and normalizer unchanged. Nov3 middle whole trial per channel reserved, no interpolation. Nov2 never used for fitting/selection.\nAdam restarted, learning rate %.7g,1200 extra updates per arm; seed113. Mixed16 old/16 new streams; old-only32 old. Same loss/history/physics.\nSelection=.5old+.5Nov3, each .7oneStep+.3halfSecond macro RMSE. Old guard<=1.10 initial old score. Initial weights remain eligible.\nInitial old/new/composite: %.9g/%.9g/%.9g mm.\nAll models frozen before Nov2 reading. No hardware/controller change.\n',initialHash,cfg.learningRate,cohorts(1),cohorts(2),composite);fclose(fid);
old=ftTrain(data,initialFile,fullfile(folder,'old_only'),cfg,'old_only',baseline);
mixed=ftTrain(data,initialFile,fullfile(folder,'mixed_replay'),cfg,'mixed_replay',baseline);
results=table(["frozen";"old_only";"mixed_replay"],[0;old.selectedUpdate;mixed.selectedUpdate], ...
    [cohorts(1);old.cohorts(1);mixed.cohorts(1)],[cohorts(2);old.cohorts(2);mixed.cohorts(2)], ...
    [composite;old.score;mixed.score],[string(initialFile);string(fullfile(old.folder,'force_model.mat'));string(fullfile(mixed.folder,'force_model.mat'))], ...
    [string(initialHash);string(old.sha256);string(mixed.sha256)], ...
    'VariableNames',{'arm','extraUpdatesSelected','oldScore_mm','nov3Score_mm','compositeScore_mm','modelPath','sha256'});
writetable(results,fullfile(folder,'internal_results.csv'));[~,winner]=min(results.compositeScore_mm);
for j=1:height(hashes),assert(strcmp(sgHash(hashes.source(j)),hashes.sha256(j)),'Fitting inputs/source changed.');end
for j=1:height(D.manifest),assert(strcmp(sgHash(D.manifest.source(j)),D.manifest.sha256(j)));end
for j=1:height(D.addedManifest)
    assert(strcmp(sgHash(D.addedManifest.source(j)),D.addedManifest.sha256(j)));
    assert(strcmp(sgHash(D.addedManifest.preparedSource(j)),D.addedManifest.preparedSHA256(j)));
end
release=fullfile(folder,'release');mkdir(release);copyfile(results.modelPath(winner),fullfile(release,'force_model.mat'));
assert(strcmp(sgHash(fullfile(release,'force_model.mat')),results.sha256(winner)));
protocol.status='FROZEN_BEFORE_NOV2_EVALUATION';protocol.selectedArm=results.arm(winner);protocol.selectedSHA256=results.sha256(winner);
protocol.frozenAtUTC=char(datetime('now','TimeZone','UTC'));
save(fullfile(folder,'FROZEN_BEFORE_EVALUATION.mat'),'protocol','cfg','results','winner','baseline');
copyfile(fullfile(data,'nov3_trial_split.csv'),fullfile(release,'nov3_trial_split.csv'));
fid=fopen(fullfile(release,'FROZEN.md'),'w');fprintf(fid,'Research checkpoint selected using old/Nov3 internal holdouts BEFORE Nov2 evaluation.\nArm %s; SHA256 %s\nNo controller/hardware qualification.\n',results.arm(winner),results.sha256(winner));fclose(fid);
disp(results(:,1:5));fprintf('FINE-TUNE COMPLETE. Selected %s. All source/data/model checks pass. Nov2 not loaded.\n',results.arm(winner));
end
function row=hashrow(path)
row=table(string(path),string(sgHash(path)),'VariableNames',{'source','sha256'});
end
