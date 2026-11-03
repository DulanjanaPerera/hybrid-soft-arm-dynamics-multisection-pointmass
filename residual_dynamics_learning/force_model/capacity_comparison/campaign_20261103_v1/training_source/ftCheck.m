function ftCheck(D,initial,cfg,folder)
%FTCHECK Split/data/gradient checks for the new fine-tuning path.
assert(isequaln(D.normalizer,initial.normalizer));
assert(height(D.addedManifest)==5&&all(strcmp(D.addedManifest.cohort,cfg.newCohort)));
assert(all(ismember(D.manifest.cohort,{'day1_root','day2_20261027'})));
for run=1:numel(D.sequences)
    S=D.sequences{run};keys=[S.channel S.trial];
    assert(isempty(intersect(keys(S.split=="train",:),keys(S.split=="selection",:),'rows')));
end
for sid=1:numel(D.segments)
    seg=D.segments(sid);S=D.sequences{seg.run};
    for j=1:numel(seg.starts)
        k=seg.starts(j);last=seg.auxEnds(j);rows=k:last;
        assert(all(S.mature(rows))&&all(S.valid(rows))&&all(S.segment(rows)==seg.sourceSegment));
        assert(all(S.split(rows)==string(seg.partition)));
        assert(S.physicalTime(last)-S.physicalTime(k)<=cfg.maxHorizon_s+1e-9);
    end
end
train=find(strcmp({D.segments.partition},'train'));jobs=zeros(8,2);
for group=1:2
    ids=train(D.runGroup([D.segments(train).run])==group);
    for j=1:4,jobs((group-1)*4+j,:)=[ids(j) 1];end
end
B=sgBatch(D,jobs,cfg);params=initial.params;
for f=fieldnames(params)',params.(f{1})=dlarray(params.(f{1}));end
assert(isequaln(sgPlain(params),initial.params));
[loss,g,~,parts]=dlfeval(@sgLoss,params,B,sgZero(cfg,8),cfg,D.normalizer);
assert(isfinite(extractdata(loss)));
for f=fieldnames(g)',assert(all(isfinite(extractdata(g.(f{1}))),'all'));end
[updated,average,squared]=adamupdate(params,g,[],[],1,cfg.learningRate);
assert(~isequaln(sgPlain(updated),initial.params));
assert(isequaln(sgPlain(params),initial.params));
proof=struct('passed',true,'normalizerPreserved',true,'noTrialLeakage',true,'noGapCrossing',true, ...
    'exactInitialWeights',true,'finiteWarmStartGradients',true,'adamStateSerializable',isstruct(sgPlain(average))&&isstruct(sgPlain(squared)), ...
    'testLoss',extractdata(loss),'testLossParts',extractdata(parts));
save(fullfile(folder,'fine_tune_checks.mat'),'proof');disp(proof);
end
