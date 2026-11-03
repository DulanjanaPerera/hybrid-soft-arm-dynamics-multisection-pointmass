function verifyForceResidualDataContract(dataFolder,releaseFolder)
%VERIFYFORCERESIDUALDATACONTRACT Independent split/normalizer/nonzero-grad audit.
D=load(fullfile(dataFolder,'training_data.mat'));a=load(fullfile(releaseFolder,'force_model.mat'));model=a.model;
assert(all(ismember(D.manifest.cohort,{'day1_root','day2_20261027'})));
assert(~any(contains(D.manifest.source,'20261102')));
rows=[];
for i=1:numel(D.sequences)
    S=D.sequences{i};take=S.mature&S.split=="train";
    rows=[rows;S.q(take,:) S.v(take,:) S.p(take,:) S.dt(take)]; %#ok<AGROW>
    for trial=unique(S.trial(S.trial>0))'
        assert(numel(unique(S.split(S.trial==trial)))==1,'Trial split leakage.');
    end
end
assert(size(rows,1)==59630);
assert(max(abs(mean(rows)'-model.normalizer.mean))<1e-12);
stds=std(rows,0,1)';stds(stds<1e-10)=1;
assert(max(abs(stds-model.normalizer.std))<1e-12);
assert(~model.externalUsedForFitOrSelection&&strcmp(model.residualUnits,'N'));
ix=D.index(D.index.partition=="selection",:);ix=ix(1:3,:);
B=forceResidualBatch(D.sequences,ix,model.cfg,4);params=model.params;
for f=fieldnames(params)',params.(f{1})=dlarray(params.(f{1}));end
[loss,g]=dlfeval(@forceResidualLoss,params,B,model.cfg,model.normalizer);
epsilon=1e-5;plus=model.params;minus=plus;
plus.outputBias(1)=plus.outputBias(1)+epsilon;minus.outputBias(1)=minus.outputBias(1)-epsilon;
fp=value(plus,B,model.cfg,model.normalizer);fm=value(minus,B,model.cfg,model.normalizer);
fd=(fp-fm)/(2*epsilon);ad=extractdata(g.outputBias(1));relative=abs(fd-ad)/max(1e-8,abs(fd));
assert(relative<1e-3,'Nonzero trained-model gradient check failed.');
Q=forceResidualRollout(model.params,B,model.cfg,model.normalizer);changed=B;changed.truth(:)=1e5;
assert(isequaln(Q,forceResidualRollout(model.params,changed,model.cfg,model.normalizer)));
report=table(size(rows,1),extractdata(loss),ad,fd,relative,'VariableNames', ...
    {'oldTrainingNormalizerRows','trainedProbeLoss','autodiff','finiteDifference','relativeGradientError'});
writetable(report,fullfile(releaseFolder,'data_contract_verification.csv'));
copyfile(mfilename('fullpath')+".m",fullfile(releaseFolder,'source_snapshot','verifyForceResidualDataContract.m'));
disp(report);fprintf('PASS: whole-trial isolation, exact old-TRAIN-only normalization, trained nonzero-network gradient and target-free propagation.\n');
end
function loss=value(p,B,cfg,n)
[Q,R]=forceResidualRollout(p,B,cfg,n);rate=(R(:,:,2:end)-R(:,:,1:end-1))./B.dt(:,:,1:end-1)/cfg.forceRateScale_Nps;
loss=mean(((Q-B.truth)/cfg.positionScale_m).^2,'all')+cfg.smoothWeight*mean(rate.^2,'all')+ ...
    cfg.magnitudeWeight*mean((R/cfg.forceBound_N).^2,'all');
end
