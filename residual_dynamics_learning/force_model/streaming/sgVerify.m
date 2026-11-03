function sgVerify(root)
%SGVERIFY Recompute saved evaluation metrics and audit immutable inputs/models.
base=fullfile(root,'residual_dynamics_learning','force_model');src=fullfile(base,'streaming');
camp=fullfile(src,'campaign_20261103_v2');folder=fullfile(src,'evaluation_20261103_v1');
F=sgRead(fullfile(camp,'finalists.csv'),'TextType','string');O=sgRead(fullfile(folder,'one_step.csv'),'TextType','string');
sources=sgRead(fullfile(camp,'source_hashes.csv'),'TextType','string');
for j=1:height(sources),assert(strcmp(sgHash(sources.source(j)),sources.sha256(j)));end
code=sgRead(fullfile(folder,'evaluation_code_hashes.csv'),'TextType','string');
for j=1:height(code),assert(strcmp(sgHash(code.source(j)),code.sha256(j)));end
source=sgRead(fullfile(folder,'sources.csv'),'TextType','string');
for j=1:height(source),assert(strcmp(sgHash(source.source(j)),source.sha256(j)));end
maxDifference=0;names={'streaming_force','nominal','constant_velocity'};
for fit=1:height(F)
    assert(strcmp(sgHash(fullfile(F.folder(fit),'force_model.mat')),F.sha256(fit)));
    for run=1:5
        a=load(fullfile(folder,sprintf('%s_seed%d_run%d_stream.mat',F.id(fit),F.seed(fit),run)));
        for p=1:3
            e=1000*(a.Q(:,:,p)-a.truth);value=sqrt(mean(sum(e.^2,1)));
            rows=O(O.id==F.id(fit)&O.seed==F.seed(fit)&O.predictor==names{p},:);
            saved=rows.vectorRMSE_mm(run);maxDifference=max(maxDifference,abs(value-saved));
        end
    end
end
assert(maxDifference<1e-8,'Saved metrics disagree with raw predictions.');
old=fullfile(base,'release_20261103_v1','force_model.mat');assert(strcmp(sgHash(old),'b5cfd7d0cd13b2b9067a35d810b8ae166e032e96421b6a30c55d703039af268f'));
D=load(fullfile(src,'data_20261103_v1','data.mat'));original=load(fullfile(base,'dataset_20261103_v1','training_data.mat'),'normalizer','sequences');
assert(isequaln(D.normalizer,original.normalizer)&&isequaln(D.sequences,original.sequences));
for run=1:numel(D.sequences)
    S=D.sequences{run};keys=[S.channel S.trial];
    assert(isempty(intersect(unique(keys(S.split=="train",:),'rows'),unique(keys(S.split=="selection",:),'rows'),'rows')));
end
I=load(fullfile(src,'inverse_20261103_v1','summary.mat'));assert(I.summary.statePurityPass&&I.summary.maxBoundViolation_bar<1e-7);
fid=fopen(fullfile(folder,'VERIFICATION.md'),'w');
fprintf(fid,'PASS: all training sources, old data, heldout sources and frozen finalist hashes unchanged.\nPASS: old normalizer and sequences exactly preserved; no shared channel/trial keys between train/selection per recording.\nPASS: all 90 saved one-step scores independently recomputed; max discrepancy %.12g mm.\nPASS: prior release SHA256 unchanged.\nPASS: inverse pressure bounds and one-state-commit check.\nNo hardware/Simulink test or measured force validation.\n',maxDifference);fclose(fid);
fprintf('VERIFICATION PASS; maximum score discrepancy %.12g mm.\n',maxDifference);
end
