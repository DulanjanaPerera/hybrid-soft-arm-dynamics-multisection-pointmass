function auditNov3ForceCoverage(root,folder)
%AUDITNOV3FORCECOVERAGE Compare new inputs with actual old TRAIN rows, no fit.
base=fullfile(root,'residual_dynamics_learning');stream=fullfile(base,'force_model','streaming');addpath(stream);
D=load(fullfile(base,'force_model','dataset_20261103_v1','training_data.mat'));X=[];
for j=1:numel(D.sequences)
    S=D.sequences{j};take=S.mature&S.split=="train";X=[X;S.q(take,:) S.v(take,:) S.p(take,:) S.dt(take)]; %#ok<AGROW>
end
assert(max(abs(mean(X,1)'-D.normalizer.mean))<1e-12);
oldMin=min(X,[],1);oldMax=max(X,[],1);T=table;ids={'circleSweep','SineWave','helix'};
for run=1:3
    a=load(fullfile(base,'figures','motion_predictions_20261103_v1',[ids{run} '_prepared.mat']));S=a.S;
    b=load(fullfile(folder,[ids{run} '_one_step.mat']),'anchors');take=b.anchors;
    V=[S.q(take,:) S.v(take,:) S.p(take,:) S.dt(take)];Z=(V-D.normalizer.mean')./D.normalizer.std';
    for f=1:8
        tol=1e-10*max(1,max(abs([oldMin(f) oldMax(f)])));
        outside=V(:,f)<oldMin(f)-tol|V(:,f)>oldMax(f)+tol;
        T=[T;table(string(ids{run}),string(D.normalizer.columns{f}),size(V,1),size(X,1), ...
            oldMin(f),oldMax(f),min(V(:,f)),max(V(:,f)),mean(outside),max(abs(Z(:,f))), ...
            'VariableNames',{'runId','feature','newRows','oldTrainingRows','oldTrainMinimum','oldTrainMaximum', ...
            'newMinimum','newMaximum','fractionOutsideOldObservedRange','maximumAbsoluteTrainingZ'})]; %#ok<AGROW>
    end
end
writetable(T,fullfile(folder,'training_input_coverage.csv'));disp(T(T.fractionOutsideOldObservedRange>0,:));
end
