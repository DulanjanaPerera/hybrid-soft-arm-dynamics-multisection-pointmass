function ftEvaluate(root)
%FTEVALUATE Read Nov2 only AFTER every continuation arm is frozen.
base=fullfile(root,'residual_dynamics_learning');force=fullfile(base,'force_model');src=fullfile(force,'fine_tuning');
addpath(base,force,src,fullfile(force,'streaming'));
camp=fullfile(src,'campaign_20261103_v1');freeze=load(fullfile(camp,'FROZEN_BEFORE_EVALUATION.mat'));
results=freeze.results;cfg=freeze.cfg;D=load(fullfile(camp,'data','data.mat'));
folder=fullfile(camp,'evaluation');assert(~isfolder(folder));mkdir(folder);
assert(strcmp(sgHash(fullfile(camp,'release','force_model.mat')),freeze.protocol.selectedSHA256));
D2=sgExternalData(root,cfg);source=D2.sources;oneStep=table;short=table;quality=table;
protocol=struct('status','evaluating_frozen_arms','newDataRole','Nov3 original whole-trial TRAIN/selection masks retained', ...
    'Nov2','previously inspected external evaluation only; never fit/select', ...
    'modelPaths',results.modelPath,'modelHashes',results.sha256,'selectedArm',freeze.protocol.selectedArm);
save(fullfile(folder,'protocol.mat'),'protocol');writetable(source,fullfile(folder,'nov2_sources.csv'));
for arm=1:height(results)
    assert(strcmp(sgHash(results.modelPath(arm)),results.sha256(arm)));
    a=load(results.modelPath(arm));model=a.model;assert(isequaln(model.normalizer,D.normalizer));
    [score,T,cohorts]=ftAssess(model.params,D,cfg);
    assert(abs(score-results.compositeScore_mm(arm))<1e-8&&max(abs(cohorts-[results.oldScore_mm(arm) results.nov3Score_mm(arm)]))<1e-8);
    groups=T.group;T=removevars(T,'group');T.variant=repmat(results.arm(arm),height(T),1);
    T.cohort=repmat("old_internal",height(T),1);T.cohort(groups==2)="nov3_internal";short=[short;T]; %#ok<AGROW>
    D2.normalizer=model.normalizer;[~,T,detail]=sgAssess(model.params,D2,cfg);
    T.variant=repmat(results.arm(arm),height(T),1);T.cohort=repmat("nov2_external",height(T),1);short=[short;T]; %#ok<AGROW>
    save(fullfile(folder,sprintf('%s_nov2_short.mat',results.arm(arm))),'detail','-v7.3');
    for day=2:3
        if day==2,sequences=D2.sequences;else,sequences=D.sequences(D.runGroup==2);end
        for run=1:numel(sequences)
            S=sequences{run};[one,anchors,Q,R]=frozenForceMeasuredReplay(S,model);
            if day==2&&arm==1
                old=load(fullfile(force,'streaming','evaluation_20261103_v1',sprintf('B_seed113_run%d_stream.mat',run)));
                assert(isequal(anchors,old.anchors)&&max(abs(Q-old.Q),[],'all')<1e-12);
            end
            if day==2,partitions={'external'};else,partitions={'selection','train'};end
            for part=partitions
                part=part{1};take=true(numel(anchors),1);
                if day==3,take=S.split(anchors)==string(part);end
                assert(any(take));truth=S.q(anchors(take)+1,:)';
                labels={'force','nominal','constant_velocity'};
                for pred=1:3
                    if pred>1&&arm>1,continue;end
                    label=results.arm(arm);if pred>1,label=string(labels{pred});end
                    E=1000*(Q(:,take,pred)-truth);bad=any(~isfinite(Q(:,take,pred))|abs(Q(:,take,pred))>cfg.coordinateGuard_m,1);
                    rm=sqrt(mean(E.^2,2));vr=sqrt(mean(sum(E.^2,1)));if any(bad),vr=Inf;end
                    cohort="nov2_external";if day==3,cohort="nov3_"+part;end
                    oneStep=[oneStep;table(cohort,string(S.runId),label,nnz(take),nnz(bad),vr,rm(1),rm(2), ...
                        prctile(abs(E(1,:)),95),prctile(abs(E(2,:)),95), ...
                        'VariableNames',{'cohort','runId','variant','anchors','guardFailures','lengthVectorRMSE_mm', ...
                        'l2RMSE_mm','l3RMSE_mm','l2p95_mm','l3p95_mm'})]; %#ok<AGROW>
                end
                quality=[quality;table(cohort,string(S.runId),results.arm(arm),nnz(take),max(abs(R(:,take)),[],'all'), ...
                    mean(abs(R(:,take))>=.99*cfg.forceBound_N,'all'), ...
                    'VariableNames',{'cohort','runId','variant','anchors','maximumForce_N','forceBoundFraction'})]; %#ok<AGROW>
            end
            save(fullfile(folder,sprintf('%s_nov%d_run%d_stream.mat',results.arm(arm),day,run)), ...
                'one','anchors','Q','R','-v7.3');
        end
    end
    writetable(oneStep,fullfile(folder,'one_step.csv'));writetable(short,fullfile(folder,'short_paths.csv'));
    writetable(quality,fullfile(folder,'force_quality.csv'));
    fprintf('FROZEN EVALUATION complete: %s; Nov3 train/holdout separated, Nov2 evaluation only.\n',results.arm(arm));
end
summary=groupsummary(oneStep,{'cohort','variant'},'mean',{'lengthVectorRMSE_mm','l2RMSE_mm','l3RMSE_mm'});
writetable(summary,fullfile(folder,'one_step_summary.csv'));
for arm=1:height(results),assert(strcmp(sgHash(results.modelPath(arm)),results.sha256(arm)));end
for j=1:height(source),assert(strcmp(sgHash(source.source(j)),source.sha256(j)));end
protocol.status='complete';save(fullfile(folder,'protocol.mat'),'protocol');
copyfile(mfilename('fullpath')+".m",fullfile(folder,'ftEvaluate.m'));disp(summary);
end
