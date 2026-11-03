function evaluateForceResidualModel(root,trainingFolder,folder)
%EVALUATEFORCERESIDUALMODEL November2 evaluation only AFTER immutable freeze.
assert(~isfolder(folder));assert(isfile(fullfile(trainingFolder,'FROZEN_BEFORE_NEW_DATA.md')));
loaded=load(fullfile(trainingFolder,'force_model.mat'));model=loaded.model;cfg=model.cfg;
assert(~model.externalUsedForFitOrSelection&&cfg.nominal.boundMode==1);
assert(isequal(model.trainingDays,{'20261019','20261027'}));
modelPath=fullfile(trainingFolder,'force_model.mat');modelHash=filehash(modelPath);
mkdir(folder);protocol=struct('status','evaluating','frozenModel',modelPath, ...
    'frozenModelSHA256',modelHash,'fitOrSelectionOnThisData',false, ...
    'runs',{{'StepWave','SineWave4','combinePressure','circleSweep','helix'}}, ...
    'horizons_s',[.5 2 5],'windowsPerRunHorizon',12,'scope','held-out evaluation; previously inspected cohort');
save(fullfile(folder,'protocol.mat'),'protocol','cfg');
prepared=fullfile(root,'residual_dynamics_learning','validations','fresh_20261102_v2');
names={'nominal','constant_velocity','force_gru'};oneStep=table;windowScores=table;windows=table;
quality=table;provenance=table;trajectories=struct([]);
for run=1:numel(protocol.runs)
    id=protocol.runs{run};file=fullfile(prepared,[id '_sequences.mat']);a=load(file);S=a.S;
    assert(strcmp(S.cohort,'day3_20261102')&&all(S.split=="fresh_test"));
    index=table(ones(numel(S.forwardAnchors),1),S.forwardAnchors,'VariableNames',{'run','anchor'});
    N=height(index);Q=NaN(2,N,3);R=NaN(2,N);truth=S.q(index.anchor+1,:)';
    for first=1:256:N
        take=first:min(first+255,N);B=forceResidualBatch({S},index(take,:),cfg,1);
        [q,r]=forceResidualRollout(model.params,B,cfg,model.normalizer);Q(:,take,3)=q;R(:,take)=r;
        [q,~]=forceResidualStep(B.q0,B.v0,B.pressure(:,:,1),zeros(2,numel(take)),B.dt(:,:,1),cfg);
        Q(:,take,1)=q;Q(:,take,2)=B.q0+B.v0.*B.dt(:,:,1);
    end
    for m=1:3
        bad=any(~isfinite(Q(:,:,m))|abs(Q(:,:,m))>cfg.coordinateGuard_m,1);
        e=1000*(Q(:,:,m)-truth);good=~bad;
        rmse=sqrt(mean(sum(e(:,good).^2,1)));coordinates=sqrt(mean(e(:,good).^2,2));
        oneStep=[oneStep;table(string(id),string(names{m}),N,nnz(bad),rmse,coordinates(1),coordinates(2), ...
            'VariableNames',{'runId','predictor','anchors','failures','successfulVectorRMSE_mm', ...
            'successfulL2RMSE_mm','successfulL3RMSE_mm'})]; %#ok<AGROW>
    end
    save(fullfile(folder,[id '_one_step.mat']),'index','Q','R','truth','names');
    quality=[quality;table(string(id),N,max(abs(R),[],'all'),mean(any(abs(R)>.95*cfg.forceBound_N,1)), ...
        'VariableNames',{'runId','anchors','maxAbsoluteResidual_N','fractionNearForceBound'})]; %#ok<AGROW>
    raw=fullfile(root,'RL_data_collection','20261102',[id '.mat']);
    for f={raw,file}
        provenance=[provenance;table(string(id),string(f{1}),string(filehash(f{1})), ...
            'VariableNames',{'runId','source','sha256'})]; %#ok<AGROW>
    end
    for horizon=protocol.horizons_s
        eligible=[];ends=[];
        for k=S.forwardAnchors'
            last=find(S.physicalTime>=S.physicalTime(k)+horizon,1);if isempty(last),continue;end
            rows=k-cfg.historyFrames:last;
            if all(S.mature(rows))&&all(S.segment(rows)==S.segment(k)),eligible(end+1)=k;ends(end+1)=last;end %#ok<AGROW>
        end
        take=unique(round(linspace(1,numel(eligible),min(protocol.windowsPerRunHorizon,numel(eligible)))));
        starts=eligible(take);lasts=ends(take);steps=lasts-starts;
        for K=unique(steps)
            group=find(steps==K);anchors=starts(group)';last=lasts(group)';count=numel(group);
            ix=table(ones(count,1),anchors,'VariableNames',{'run','anchor'});
            B=forceResidualBatch({S},ix,cfg,K);
            [predicted,residual]=forceResidualRollout(model.params,B,cfg,model.normalizer);
            nominal=zeros(2,count,K);q=B.q0;v=B.v0;velocity=nominal;
            elapsed=cumsum(B.dt,3);
            for step=1:K
                [q,v]=forceResidualStep(q,v,B.pressure(:,:,step),zeros(2,count),B.dt(:,:,step),cfg);
                nominal(:,:,step)=q;velocity(:,:,step)=B.q0+B.v0.*elapsed(:,:,step);
            end
            predictions={nominal,velocity,predicted};
            for w=1:count
                actualDuration=S.physicalTime(last(w))-S.physicalTime(anchors(w));
                windows=[windows;table(string(id),horizon,anchors(w),last(w),S.trial(anchors(w)), ...
                    actualDuration,'VariableNames',{'runId','horizon_s','anchor','endRow','trial','actualPhysicalHorizon_s'})]; %#ok<AGROW>
                for m=1:3
                    z=reshape(predictions{m}(:,w,:),2,K);truthWindow=reshape(B.truth(:,w,:),2,K);
                    failed=any(~isfinite(z)|abs(z)>cfg.coordinateGuard_m,'all');
                    path=NaN;final=NaN;
                    if ~failed,e=1000*(z-truthWindow);path=sqrt(mean(sum(e.^2,1)));final=norm(e(:,end));end
                    windowScores=[windowScores;table(string(id),string(names{m}),horizon,anchors(w), ...
                        failed,path,final,'VariableNames',{'runId','predictor','horizon_s','anchor', ...
                        'failed','successfulPathVectorRMSE_mm','successfulFinalVectorError_mm'})]; %#ok<AGROW>
                    tr=struct('runId',id,'predictor',names{m},'horizon_s',horizon,'anchor',anchors(w), ...
                        'elapsed_s',reshape(elapsed(:,w,:),1,K),'predicted',z,'truth',truthWindow, ...
                        'pressure',reshape(B.pressure(:,w,:),3,K),'failed',failed,'residual_N',zeros(2,K));
                    if m==3,tr.residual_N=reshape(residual(:,w,:),2,K);end
                    trajectories=[trajectories;tr]; %#ok<AGROW>
                end
            end
        end
    end
    writetable(oneStep,fullfile(folder,'one_step_by_run.csv'));
    writetable(windowScores,fullfile(folder,'rollout_by_window.csv'));
    writetable(windows,fullfile(folder,'rollout_windows.csv'));
    writetable(quality,fullfile(folder,'residual_output_quality.csv'));
    writetable(provenance,fullfile(folder,'heldout_sources.csv'));
    fprintf('Force model external evaluation complete: %s, %d anchors.\n',id,N);
end
save(fullfile(folder,'rollout_predictions.mat'),'trajectories','protocol','-v7.3');
summary=table;
for horizon=protocol.horizons_s
    for m=1:3
        runs=[];failures=0;count=0;
        for id=protocol.runs
            take=windowScores.horizon_s==horizon&windowScores.predictor==names{m}&windowScores.runId==id{1};
            a=windowScores(take,:);failures=failures+nnz(a.failed);count=count+height(a);
            runs(end+1)=sqrt(mean(a.successfulPathVectorRMSE_mm.^2,'omitnan')); %#ok<AGROW>
        end
        summary=[summary;table(string(names{m}),horizon,count,failures,mean(runs,'omitnan'), ...
            'VariableNames',{'predictor','horizon_s','windows','failures','successfulMacroPathVectorRMSE_mm'})]; %#ok<AGROW>
    end
end
writetable(summary,fullfile(folder,'rollout_summary.csv'));
assert(strcmp(filehash(modelPath),modelHash),'Frozen model changed during external evaluation.');
for j=1:height(provenance),assert(strcmp(filehash(provenance.source(j)),provenance.sha256(j)));end
protocol.status='complete';save(fullfile(folder,'protocol.mat'),'protocol','cfg');
copyfile(mfilename('fullpath')+".m",fullfile(folder,'evaluateForceResidualModel.m'));
plots(folder,trainingFolder,oneStep,windowScores,protocol);
disp(oneStep);disp(summary);
fprintf('PASS: frozen model unchanged; November2 used only for evaluation; no normalizer refit or model selection.\n');
end
function plots(folder,trainingFolder,oneStep,windowScores,protocol)
f=figure('Visible','off','Color','w','Position',[80 80 1400 800]);tiledlayout(2,2,'TileSpacing','compact');
history=readtable(fullfile(trainingFolder,'training_history.csv'));
nexttile;plot(history.epoch,history.oldSelectionMacroRMSE_mm,'-o');grid on;xlabel('Epoch');ylabel('Old-trial selection length RMSE (mm)');title('Selection uses older days only');
names={'nominal','constant_velocity','force_gru'};labels={'Nominal, zero residual','Constant velocity','Force-residual GRU'};
for panel=1:3
    nexttile;Y=NaN(5,3);
    for r=1:5
        for m=1:3
            if panel==1
                a=oneStep(oneStep.runId==protocol.runs{r}&oneStep.predictor==names{m},:);Y(r,m)=a.successfulVectorRMSE_mm;
            else
                h=2;if panel==3,h=5;end
                a=windowScores(windowScores.runId==protocol.runs{r}&windowScores.predictor==names{m}&windowScores.horizon_s==h,:);
                if ~any(a.failed),Y(r,m)=sqrt(mean(a.successfulPathVectorRMSE_mm.^2));end
            end
        end
    end
    bar(Y);set(gca,'YScale','log','XTick',1:5,'XTickLabel',protocol.runs);xtickangle(20);grid on;
    ylabel('Length vector RMSE (mm)');
    if panel==1,title('November2: one-step prediction');elseif panel==2,title('November2: 2 s free-running paths');else,title('November2: 5 s free-running paths');end
end
lg=legend(labels,'Location','eastoutside');lg.Layout.Tile='east';
sgtitle('New generalized-force residual model: older-day training, November2 evaluation');drawnow;
savefig(f,fullfile(folder,'force_model_evaluation.fig'));exportgraphics(f,fullfile(folder,'force_model_evaluation.png'),'Resolution',300);
exportgraphics(f,fullfile(folder,'force_model_evaluation.pdf'),'ContentType','vector');close(f);
end
function h=filehash(path)
f=fopen(path,'rb');assert(f>=0);c=onCleanup(@()fclose(f)); %#ok<NASGU>
bytes=fread(f,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
