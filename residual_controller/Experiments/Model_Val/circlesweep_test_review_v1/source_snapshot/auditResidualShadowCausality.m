function report=auditResidualShadowCausality(root,folder)
%AUDITRESIDUALSHADOWCAUSALITY Offline causality/future-target falsification.
% Saved measurements/commands only. No training, SLX changes or device access.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,folder=fullfile(root,'residual_controller','Experiments','Model_Val','circlesweep_test_review_v1');end
target=fullfile(folder,'causality_audit');assert(~isfolder(target),'Choose a fresh audit folder.');mkdir(target);
addpath(root,fullfile(root,'residual_controller'),fullfile(root,'passivity_adaptive_control'), ...
    fullfile(root,'residual_dynamics_learning','force_model'));
A=load(fullfile(folder,'analysis.mat'));S=A.S;D=A.D;C=A.C;t=A.t;N=numel(t);
q=S.meas_length(:,2:3);p=A.pressure;physical=logical(S.orientationValid)&logical(S.baseOreintationValid);
options=residualShadowOptions();options(1)=unique(D.units);options(2)=unique(C.horizon_s(logical(C.valid)));
network=residualControllerNetwork(options(1));
report.source=fullfile(root,'residual_controller','Experiments','Model_Val','circlesweep_test.mat');
report.replayOptions=options';report.replayOptionsSource='Current source defaults; width/horizon verified from log. Full options vector was not saved.';
state=residualShadowInitialState();P=nan(N,6);R=nan(N,6);F=nan(N,2);status=zeros(N,1);ready=false(N,1);
selected=find(t(1:end-1)>=30&vecnorm(D.velocity_mPerS(1:end-1,:),2,2)>.001& ...
    logical(C.valid(1:end-1))&A.scored(2:end),1);
assert(~isempty(selected));chosenAfter=residualShadowInitialState();
for k=1:N
    [pred,matched,~,d,~,state]=callStep(k,state,S,q,p,t,D,options,network);
    P(k,:)=pred(:)';R(k,:)=matched(:)';F(k,:)=d.rawForce_N';status(k)=d.status;ready(k)=d.ready;
    if k==selected,chosenAfter=state;end
end
originalP=reshape(A.issued,N,6);originalR=reshape(A.Q,N,6);
report.replayIssuedDifference_m=max(abs(P-originalP),[],'all');
report.replayMatchedDifference_m=max(abs(R-originalR),[],'all');
report.replayForceDifference_N=max(abs(F-D.rawForce_N),[],'all');
report.replayStatusMismatches=nnz(status~=D.status);report.replayReadyMismatches=nnz(ready~=logical(D.ready));
report.replayIssuedFiniteMaskMismatches=nnz(isfinite(P)~=isfinite(originalP));
assert(report.replayIssuedDifference_m<1e-10&&report.replayMatchedDifference_m<1e-10&& ...
    report.replayForceDifference_N<1e-8&&report.replayStatusMismatches==0&& ...
    report.replayReadyMismatches==0&&report.replayIssuedFiniteMaskMismatches==0);
% Change ONLY the next observed target. The pending prediction must stay fixed;
% error must change by minus the injected offset. Raw MAT data stay untouched.
k=selected+1;changed=q;offset=[.005;-.003];changed(k,:)=changed(k,:)+offset';
[nextOriginal,matchedOriginal,~,dOriginal]=callStep(k,chosenAfter,S,q,p,t,D,options,network);
[nextChanged,matchedChanged,~,dChanged]=callStep(k,chosenAfter,S,changed,p,t,D,options,network);
report.counterfactualOrigin_sim_s=t(selected);report.counterfactualTarget_sim_s=t(k);
report.futureMeasurementOffset_mm=1000*offset';
report.previousPredictionChange_m=max(abs(matchedChanged-matchedOriginal),[],'all');
report.targetErrorChange_mm=1000*(dChanged.errorAtMeasurement_m(:,1)-dOriginal.errorAtMeasurement_m(:,1))';
report.newForecastChange_mm=1000*(nextChanged(:,1)-nextOriginal(:,1))';
assert(report.previousPredictionChange_m==0&&all(abs(report.targetErrorChange_mm+1000*offset')<1e-8));

% Persistence (simply repeat the latest observed length) is another strong
% short-horizon reference. The target itself is never used as its prediction.
origins=A.sourceRows;targets=find(A.scored);baseline=q(origins,:);
e=1000*(baseline-q(targets,:));report.persistenceOneStepRMSE_mm=sqrt(mean(e.^2,1));
moving=vecnorm(D.velocity_mPerS(targets,:),2,2)>.001;
report.persistenceMovingRMSE_mm=sqrt(mean(e(moving,:).^2,1));

% Conditional recorded-input rollouts. Future COMMANDS and elapsed intervals
% are known offline. Future NDI q/v are used ONLY for scoring after propagation.
steps=[1 2 5 10];M=max(steps);eligible=[];
for origin=find(logical(C.valid)&S.trialId>0)'
    if origin+M>N,continue;end
    rows=origin:origin+M;dt=diff(S.pollStart_s(rows));
    if all(physical(rows))&&all(isfinite(p(rows,:)),'all')&&all(dt>0&dt<=options(6))&& ...
            all(diff(t(rows))>0&diff(t(rows))<=options(7)+1e-9)&&all(S.trialId(rows)==S.trialId(origin))
        eligible(end+1,1)=origin; %#ok<AGROW>
    end
end
anchors=eligible(1:20:end);assert(~isempty(anchors));count=numel(anchors);
predictions=nan(count,2,4,numel(steps));truth=nan(count,2,numel(steps));
failed=false(count,4,numel(steps));elapsed=zeros(count,numel(steps));
for j=1:count
    origin=anchors(j);ctx=residualShadowContext();ctx.valid=true;ctx.measured=C.measured(origin,:)';
    ctx.historyPre=C.historyPre(origin,:)';ctx.featureDt_s=C.featureDt_s(origin);
    commands=p(origin:origin+M-1,:)';dt=diff(S.pollStart_s(origin:origin+M));
    [pred,bad]=rollout(ctx,commands,dt,network,steps,options(10));
    predictions(j,:,:,:)=reshape(pred,1,2,4,numel(steps));failed(j,:,:)=reshape(bad,1,4,numel(steps));
    for h=1:numel(steps)
        truth(j,:,h)=q(origin+steps(h),:);elapsed(j,h)=sum(dt(1:steps(h)));
    end
end
metrics=table;labels=["residual_plus_nominal","nominal","constant_velocity","persistence"];
anchorMoving=vecnorm(C.measured(anchors,3:4),2,2)>.001;
for h=1:numel(steps)
    for subset=["all_anchors","moving_at_origin"]
        select=true(count,1);if subset=="moving_at_origin",select=anchorMoving;end
        for model=1:4
            valid=select&~failed(:,model,h);e=1000*(predictions(valid,:,model,h)-truth(valid,:,h));rmse=sqrt(mean(e.^2,1));
            metrics=[metrics;table(steps(h)*options(2),median(elapsed(select,h)),subset,labels(model), ...
                nnz(select),nnz(valid),nnz(select&failed(:,model,h)),rmse(1),rmse(2), ...
                'VariableNames',{'nominalHorizon_s','medianActualHorizon_s','subset','predictor', ...
                'eligibleAnchors','scoredAnchors','guardFailures','l2RMSE_mm','l3RMSE_mm'})]; %#ok<AGROW>
        end
    end
end
report.rolloutAnchors=count;report.rolloutHorizonSteps=steps;
report.rolloutProtocol='Every20th eligible positive-trial anchor. Known recorded future commands/elapsed intervals; no future measured state correction. Guard .05m; failures counted, never clipped.';
report.rolloutFutureMeasurementsUsedAsInputs=false;
writetable(metrics,fullfile(target,'conditional_rollout_metrics.csv'));
f=fopen(fullfile(target,'causality_summary.json'),'w');fwrite(f,jsonencode(report,PrettyPrint=true));fclose(f);
save(fullfile(target,'audit.mat'),'report','metrics','P','R','F','status','ready','selected','offset', ...
    'matchedOriginal','matchedChanged','nextOriginal','nextChanged','dOriginal','dChanged', ...
    'anchors','predictions','truth','failed','elapsed');
fig=figure('Visible','off','Color','w','Position',[30 30 1150 650]);tiledlayout(1,2);
for coordinate=1:2
    nexttile;hold on;
    for model=1:4
        use=metrics.subset=="moving_at_origin"&metrics.predictor==labels(model);
        y=metrics.l2RMSE_mm(use);if coordinate==2,y=metrics.l3RMSE_mm(use);end
        plot(metrics.nominalHorizon_s(use),y,'-o','LineWidth',1.4);
    end
    xlabel('Nominal horizon (s)');ylabel(sprintf('l%d RMSE (mm)',coordinate+1));grid on;
end
legend('Residual + nominal','Nominal','Constant velocity','Persistence','Location','eastoutside');
sgtitle('No measurement resets: conditional recorded-input rollout, moving anchors');
exportgraphics(fig,fullfile(target,'conditional_rollouts.png'),'Resolution',180);
exportgraphics(fig,fullfile(target,'conditional_rollouts.pdf'),'ContentType','vector');savefig(fig,fullfile(target,'conditional_rollouts.fig'));close(fig);
disp(report);disp(metrics);
end

function [prediction,matched,retimed,d,context,state]=callStep(k,state,S,q,p,t,D,options,network)
stamp=t(k)-D.commandAge_sim_s(k);if ~isfinite(stamp),stamp=t(k);end
lengths=[-sum(q(k,:));q(k,:)'];
[prediction,matched,retimed,d,context,state]=residualShadowStep(lengths,logical(S.orientationValid(k)), ...
    logical(S.baseOreintationValid(k)),S.pollStart_s(k),S.frame(k,2),p(k,:)',stamp,t(k),false,options,state,network);
end

function [prediction,failed]=rollout(context,commands,dt,network,steps,guard)
q0=context.measured(1:2);v0=context.measured(3:4);qa=q0;va=v0;qn=q0;vn=v0;
h=context.historyPre;featureDt=context.featureDt_s;cfg=residualShadowPhysicsConfig();
prediction=nan(2,4,numel(steps));failed=false(4,numel(steps));badA=false;badN=false;elapsed=0;
for j=1:max(steps)
    if ~badA
        [h,force]=residualGruCell(network,[qa;va;commands(:,j);featureDt],h);
        [qa,va]=forceResidualStep(qa,va,commands(:,j),force,dt(j),cfg);
        badA=any(~isfinite([qa;va;h]))||any(abs(qa)>guard);
    end
    if ~badN
        [qn,vn]=forceResidualStep(qn,vn,commands(:,j),zeros(2,1),dt(j),cfg);
        badN=any(~isfinite([qn;vn]))||any(abs(qn)>guard);
    end
    elapsed=elapsed+dt(j);featureDt=dt(j);horizon=find(steps==j,1);
    if ~isempty(horizon)
        candidate=[qa qn q0+elapsed*v0 q0];bad=[badA badN any(abs(candidate(:,3))>guard) false];
        failed(:,horizon)=bad';candidate(:,bad)=NaN;prediction(:,:,horizon)=candidate;
    end
end
end
