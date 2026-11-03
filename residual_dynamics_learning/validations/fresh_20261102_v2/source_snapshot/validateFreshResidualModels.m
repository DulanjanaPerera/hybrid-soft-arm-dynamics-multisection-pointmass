function validateFreshResidualModels(root,outFolder)
%VALIDATEFRESHRESIDUALMODELS Frozen 20261101 models on new 20261102 exports.
% Offline only. No fitting, model selection, Simulink, serial or DAQ calls.
assert(~isfolder(outFolder),'Use a fresh output directory.');
addpath(root,fullfile(root,'residual_dynamics_learning'));
dataset=fullfile(root,'residual_dynamics_learning','datasets','two_day_v2');
screen=fullfile(root,'residual_dynamics_learning','benchmarks','screen_20261101_v2');
old=load(fullfile(dataset,'preparation.mat'));n=load(fullfile(dataset,'normalizer.mat'));
f=load(fullfile(screen,'forward','frozen_selection.mat'));
i=load(fullfile(screen,'inverse','frozen_selection.mat'));
fa=load(fullfile(screen,'ablations','forward_frozen.mat'));
ia=load(fullfile(screen,'ablations','inverse_frozen.mat'));
cfg=old.cfg;
if ~isfield(cfg,'boundMode'),cfg.boundMode=0;end % Original saved preparation predates this explicit field.
assert(cfg.boundMode==0);
protocol=struct('status','frozen_before_fresh_evaluation','fit',false, ...
    'runs',{{'StepWave','SineWave4','combinePressure','circleSweep','helix'}}, ...
    'horizons_s',[.5 2 5],'maxRolloutsPerRunHorizon',12, ...
    'coordinateGuard_m',.05,'historicalSelectedForward',f.selected, ...
    'baseline','legacy mode0 context; corrected mode1 nominal reported separately', ...
    'inputApproximation','NDI-row pressure ZOH; approximate command/physical clock mapping', ...
    'stateUpdate','predicted q; causal velocity update tau=.20s, no future-motion resets', ...
    'partition','all five recordings fresh_test; no normalization or fitting');
mkdir(outFolder);save(fullfile(outFolder,'protocol.mat'),'protocol','cfg');
files={fullfile(dataset,'normalizer.mat'),fullfile(dataset,'preparation.mat'), ...
    fullfile(screen,'forward','frozen_selection.mat'), ...
    fullfile(screen,'inverse','frozen_selection.mat'), ...
    fullfile(screen,'ablations','forward_frozen.mat'), ...
    fullfile(screen,'ablations','inverse_frozen.mat')};
sources={'validateFreshResidualModels.m','residualBatch.m','residualPredict.m', ...
    'residualFrozenStep.m','residualNominal.m','armS_single_boundForce.m', ...
    'armS_single_core.m','armS_single_dynamics.m','integratedPosition_nume.m', ...
    'integratedJacobianProduct_compact.m','christoffelSymbol.m'};
mkdir(fullfile(outFolder,'source_snapshot'));
for j=1:numel(sources)
    p=which(sources{j});assert(~isempty(p),'Missing dependency %s',sources{j});
    files{end+1}=p;copyfile(p,fullfile(outFolder,'source_snapshot',sources{j})); %#ok<AGROW>
end
locked=string(files(:));before=cellfun(@filehash,cellstr(locked),'UniformOutput',false);
writetable(table(locked,string(before),'VariableNames',{'source','sha256'}), ...
    fullfile(outFolder,'frozen_inputs.csv'));
quality=table;forward=table;inverse=table;rollouts=table;traces=table;provenance=table;
for run=1:numel(protocol.runs)
    id=protocol.runs{run};folder=fullfile(root,'RL_data_collection','20261102','exports',id);
    assert(isfolder(folder));meta=jsondecode(fileread(fullfile(folder,'trial_manifest.json')));
    raw=meta.source_mat;assert(isfile(raw));
    O=readtable(fullfile(folder,'observations.csv'),'VariableNamingRule','preserve');
    S=prepare(O,cfg,id);save(fullfile(outFolder,[id '_sequences.mat']),'S','cfg');
    idx=S.forwardAnchors;assert(~isempty(idx),'No eligible windows for %s',id);
    rows=table(S.t,S.physicalTime,S.trial,S.channel,S.phase,S.valid,S.segment,S.mature, ...
        'VariableNames',{'time_sim_s','time_physical_s','trial','channel','phase','valid','segment','mature'});
    writetable(rows,fullfile(outFolder,[id '_rows.csv']));
    writetable(table(idx,S.trial(idx),S.segment(idx),'VariableNames', ...
        {'anchorRow','trial','segment'}),fullfile(outFolder,[id '_anchors.csv']));
    values=[S.q(S.mature,:) S.v(S.mature,:) S.p(S.mature,:)];
    zz=(values-n.normalizer.mean)./n.normalizer.std;
    quality=[quality;table(string(id),height(O),nnz(O.orientation_valid~=1), ...
        nnz(S.valid),max(S.segment),nnz(S.mature),numel(idx),median(S.dt(2:end)), ...
        nnz(S.dt>.10),isfield(meta.quality,'ndi_diagnostics_logged')&&meta.quality.ndi_diagnostics_logged, ...
        max(abs(zz),[],'all'),mean(any(abs(zz)>5,2)), ...
        'VariableNames',{'runId','rows','orientationInvalidRows','validRows','segments', ...
        'matureRows','anchors','medianPhysicalDt_s','pollStartGapsAbovePoint1_s', ...
        'diagnosticsLogged','maxAbsoluteTrainingZ','fractionRowsAnyTrainingZAbove5'})]; %#ok<AGROW>
    for source={raw,fullfile(folder,'observations.csv'),fullfile(folder,'commands.csv'), ...
            fullfile(folder,'trial_manifest.json')}
        provenance=[provenance;table(string(id),string(source{1}),string(filehash(source{1})), ...
            'VariableNames',{'runId','source','sha256'})]; %#ok<AGROW>
    end
    % All indices preserve the exact original history/filter/maturity contract.
    B=residualBatch(S,idx,'forward',cfg,n.normalizer);N=numel(idx);
    D=struct('Seq',single(B.Xsequence),'Context',zeros(N,3));base=zeros(N,2);corrected=base;
    for j=1:N
        k=idx(j);z=residualFrozenStep(S.q(k,:)',S.v(k,:)',S.p(k,:)',B.dtNext(j),0);
        base(j,:)=z(1:2)';z=residualFrozenStep(S.q(k,:)',S.v(k,:)',S.p(k,:)',B.dtNext(j),1);
        corrected(j,:)=z(1:2)';
    end
    D.Context=[B.dtNext base];D.Base=base;D.Truth=B.Y;D.Y=1e9*ones(N,2);
    preds={base,corrected,S.q(idx,:),S.q(idx,:)+S.v(idx,:).*B.dtNext};
    names={'nominal_legacy','nominal_corrected','persistence','constant_velocity'};
    for j=1:numel(f.bank)
        names{end+1}=f.bank(j).name;preds{end+1}=base+residualPredict(f.bank(j),D)/1000; %#ok<AGROW>
    end
    DK=D;DK.Base=preds{4};DK.Context=[B.dtNext DK.Base];
    names{end+1}='kinematic_ridge';preds{end+1}=DK.Base+residualPredict(fa.bank(1),DK)/1000;
    % Confirm inference cannot read targets; future motion cannot alter inputs.
    small=D;small.Seq=D.Seq(:,:,1:3);small.Context=D.Context(1:3,:);small.Base=D.Base(1:3,:);
    small.Truth=D.Truth(1:3,:);small.Y=D.Y(1:3,:);
    for j=1:numel(f.bank)
        a=residualPredict(f.bank(j),small);small.Truth(:)=-1e9;small.Y(:)=NaN;
        assert(isequaln(a,residualPredict(f.bank(j),small)));
    end
    changed=S;k=idx(1);changed.q(k+1:k+cfg.previewFrames,:)=1e6;
    a=residualBatch(S,k,'forward',cfg,n.normalizer);b=residualBatch(changed,k,'forward',cfg,n.normalizer);
    assert(isequaln(a.Xsequence,b.Xsequence)&&all(isnan(a.preview(:))));
    for j=1:numel(names),forward=[forward;metrics(names{j},id,S,idx,preds{j},B.Y,1000)];end %#ok<AGROW>
    % Conditional inverse uses observed motion preview ONLY as offline requirement.
    BI=residualBatch(S,idx,'inverse',cfg,n.normalizer);dn=S.inverseDifferentialNominal(idx,:);
    DI=struct('Seq',single(BI.Xsequence),'Context',[reshape(BI.preview,[],N)' BI.p0 dn], ...
        'Base',dn);map=[-1/3 -1/3;2/3 -1/3;-1/3 2/3];
    inames={'nominal_inverse','previous_differential'};
    ipreds={BI.p0+dn*map',BI.p0+S.d(idx-1,:)*map'};
    ibank=[i.bank(:);ia.bank(:)];
    for j=1:numel(ibank)
        inames{end+1}=ibank(j).name;ipreds{end+1}=BI.p0+(dn+residualPredict(ibank(j),DI))*map'; %#ok<AGROW>
    end
    changed=S;changed.p(k,:)=100;
    a=residualBatch(S,k,'inverse',cfg,n.normalizer);b=residualBatch(changed,k,'inverse',cfg,n.normalizer);
    assert(isequaln(a.Xsequence,b.Xsequence));
    for j=1:numel(inames)
        inverse=[inverse;metrics(inames{j},id,S,idx,ipreds{j},S.p(idx,:),1)]; %#ok<AGROW>
    end
    save(fullfile(outFolder,[id '_one_step.mat']),'idx','names','preds','inames','ipreds','B');
    [rr,tt]=rollout(S,cfg,n.normalizer,f.bank,fa.bank,protocol);
    rollouts=[rollouts;rr];traces=[traces;tt]; %#ok<AGROW>
    writetable(quality,fullfile(outFolder,'quality.csv'));
    writetable(forward,fullfile(outFolder,'forward_by_trial.csv'));
    writetable(inverse,fullfile(outFolder,'inverse_by_trial.csv'));
    writetable(rollouts,fullfile(outFolder,'rollouts.csv'));
    writetable(traces,fullfile(outFolder,'rollout_traces.csv'));
    writetable(provenance,fullfile(outFolder,'fresh_sources.csv'));
    fprintf('Fresh validation complete: %s, %d guarded anchors.\n',id,N);
end
after=cellfun(@filehash,cellstr(locked),'UniformOutput',false);
assert(isequal(before,after),'Frozen models, normalizer or source changed.');
for j=1:height(provenance),assert(strcmp(filehash(provenance.source(j)),provenance.sha256(j)));end
summary=aggregate(forward);writetable(summary,fullfile(outFolder,'forward_summary.csv'));
inverseSummary=aggregate(inverse);writetable(inverseSummary,fullfile(outFolder,'inverse_summary.csv'));
figures(outFolder,forward,rollouts,traces,protocol.runs);
protocol.status='complete';save(fullfile(outFolder,'protocol.mat'),'protocol','cfg');
disp(summary);disp(inverseSummary);
fprintf('PASS: same-segment mature windows; forward target and inverse-current-pressure isolation; frozen-file/raw/export hashes unchanged; no fitting.\n');
end

function S=prepare(O,cfg,id)
S=struct;S.t=O.time_sim_s;S.physicalTime=O.sensor_poll_start_host_s;
S.dt=[NaN;diff(S.physicalTime)];assert(all(diff(S.t)>0));
S.q=O{:,{'length_inferred_m_2','length_inferred_m_3'}};
S.p=O{:,{'p1_command_bar','p2_command_bar','p3_command_bar'}};
S.trial=O.trialId;S.phase=O.phaseId;S.channel=zeros(height(O),1);
if ismember('channelId',O.Properties.VariableNames),S.channel=O.channelId;end
S.split=repmat("fresh_test",height(O),1);S.runId=id;S.cohort='day3_20261102';
ages=O{:,{'p1_command_bar_last_log_age_sim_s','p2_command_bar_last_log_age_sim_s','p3_command_bar_last_log_age_sim_s'}};
S.valid=O.dynamics_observation_valid==1&O.base_orientation_valid==1& ...
    all(isfinite([S.q S.p]),2)&all(ages>=-1e-9&ages<=cfg.maxCommandAge_s,2)& ...
    all(S.p>=-1e-9&S.p<=3+1e-9,2)&isfinite(S.physicalTime)&isfinite(S.trial)&S.trial>0;
N=height(O);S.segment=zeros(N,1);S.mature=false(N,1);S.v=NaN(N,2);S.a=S.v;seg=0;start=0;
for k=1:N
    if ~S.valid(k),continue;end
    continuous=k>1&&S.valid(k-1)&&S.trial(k)==S.trial(k-1)&&S.channel(k)==S.channel(k-1)&& ...
        S.dt(k)>0&&S.dt(k)<=cfg.maxPhysicalGap_s&&S.t(k)-S.t(k-1)<=cfg.maxSimulationGap_s;
    if ~continuous,seg=seg+1;start=k;S.v(k,:)=0;S.a(k,:)=0;
    else
        dt=S.dt(k);raw=(S.q(k,:)-S.q(k-1,:))/dt;
        S.v(k,:)=S.v(k-1,:)+dt/(cfg.velocityTau_s+dt)*(raw-S.v(k-1,:));
        raw=(S.v(k,:)-S.v(k-1,:))/dt;
        S.a(k,:)=S.a(k-1,:)+dt/(cfg.accelerationTau_s+dt)*(raw-S.a(k-1,:));
    end
    S.segment(k)=seg;S.mature(k)=S.physicalTime(k)-S.physicalTime(start)>=cfg.warmup_s;
end
S.inverseDifferentialNominal=NaN(N,2);
for k=find(S.mature)'
    [~,~,d]=residualNominal(S.q(k,:)',S.v(k,:)',S.a(k,:)',S.p(k,:)',0);
    S.inverseDifferentialNominal(k,:)=d';
end
S.p0=mean(S.p,2);S.d=[S.p(:,2)-S.p(:,1),S.p(:,3)-S.p(:,1)];
S.inverseResidualCandidate=S.d-S.inverseDifferentialNominal;
idx=[];
for k=cfg.historyFrames+1:N-cfg.previewFrames
    rows=k-cfg.historyFrames:k+cfg.previewFrames;
    if S.segment(k)>0&&all(S.mature(rows))&&all(S.segment(rows)==S.segment(k)),idx(end+1,1)=k;end %#ok<AGROW>
end
S.forwardAnchors=idx;S.inverseAnchors=idx;
end

function T=metrics(name,id,S,idx,pred,truth,scale)
T=table;
for trial=[0;unique(S.trial(idx))]'
    keep=true(numel(idx),1);if trial~=0,keep=S.trial(idx)==trial;end
    e=scale*(pred(keep,:)-truth(keep,:));assert(all(isfinite(e),'all'));
    rmse=sqrt(mean(e.^2,1));vector=sqrt(mean(sum(e.^2,2)));third=NaN;
    if size(e,2)==3,third=rmse(3);end
    gap=NaN;frac=NaN;
    if scale==1,g=max(max(-pred(keep,:),pred(keep,:)-3),0);gap=max(g,[],'all');frac=mean(any(g>.01,2));end
    T=[T;table(string(name),string(id),trial,nnz(keep),rmse(1),rmse(2),third,vector, ...
        gap,frac,'VariableNames',{'model','runId','trial','anchors','coordinate1RMSE', ...
        'coordinate2RMSE','coordinate3RMSE','vectorRMSE','maximumPressureViolation_bar','pressureViolationAbovePoint01Fraction'})]; %#ok<AGROW>
end
end

function T=aggregate(detail)
R=detail(detail.trial==0,:);T=table;
for name=unique(R.model,'stable')'
    a=R(R.model==name,:);T=[T;table(name,height(a),sum(a.anchors),mean(a.vectorRMSE), ...
        sqrt(sum(a.anchors.*a.vectorRMSE.^2)/sum(a.anchors)), ...
        'VariableNames',{'model','runs','anchors','macroVectorRMSE','pooledVectorRMSE'})]; %#ok<AGROW>
end
end

function [report,traces]=rollout(S,cfg,n,bank,kin,protocol)
report=table;traces=table;H=cfg.historyFrames;
names=[{'nominal_legacy','nominal_corrected','persistence','constant_velocity'}, ...
    reshape(arrayfun(@(m)m.name,bank,'UniformOutput',false),1,[]),{'kinematic_ridge'}];
names=reshape(names,1,[]);
for horizon=protocol.horizons_s
    eligible=[];ends=[];
    for k=S.forwardAnchors'
        last=find(S.physicalTime>=S.physicalTime(k)+horizon,1);
        if isempty(last),continue;end
        rows=k-H:last;
        if all(S.mature(rows))&&all(S.segment(rows)==S.segment(k)),eligible(end+1)=k;ends(end+1)=last;end %#ok<AGROW>
    end
    if isempty(eligible),continue;end
    take=unique(round(linspace(1,numel(eligible),min(protocol.maxRolloutsPerRunHorizon,numel(eligible)))));
    anchors=eligible(take);lasts=ends(take);N=numel(anchors);steps=lasts-anchors;
    for mi=1:numel(names)
        name=names{mi};Q=zeros(2,H,N);V=Q;errors=cell(N,1);failed=false(N,1);
        for j=1:N,rows=anchors(j)-H+1:anchors(j);Q(:,:,j)=S.q(rows,:)';V(:,:,j)=S.v(rows,:)';errors{j}=NaN(steps(j),2);end
        for step=1:max(steps)
            active=find(~failed&step<=steps(:));if isempty(active),break;end
            P=struct('Seq',zeros(8,H,numel(active),'single'),'Context',zeros(numel(active),3));bases=zeros(numel(active),2);
            for a=1:numel(active)
                j=active(a);k=anchors(j)+step-1;dt=S.dt(k+1);rows=k-H+1:k;
                mode=double(strcmp(name,'nominal_corrected'));
                z=residualFrozenStep(Q(:,end,j),V(:,end,j),S.p(k,:)',dt,mode);bases(a,:)=z(1:2)';
                if strcmp(name,'kinematic_ridge'),bases(a,:)=(Q(:,end,j)+V(:,end,j)*dt)';end
                values=[Q(:,:,j)' V(:,:,j)' S.p(rows,:)];
                P.Seq(:,:,a)=[(values-n.mean)./n.std S.dt(rows)]';P.Context(a,:)=[dt bases(a,:)];
            end
            R=zeros(numel(active),2);
            if mi>4&&mi<=4+numel(bank),R=residualPredict(bank(mi-4),P);
            elseif strcmp(name,'kinematic_ridge'),R=residualPredict(kin(1),P);end
            for a=1:numel(active)
                j=active(a);k=anchors(j)+step-1;dt=S.dt(k+1);
                if strcmp(name,'persistence'),next=Q(:,end,j);
                elseif strcmp(name,'constant_velocity'),next=Q(:,end,j)+V(:,end,j)*dt;
                else,next=bases(a,:)'+R(a,:)'/1000;end
                if any(~isfinite(next))||any(abs(next)>protocol.coordinateGuard_m),failed(j)=true;continue;end
                velocity=V(:,end,j)+dt/(cfg.velocityTau_s+dt)*((next-Q(:,end,j))/dt-V(:,end,j));
                errors{j}(step,:)=1000*(next'-S.q(k+1,:));
                if j==1
                    tr=table(string(S.runId),string(name),horizon,anchors(j),k+1, ...
                        S.physicalTime(k+1)-S.physicalTime(anchors(j)), ...
                        S.q(k+1,1)*1000,S.q(k+1,2)*1000,next(1)*1000,next(2)*1000, ...
                        'VariableNames',{'runId','model','horizon_s','anchorRow','targetRow', ...
                        'elapsed_s','trueQ2_mm','trueQ3_mm','predQ2_mm','predQ3_mm'});
                    traces=[traces;tr]; %#ok<AGROW>
                end
                Q(:,:,j)=[Q(:,2:end,j) next];V(:,:,j)=[V(:,2:end,j) velocity];
            end
        end
        path=NaN;final=NaN;success=errors(~failed);
        if ~isempty(success)
            e=vertcat(success{:});path=sqrt(mean(sum(e.^2,2)));
            e=cellfun(@(a)a(end,:),success,'UniformOutput',false);e=vertcat(e{:});final=sqrt(mean(sum(e.^2,2)));
        end
        duration=S.physicalTime(lasts)-S.physicalTime(anchors);
        report=[report;table(string(name),string(S.runId),horizon,N,nnz(failed), ...
            min(duration),max(duration),path,final,'VariableNames',{'model','runId', ...
            'horizon_s','rollouts','failures','minActualHorizon_s','maxActualHorizon_s', ...
            'successfulPathVectorRMSE_mm','successfulFinalVectorRMSE_mm'})]; %#ok<AGROW>
    end
end
end

function figures(folder,F,R,T,ids)
f=figure('Visible','off','Color','w','Position',[80 80 1200 650]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;A=F(F.trial==0,:);names=unique(A.model,'stable');Y=NaN(numel(ids),numel(names));
for j=1:numel(ids),for k=1:numel(names),a=A(A.runId==ids{j}&A.model==names(k),:);Y(j,k)=a.vectorRMSE;end,end
bar(Y);set(gca,'YScale','log','XTick',1:numel(ids),'XTickLabel',ids);xtickangle(25);grid on;
ylabel('One-step length vector RMSE (mm)');title('Frozen models: fresh recordings');
legend(strrep(cellstr(names),'_',' '),'Location','southoutside','NumColumns',2);
nexttile;show={'nominal_legacy','constant_velocity','ridge','mlp','gru','kinematic_ridge'};hold on;
for j=1:numel(show)
    a=R(R.model==show{j},:);hs=unique(a.horizon_s);vals=NaN(size(hs));
    for k=1:numel(hs),v=a(a.horizon_s==hs(k),:);if all(v.failures==0),vals(k)=mean(v.successfulPathVectorRMSE_mm);end,end
    plot(hs,vals,'-o','DisplayName',strrep(show{j},'_',' '));
end
set(gca,'YScale','log');grid on;xlabel('Physical prediction horizon (s)');ylabel('Macro path vector RMSE (mm)');
title('No measured-state resets; failure counts in CSV');legend('Location','southoutside','NumColumns',2);
savefig(f,fullfile(folder,'fresh_prediction_summary.fig'));exportgraphics(f,fullfile(folder,'fresh_prediction_summary.png'),'Resolution',300);
exportgraphics(f,fullfile(folder,'fresh_prediction_summary.pdf'),'ContentType','vector');close(f);
f=figure('Visible','off','Color','w','Position',[80 80 1400 750]);tiledlayout(2,3,'TileSpacing','compact');
for j=1:numel(ids)
    nexttile;hold on;A=T(T.runId==ids{j}&T.horizon_s==5,:);
    first=A(A.model=="constant_velocity",:);plot(first.elapsed_s,first.trueQ3_mm,'k-','LineWidth',1.6,'DisplayName','Measured l3');
    for name={'constant_velocity','ridge','kinematic_ridge'}
        a=A(A.model==name{1},:);plot(a.elapsed_s,a.predQ3_mm,'DisplayName',strrep(name{1},'_',' '));
    end
    grid on;title(ids{j});xlabel('Physical elapsed time (s)');ylabel('l3 (mm)');
end
nexttile;hold on;plot(NaN,NaN,'k-','LineWidth',1.6);plot(NaN(2,3));axis off;
legend('Measured l3','constant velocity','ridge','kinematic ridge','Location','northwest');
sgtitle('First deterministic eligible 5 s window per recording; full q2/q3 traces in CSV');
savefig(f,fullfile(folder,'fresh_rollout_examples.fig'));exportgraphics(f,fullfile(folder,'fresh_rollout_examples.png'),'Resolution',300);
exportgraphics(f,fullfile(folder,'fresh_rollout_examples.pdf'),'ContentType','vector');close(f);
end

function h=filehash(path)
fid=fopen(path,'rb');assert(fid>=0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(bytes,'int8'));h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
