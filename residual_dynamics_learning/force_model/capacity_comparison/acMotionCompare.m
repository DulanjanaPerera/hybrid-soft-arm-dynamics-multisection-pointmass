function acMotionCompare(root,renderOnly)
%ACMOTIONCOMPARE Fixed Nov2 examples and exclusively held-out Nov3 trials.
base=fullfile(root,'residual_dynamics_learning');force=fullfile(base,'force_model');src=fullfile(force,'capacity_comparison');
addpath(base,force,src,fullfile(force,'streaming'),fullfile(force,'fine_tuning'));
camp=fullfile(src,'campaign_20261103_v1');F=load(fullfile(camp,'FROZEN_BEFORE_EVALUATION.mat'));F.results=F.representatives;
folder=fullfile(camp,'motion_comparison');
if nargin>=2&&renderOnly
    saved=load(fullfile(folder,'predictions.mat'));renderBundle(camp,saved.trajectories,saved.cfg,F);
    copyfile(mfilename('fullpath')+".m",fullfile(folder,'acMotionCompare.m'));return;
end
assert(~isfolder(folder));mkdir(folder);
D=load(fullfile(force,'fine_tuning','campaign_20261103_v1','data','data.mat'));models=cell(height(F.results),1);
for a=1:numel(models)
    assert(strcmp(sgHash(F.results.modelPath(a)),F.results.sha256(a)));
    b=load(F.results.modelPath(a));models{a}=b.model;
end
scores=table;selection=table;trajectories=struct([]);cfg=models{1}.cfg;
for day=2:3
    if day==2
        p=load(fullfile(base,'figures','motion_predictions_20261102_v1','predictions.mat'));W=p.windows;
        circle=load(fullfile(base,'validations','fresh_20261102_v2','circleSweep_sequences.mat'));
        sine=load(fullfile(base,'validations','fresh_20261102_v2','SineWave4_sequences.mat'));
        helix=load(fullfile(base,'validations','fresh_20261102_v2','helix_sequences.mat'));
        sequences={circle.S,sine.S,helix.S};
    else
        p=load(fullfile(base,'figures','motion_predictions_20261103_v1','frozen_model_predictions_v1','predictions.mat'));W=p.windows;
        sequences=D.sequences(11:13);
        % Original Nov3 illustration used TRAIN helix trial1. Here use the
        % declared whole-trial holdout2 for every candidate, including baseline.
        S=sequences{3};O=sgRead(fullfile(base,'figures','motion_predictions_20261103_v1','exports','helix','observations.csv'));
        candidates=find(S.split=="selection"&S.phase==6&S.mature);first=[];
        for k=candidates'
            mature=find(S.mature&S.segment==S.segment(k),1);
            if S.physicalTime(k)-S.physicalTime(mature)>=cfg.commonBurnIn_s,first=k;break;end
        end
        assert(~isempty(first));last=find(S.trial==S.trial(first)&S.phase==7,1,'last');rows=(first:last)';
        A=W(5);A.rows=rows;A.t=S.t(rows)-S.t(first);A.physicalElapsed=S.physicalTime(rows)-S.physicalTime(first);
        A.truth=S.q(rows,:);A.truth(~S.valid(rows),:)=NaN;A.pressure=S.p(rows,:);A.missing=~S.valid(rows);A.trial=S.trial(first);
        A.thetaReference=O.theta_reference_rad(rows);A.phiReference=O.phi_reference_rad(rows);A.thetaObservable=O.theta_valid(rows)>.5;
        A.description='first eligible loading/unloading remainder in declared Nov3 heldout helix trial';W(5)=A;
    end
    % Keep only observed window metadata; legacy predictions belong to the
    % original example, especially Nov3 helix trial1 replaced by holdout2 here.
    W=rmfield(W,intersect(fieldnames(W),{'free','one','residual_N','completed','stopAt'}));
    for w=1:numel(W)
        A=W(w);S=sequences{A.kind};first=A.rows(1);last=A.rows(end);N=numel(A.rows);variants=height(F.results);
        assert(all(S.dt(first+1:last)>0&S.dt(first+1:last)<=cfg.maxPhysicalGap_s));
        assert(all(diff(S.t(first:last))>0&diff(S.t(first:last))<=.075));
        if day==3,assert(all(S.split(first:last)=="selection"),'All new comparison windows must be held out.');end
        q=repmat(S.q(first,:)',1,variants);v=repmat(S.v(first,:)',1,variants);H=cell(variants,1);
        for a=1:variants,H{a}=preState(models{a},S,first);end
        predicted=NaN(N,2,variants);predicted(1,:,:)=reshape(q,1,2,variants);complete=true(variants,1);
        stopped=NaN(variants,1);residual=NaN(N,2,variants);
        for k=first:last-1
            r=zeros(2,variants);pressure=S.p(k,:)';dt=S.dt(k+1);
            for a=1:variants
                if ~complete(a),q(:,a)=S.q(first,:)';v(:,a)=0;continue;end
                x=[q(:,a);v(:,a);pressure;S.dt(k)];m=models{a};
                [H{a},r(:,a)]=sgCell(m.params,(x-m.normalizer.mean)./m.normalizer.std,H{a},m.cfg);
            end
            [q,v]=sgFlow(q,v,repmat(pressure,1,variants),r,repmat(dt,1,variants),cfg);
            for a=1:variants
                if ~complete(a),continue;end
                if any(~isfinite(q(:,a)))||any(abs(q(:,a))>cfg.coordinateGuard_m)
                    complete(a)=false;stopped(a)=S.t(k+1);continue;
                end
                predicted(k-first+2,:,a)=q(:,a)';residual(k-first+1,:,a)=r(:,a)';
            end
        end
        for a=1:variants
            q0=S.q(first,:)';v0=S.v(first,:)';m=models{a};h=preState(m,S,first);
            [~,r]=sgCell(m.params,([q0;v0;S.p(first,:)';S.dt(first)]-m.normalizer.mean)./m.normalizer.std,h,m.cfg);
            expected=sgFlow(q0,v0,S.p(first,:)',r,S.dt(first+1),m.cfg);
            assert(max(abs(expected'-predicted(2,:,a)))<1e-10);
            good=all(isfinite(A.truth),2)&all(isfinite(predicted(:,:,a)),2);e=1000*(predicted(good,:,a)-A.truth(good,:));
            rm=sqrt(mean(e.^2,1));vr=sqrt(mean(sum(e.^2,2)));if ~complete(a),rm(:)=NaN;vr=NaN;end
            scores=[scores;table(day,string(A.id),A.channel,A.trial,F.results.arm(a),complete(a),nnz(good),N,vr,rm(1),rm(2), ...
                'VariableNames',{'day','runId','channel','trial','variant','complete','scoredSamples','windowSamples', ...
                'lengthVectorRMSE_mm','l2RMSE_mm','l3RMSE_mm'})]; %#ok<AGROW>
        end
        label="nov2_external";if day==3,label="nov3_internal_holdout";end
        selection=[selection;table(day,label,string(A.id),A.channel,A.trial,first,last,S.t(first),S.t(last),A.physicalElapsed(end),nnz(A.missing), ...
            'VariableNames',{'day','role','runId','channel','trial','firstRow','lastRow','startSim_s','endSim_s','durationPhysical_s','missingTruth'})]; %#ok<AGROW>
        tr=struct('day',day,'role',label,'window',A,'predicted',predicted,'residual_N',residual,'complete',complete,'stopAt',stopped);
        trajectories=[trajectories;tr]; %#ok<AGROW>
        fprintf('LONG WINDOW Nov%d %s channel%d trial%d complete %s.\n',day,A.id,A.channel,A.trial,mat2str(complete'));
    end
end
writetable(scores,fullfile(folder,'scores.csv'));writetable(selection,fullfile(folder,'windows.csv'));
save(fullfile(folder,'predictions.mat'),'trajectories','cfg','-v7.3');
renderBundle(camp,trajectories,cfg,F);
for a=1:height(F.results),assert(strcmp(sgHash(F.results.modelPath(a)),F.results.sha256(a)));end
copyfile(mfilename('fullpath')+".m",fullfile(folder,'acMotionCompare.m'));
disp(scores);fprintf('Capacity motion comparison complete; Nov3 windows all heldout.\n');
end
function renderBundle(camp,trajectories,cfg,F)
folder=fullfile(camp,'motion_comparison');scores=sgRead(fullfile(folder,'scores.csv'),'TextType','string');
for day=2:3
    take=trajectories([trajectories.day]==day);
    motionFigure(take(1),folder,sprintf('nov%d_circle_before_after',day),cfg,F.results.arm);
    sineFigure(take(2:4),folder,sprintf('nov%d_sine_before_after',day),cfg,F.results.arm);
    motionFigure(take(5),folder,sprintf('nov%d_helix_before_after',day),cfg,F.results.arm);
end
% Architecture overview rendered separately from saved tables.
end
function h=preState(m,S,k)
first=find(S.mature&S.segment==S.segment(k),1);start=find(S.physicalTime>=S.physicalTime(first)+m.cfg.commonBurnIn_s,1);
h=sgWarm(m.params,S,struct('starts',start,'firstMature',first),k,m.cfg,m.normalizer);
end
function motionFigure(T,folder,name,cfg,ids)
A=T.window;f=figure('Visible','off','Color','w','Position',[60 60 1350 850]);tiledlayout(2,2,'TileSpacing','compact');
colors=[.85 .325 .098;0 .447 .741;.466 .674 .188];styles={'--','-',':'};
for panel=1:3
    nexttile;
    if panel<3,truth=1000*A.truth(:,panel);unit=sprintf('l%d (mm)',panel+1);else,truth=phi(A.truth,cfg);unit='phi (rad)';end
    plot(A.t,truth,'k','LineWidth',1.7);hold on;
    for a=1:3
        if panel<3,z=1000*T.predicted(:,panel,a);else,z=phi(T.predicted(:,:,a),cfg);end
        plot(A.t,z,styles{a},'Color',colors(a,:),'LineWidth',1.3);
    end
    ylabel(unit);xlabel('Elapsed simulation time (s)');grid on;
    if panel==1,legend('Measured','GRU 1x64','GRU 1x128','GRU 2x64','Location','best');end
end
nexttile;z=xyz(A.truth,cfg);plot3(1000*z(:,1),1000*z(:,2),1000*z(:,3),'k','LineWidth',1.7);hold on;
for a=1:3,z=xyz(T.predicted(:,:,a),cfg);plot3(1000*z(:,1),1000*z(:,2),1000*z(:,3),styles{a},'Color',colors(a,:),'LineWidth',1.2);end
axis equal;grid on;xlabel('PCC X (mm)');ylabel('PCC Y (mm)');zlabel('PCC Z (mm)');view(35,25);
sgtitle(sprintf('Nov%d %s trial%d | %s | free prediction, one measured initialization',T.day,A.id,A.trial,strrep(T.role,'_',' ')));finish(f,folder,name);
end
function sineFigure(T,folder,name,cfg,ids)
f=figure('Visible','off','Color','w','Position',[60 60 1450 1000]);tiledlayout(3,3,'TileSpacing','compact');
colors=[.85 .325 .098;0 .447 .741;.466 .674 .188];styles={'--','-',':'};
for ch=1:3
    A=T(ch).window;
    for panel=1:3
        nexttile;
        if panel<3,truth=1000*A.truth(:,panel);unit=sprintf('l%d (mm)',panel+1);else,truth=phi(A.truth,cfg);unit='phi (rad)';end
        plot(A.t,truth,'k','LineWidth',1.7);hold on;
        for a=1:3
            if panel<3,z=1000*T(ch).predicted(:,panel,a);else,z=phi(T(ch).predicted(:,:,a),cfg);end
            plot(A.t,z,styles{a},'Color',colors(a,:),'LineWidth',1.2);
        end
        ylabel(unit);xlabel('Elapsed simulation time (s)');title(sprintf('Channel%d trial%d',ch,A.trial));grid on;
        if ch==1&&panel==1,legend('Measured','GRU 1x64','GRU 1x128','GRU 2x64','Location','best');end
    end
end
sgtitle(sprintf('Nov%d sine | %s | fixed windows, no measured-state resets',T(1).day,strrep(T(1).role,'_',' ')));finish(f,folder,name);
end
function p=phi(q,cfg)
p=hypot((q(:,1)+q(:,2))/cfg.nominal.r,(q(:,2)-q(:,1))/(sqrt(3)*cfg.nominal.r));
end
function z=xyz(q,cfg)
p=phi(q,cfg);angle=atan2(q(:,2)-q(:,1),sqrt(3)*(q(:,1)+q(:,2)));
rho=cfg.nominal.L*(1-cos(p))./p;zz=cfg.nominal.L*sin(p)./p;small=p<1e-8;rho(small)=cfg.nominal.L*p(small)/2;zz(small)=cfg.nominal.L;
z=[rho.*cos(angle) rho.*sin(angle) zz];
end
function finish(f,folder,name)
drawnow;savefig(f,fullfile(folder,[name '.fig']));exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',300);
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');close(f);
end
