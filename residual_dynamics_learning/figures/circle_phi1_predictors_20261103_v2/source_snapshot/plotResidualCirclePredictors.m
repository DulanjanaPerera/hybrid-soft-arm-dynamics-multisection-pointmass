function plotResidualCirclePredictors(root,folder)
%PLOTRESIDUALCIRCLEPREDICTORS Requested phi=1 section, five frozen predictors.
% Offline recorded-input predictions, not new closed-loop arm experiments.
assert(~isfolder(folder),'Use a fresh output folder.');
addpath(root,fullfile(root,'residual_dynamics_learning'));
evaluation=fullfile(root,'residual_dynamics_learning','validations','fresh_20261102_v2');
screen=fullfile(root,'residual_dynamics_learning','benchmarks','screen_20261101_v2');
dataset=fullfile(root,'residual_dynamics_learning','datasets','two_day_v2');
export=fullfile(root,'RL_data_collection','20261102','exports','circleSweep');
z=load(fullfile(evaluation,'circleSweep_sequences.mat'));S=z.S;cfg=z.cfg;
saved=load(fullfile(evaluation,'circleSweep_one_step.mat'));
f=load(fullfile(screen,'forward','frozen_selection.mat'));
kin=load(fullfile(screen,'ablations','forward_frozen.mat'));
n=load(fullfile(dataset,'normalizer.mat'));n=n.normalizer;
O=readtable(fullfile(export,'observations.csv'),'VariableNamingRule','preserve');
assert(isequal(S.t,O.time_sim_s));
selected=O.phaseId==5&abs(O.phi_reference_rad-1)<1e-10;
trials=unique(O.trialId(selected));assert(isscalar(trials));
cycles=unique(O.cycleId(selected));assert(numel(cycles)==5);
% First chronological complete circle, chosen without inspecting model errors.
rows=find(selected&O.cycleId==cycles(1));start=rows(1);last=rows(end)+1;
assert(last<=height(O)&&abs(O.phi_reference_rad(last)-1)<1e-10);
check=start-cfg.historyFrames:last;
assert(all(S.mature(check))&&all(S.segment(check)==S.segment(start)), ...
    'The first complete circle crosses a gap; do not silently bridge it.');
protocol=struct('requestedPhi_rad',1,'trialId',trials,'cycleId',cycles(1), ...
    'startSim_s',S.t(start),'endSim_s',S.t(last), ...
    'durationPhysical_s',S.physicalTime(last)-S.physicalTime(start), ...
    'initialHistoryFrames',cfg.historyFrames,'measuredStateResets',0, ...
    'boundMode',0,'radius_m',.013,'moduleLength_m',.17411, ...
    'interpretation','Recorded-input forward prediction; no controller/hardware experiment');
mkdir(folder);save(fullfile(folder,'protocol.mat'),'protocol','cfg');
locked={fullfile(screen,'forward','frozen_selection.mat'), ...
    fullfile(screen,'ablations','forward_frozen.mat'),fullfile(dataset,'normalizer.mat'), ...
    fullfile(evaluation,'circleSweep_sequences.mat'),fullfile(evaluation,'circleSweep_one_step.mat'), ...
    fullfile(export,'observations.csv'),fullfile(root,'RL_data_collection','20261102','circleSweep.mat')};
hashes=cellfun(@filehash,locked,'UniformOutput',false);
writetable(table(string(locked(:)),string(hashes(:)),'VariableNames',{'source','sha256'}),fullfile(folder,'sources.csv'));
ids={'nominal_legacy','constant_velocity','ridge','kinematic_ridge','gru'};
labels={'Nominal model','Constant velocity','Nominal + ridge','Kinematic + ridge','Nominal + GRU'};
N=last-start+1;Qpred=NaN(N,2,numel(ids));completed=false(numel(ids),1);stopAt=NaN(numel(ids),1);
for m=1:numel(ids)
    Q=S.q(start-cfg.historyFrames+1:start,:)';V=S.v(start-cfg.historyFrames+1:start,:)';
    Qpred(1,:,m)=Q(:,end)';
    for row=start:last-1
        dt=S.dt(row+1);history=row-cfg.historyFrames+1:row;
        if ismember(ids{m},{'constant_velocity','kinematic_ridge'})
            base=Q(:,end)+V(:,end)*dt;
        else
            state=residualFrozenStep(Q(:,end),V(:,end),S.p(row,:)',dt,0);base=state(1:2);
        end
        next=base;
        if ismember(ids{m},{'ridge','kinematic_ridge','gru'})
            values=[Q' V' S.p(history,:)];
            D=struct('Seq',single([(values-n.mean)./n.std S.dt(history)]'), ...
                'Context',[dt base']);
            if strcmp(ids{m},'kinematic_ridge'),model=kin.bank(1);
            else,model=f.bank(strcmp({f.bank.name},ids{m}));end
            next=base+residualPredict(model,D)'/1000;
        end
        if any(~isfinite(next))||any(abs(next)>.05),stopAt(m)=S.t(row+1);break;end
        velocity=V(:,end)+dt/(cfg.velocityTau_s+dt)*((next-Q(:,end))/dt-V(:,end));
        Q=[Q(:,2:end) next];V=[V(:,2:end) velocity];Qpred(row-start+2,:,m)=next';
    end
    completed(m)=all(isfinite(Qpred(:,:,m)),'all');
    fprintf('Circle prediction %s: complete=%d\n',labels{m},completed(m));
end
% Use exactly the same PCC map on measured inferred q and predicted q.
% This XY is an orientation/length-derived visualization, not raw sensor XYZ.
actual=(start:last)';[truthXY,truthTheta,truthPhi]=geometry(S.q(actual,:),protocol);
referenceTheta=O.theta_reference_rad(actual);referencePhi=O.phi_reference_rad(actual);
qref=protocol.radius_m*referencePhi.*[cos(referenceTheta)/2-sqrt(3)*sin(referenceTheta)/2, ...
    cos(referenceTheta)/2+sqrt(3)*sin(referenceTheta)/2];
[referenceXY,~,~]=geometry(qref,protocol);
phiError=max(abs(truthPhi-O.phi_measured_rad(actual)));
thetaError=max(abs(atan2(sin(truthTheta-O.theta_measured_rad(actual)),cos(truthTheta-O.theta_measured_rad(actual)))));
assert(phiError<1e-10&&thetaError<1e-10,'Geometry convention mismatch.');
stats=table;series=table;XY=NaN(N,2,numel(ids));TH=NaN(N,numel(ids));PH=TH;
elapsed=S.physicalTime(actual)-S.physicalTime(start);relativeTime=S.t(actual)-S.t(start);
for m=1:numel(ids)
    [XY(:,:,m),TH(:,m),PH(:,m)]=geometry(Qpred(:,:,m),protocol);
    good=all(isfinite(Qpred(:,:,m)),2);err=1000*(Qpred(good,:,m)-S.q(actual(good),:));
    eXY=1000*(XY(good,:,m)-truthXY(good,:));
    stats=[stats;table(string(labels{m}),completed(m),stopAt(m),nnz(good), ...
        sqrt(mean(sum(err.^2,2))),sqrt(mean(sum(eXY.^2,2))), ...
        mean(abs(PH(good,m)-truthPhi(good))), ...
        'VariableNames',{'predictor','completedCircle','guardStopSim_s','samples', ...
        'lengthVectorRMSE_mm','PccXYVectorRMSE_mm','phiMAE_rad'})]; %#ok<AGROW>
    series=[series;table(repmat(string(labels{m}),N,1),actual,relativeTime,elapsed, ...
        S.q(actual,1),S.q(actual,2),Qpred(:,1,m),Qpred(:,2,m), ...
        truthTheta,truthPhi,TH(:,m),PH(:,m),referenceTheta,referencePhi, ...
        truthXY(:,1),truthXY(:,2),XY(:,1,m),XY(:,2,m),referenceXY(:,1),referenceXY(:,2), ...
        'VariableNames',{'predictor','sourceRow','elapsedSim_s','elapsedPhysical_s', ...
        'measuredL2_m','measuredL3_m','predictedL2_m','predictedL3_m', ...
        'measuredTheta_rad','measuredPhi_rad','predictedTheta_rad','predictedPhi_rad', ...
        'referenceTheta_rad','referencePhi_rad','measuredPccX_m','measuredPccY_m', ...
        'predictedPccX_m','predictedPccY_m','referencePccX_m','referencePccY_m'})]; %#ok<AGROW>
end
writetable(stats,fullfile(folder,'full_circle_scores.csv'));writetable(series,fullfile(folder,'full_circle_predictions.csv'));
colors=[0 .447 .741;.85 .325 .098;.466 .674 .188;.494 .184 .556;.301 .745 .933];
fig=figure('Visible','off','Color','w','Position',[60 60 1400 860]);tiledlayout(2,3,'TileSpacing','compact');
limit=1.12*1000*max(abs([referenceXY;truthXY;reshape(permute(XY,[1 3 2]),[],2)]),[],'all');
for m=1:numel(ids)
    nexttile;plot(1000*referenceXY(:,1),1000*referenceXY(:,2),'r--','LineWidth',1.3);hold on;
    plot(1000*truthXY(:,1),1000*truthXY(:,2),'k-','LineWidth',1.8);
    plot(1000*XY(:,1,m),1000*XY(:,2,m),'Color',colors(m,:),'LineWidth',1.5);
    plot(1000*truthXY(1,1),1000*truthXY(1,2),'ko','MarkerFaceColor','k');
    axis equal;xlim([-limit limit]);ylim([-limit limit]);grid on;xlabel('PCC X (mm)');ylabel('PCC Y (mm)');
    title(labels{m});
    if ~completed(m),text(.04,.06,'Stopped at numerical guard','Units','normalized','Color',[.8 0 0]);end
end
nexttile;hold on;plot(NaN,NaN,'r--','LineWidth',1.3);plot(NaN,NaN,'k-','LineWidth',1.8);plot(NaN,NaN,'b-','LineWidth',1.5);
axis off;legend('Requested circle','Measured response (PCC)','Predictor output','Location','northwest');
text(.02,.52,sprintf('Requested phi = 1 rad\nCircle 1: %.0f-%.0f s simulation\n%.2f s physical duration\nMeasured initialization only\nNo state resets during circle', ...
    S.t(start),S.t(last),protocol.durationPhysical_s),'Units','normalized','FontSize',12);
sgtitle('Full-circle forward predictions under recorded pressures; no feedback control');finish(fig,folder,'full_circle_xy');
fig=figure('Visible','off','Color','w','Position',[60 60 1350 800]);tiledlayout(2,1,'TileSpacing','compact');
nexttile;plot(relativeTime,referencePhi,'r--','LineWidth',1.5);hold on;plot(relativeTime,truthPhi,'k-','LineWidth',2);
for m=1:numel(ids),plot(relativeTime,PH(:,m),'Color',colors(m,:),'LineWidth',1.2);end
grid on;ylabel('phi (rad)');title('One full circle: initialized once, then predicted motion only');
legend([{'Requested','Measured'},labels],'Location','eastoutside');
nexttile;plot(relativeTime,referenceTheta-referenceTheta(1),'r--','LineWidth',1.5);hold on;
plot(relativeTime,lift(truthTheta,referenceTheta(1))-referenceTheta(1),'k-','LineWidth',2);
for m=1:numel(ids),plot(relativeTime,lift(TH(:,m),referenceTheta(1))-referenceTheta(1),'Color',colors(m,:),'LineWidth',1.2);end
grid on;ylabel('theta relative to reference start (rad)');xlabel('Elapsed simulation time (s)');finish(fig,folder,'full_circle_angles');
% Existing teacher-forced one-step predictions across all five circles.
target=saved.idx+1;choose=selected(target);target=target(choose);t=S.t(target)-S.t(find(selected,1));
oneStep=table;fig=figure('Visible','off','Color','w','Position',[60 60 1350 800]);tiledlayout(2,1,'TileSpacing','compact');
ax1=nexttile;plot(t,O.phi_reference_rad(target),'r--','LineWidth',1.5);hold on;plot(t,O.phi_measured_rad(target),'k-','LineWidth',2);ylabel('phi (rad)');grid on;
ax2=nexttile;origin=O.theta_reference_rad(target(1));plot(t,O.theta_reference_rad(target)-origin,'r--','LineWidth',1.5);hold on;
plot(t,lift(O.theta_measured_rad(target),origin)-origin,'k-','LineWidth',2);ylabel('theta relative to reference start (rad)');xlabel('Elapsed simulation time (s)');grid on;
for m=1:numel(ids)
    index=find(strcmp(saved.names,ids{m}));assert(isscalar(index));q=saved.preds{index}(choose,:);
    [xy,theta,phi]=geometry(q,protocol);plot(ax1,t,phi,'Color',colors(m,:),'LineWidth',1.2);
    plot(ax2,t,lift(theta,origin)-origin,'Color',colors(m,:),'LineWidth',1.2);
    oneStep=[oneStep;table(repmat(string(labels{m}),numel(target),1),target,t,q(:,1),q(:,2),theta,phi,xy(:,1),xy(:,2), ...
        'VariableNames',{'predictor','sourceRow','elapsedSim_s','predictedL2_m','predictedL3_m', ...
        'predictedTheta_rad','predictedPhi_rad','predictedPccX_m','predictedPccY_m'})]; %#ok<AGROW>
end
title(ax1,'One-step predictions: measured history supplied again at every sample');
legend(ax1,[{'Requested','Measured'},labels],'Location','eastoutside');finish(fig,folder,'one_step_five_circles_angles');
writetable(oneStep,fullfile(folder,'one_step_predictions.csv'));
save(fullfile(folder,'predictions.mat'),'protocol','ids','labels','Qpred','actual','stats','selected');
mkdir(fullfile(folder,'source_snapshot'));
for name={'plotResidualCirclePredictors.m','residualFrozenStep.m','residualPredict.m','ndiSupportFnc.m'}
    copyfile(which(name{1}),fullfile(folder,'source_snapshot',name{1}));
end
assert(isequal(hashes,cellfun(@filehash,locked,'UniformOutput',false)));
disp(stats);fprintf('PASS: logged phi1 selection, complete contiguous first circle, angle/length roundtrip, unchanged frozen files, no fitting/hardware.\n');
end

function [xy,theta,phi]=geometry(q,p)
u=(q(:,1)+q(:,2))/p.radius_m;v=(q(:,2)-q(:,1))/(sqrt(3)*p.radius_m);
phi=hypot(u,v);theta=atan2(v,u);rho=zeros(size(phi));nonzero=phi>1e-8;
rho(nonzero)=p.moduleLength_m*(1-cos(phi(nonzero)))./phi(nonzero);
rho(~nonzero)=p.moduleLength_m*phi(~nonzero)/2;
xy=rho.*[cos(theta) sin(theta)];
end
function theta=lift(theta,reference)
good=isfinite(theta);theta(good)=unwrap(theta(good));first=find(good,1);
theta(good)=theta(good)+2*pi*round((reference-theta(first))/(2*pi));
end
function finish(fig,folder,name)
savefig(fig,fullfile(folder,[name '.fig']));exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',300);
exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector');close(fig);
end
function h=filehash(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
