function sgFigures(root)
base=fullfile(root,'residual_dynamics_learning','force_model','streaming');campaign=fullfile(base,'campaign_20261103_v2');
folder=fullfile(base,'evaluation_20261103_v1');S=sgRead(fullfile(campaign,'screen.csv'),'TextType','string');
F=sgRead(fullfile(campaign,'finalists.csv'),'TextType','string');
P=sgRead(fullfile(folder,'short_paths.csv'),'TextType','string');
f=figure('Visible','off','Color','w','Position',[80 80 1350 800]);tiledlayout(2,2,'TileSpacing','compact');
nexttile;bar(categorical(S.id),S.oldScore_mm);ylabel('Old holdout weighted score (mm)');grid on;
title('Six-model screen, seed113; lower is better');
nexttile;hold on;
for j=1:height(S)
    H=sgRead(fullfile(S.folder(j),'training_history.csv'));plot(H.updates,H.oldScore_mm,'-o','DisplayName',S.id(j));
end
xlabel('Optimizer updates');ylabel('Old holdout weighted score (mm)');legend('Location','northeast');grid on;title('Equal update budgets, old data only');
ids=unique(F.id,'stable');nexttile;hold on;
for j=1:numel(ids)
    a=F.oldScore_mm(F.id==ids(j));errorbar(j,mean(a),std(a),'o','LineWidth',1.5);
end
xlim([.5 numel(ids)+.5]);xticks(1:numel(ids));xticklabels(ids);ylabel('Old holdout score mean +/- SD (mm)');grid on;title('Three seeds per finalist');
nexttile;hold on;hz=[.1 .25 .5 1];
for j=1:numel(ids)
    scores=zeros(3,4);seeds=unique(F.seed);
    for k=1:3
        for h=1:4,scores(k,h)=mean(P.pathVectorRMSE_mm(P.id==ids(j)&P.seed==seeds(k)&P.horizon_s==hz(h)));end
    end
    errorbar(hz,mean(scores),std(scores),'-o','DisplayName',ids(j));
end
baseline=sgRead(fullfile(folder,'short_baselines.csv'),'TextType','string');
for name=["nominal","constant_velocity"]
    y=zeros(1,4);for h=1:4,y(h)=mean(baseline.pathVectorRMSE_mm(baseline.predictor==name&baseline.horizon_s==hz(h)));end
    plot(hz,y,'--','DisplayName',strrep(name,'_',' '));
end
set(gca,'YScale','log');
xlabel('Maximum physical horizon (s)');ylabel('Nov2 macro path vector RMSE (mm)');legend('Location','best');grid on;title('Frozen finalist evaluation, mean +/- SD across seeds');
sgtitle('Generalized-force GRU: streaming measurements and separate free prediction');bundle(f,folder,'streaming_comparison');
% Predeclared representative: best MEAN OLD holdout configuration, fixed seed113.
G=groupsummary(F,'id','mean','oldScore_mm');[~,j]=min(G.mean_oldScore_mm);id=G.id(j);seed=113;
file=fullfile(folder,sprintf('%s_seed%d_run4_stream.mat',id,seed));a=load(file);
raw=load(fullfile(root,'residual_dynamics_learning','validations','fresh_20261102_v2','circleSweep_sequences.mat'));S=raw.S;
first=1;elapsed=S.physicalTime(a.anchors)';keep=find(elapsed-elapsed(first)<=10 & a.segmentIds==a.segmentIds(first));
t=S.physicalTime(a.anchors(keep)+1)-S.physicalTime(a.anchors(first));
f=figure('Visible','off','Color','w','Position',[80 80 1300 850]);tiledlayout(3,1,'TileSpacing','compact');
for d=1:2
    nexttile;plot(t,1000*a.truth(d,keep),'k',t,1000*a.Q(d,keep,1),'-',t,1000*a.Q(d,keep,3),'--');
    ylabel(sprintf('l%d (mm)',d+1));grid on;if d==1,legend('Observed','Force GRU one-step','Constant velocity one-step','Location','best');end
end
nexttile;plot(t,a.R(:,keep)');ylabel('Residual force (N)');xlabel('Physical elapsed time (s)');grid on;legend('r2','r3');
sgtitle(sprintf('Circle recording, first eligible 10 s: %s seed113; one-step prediction, not tracking',id));bundle(f,folder,'circle_streaming_example');
% First chronologically selected circle free window, SAME representative.
a=load(fullfile(folder,sprintf('%s_seed%d_short.mat',id,seed)));d=a.detail;i=find(d.index.run==4,1);mask=logical(reshape(d.batch.auxMask(1,i,:),1,[]));
t=cumsum(reshape(d.batch.auxDt(1,i,mask),1,[]));q=reshape(d.Q(:,i,mask),2,[]);truth=reshape(d.batch.auxTruth(:,i,mask),2,[]);
f=figure('Visible','off','Color','w','Position',[80 80 1100 700]);tiledlayout(2,1);
for k=1:2
    nexttile;plot(t,1000*truth(k,:),'k',t,1000*q(k,:),'-');ylabel(sprintf('l%d (mm)',k+1));grid on;legend('Observed','Free integrated force GRU','Location','best');
end
xlabel('Physical elapsed time (s)');sgtitle(sprintf('Circle recording, first eligible window: %s seed113; free prediction, no resets',id));bundle(f,folder,'circle_short_free_example');
copyfile(mfilename('fullpath')+".m",fullfile(folder,'sgFigures.m'));
end
function bundle(f,folder,name)
drawnow;savefig(f,fullfile(folder,[name '.fig']));exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',300);
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');close(f);
end
