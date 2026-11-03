function plotNdiShortGapAudit(root)
%PLOTNDISHORTGAPAUDIT Candidate preview, not filled training or validation truth.
addpath(fullfile(root,'residual_dynamics_learning','force_model','streaming'));
folder=fullfile(root,'residual_dynamics_learning','figures','motion_predictions_20261103_v1','short_gap_audit');
G=sgRead(fullfile(folder,'gap_candidates.csv'),'TextType','string');T=sgRead(fullfile(folder,'reconstruction_summary.csv'),'TextType','string');
ids=unique(G.runId,'stable');counts=zeros(numel(ids),2);errors=counts;p95=counts;
for j=1:numel(ids)
    counts(j,:)=[nnz(G.runId==ids(j)&G.missingSamples==1),nnz(G.runId==ids(j)&G.missingSamples==2)];
    for n=1:2,a=T(T.runId==ids(j)&T.hiddenFrames==n,:);errors(j,n)=a.vectorRMSE_mm;p95(j,n)=a.vectorP95_mm;end
end
f=figure('Visible','off','Color','w','Position',[60 60 1400 860]);tiledlayout(2,2,'TileSpacing','compact');
nexttile;bar(counts,'stacked');xticks(1:numel(ids));xticklabels(ids);xtickangle(20);ylabel('Actual NDI gap episodes');grid on;legend('One frame','Two frames','Location','best');
title('November3: 57 missing samples in 55 episodes');
nexttile;histogram(G.bracketingSpan_s,15);hold on;xline(.12,'r--','Preview cap .12s');xlabel('Valid-to-valid physical span (s)');ylabel('Actual gap episodes');grid on;
nexttile;bar(errors);xticks(1:numel(ids));xticklabels(ids);xtickangle(20);ylabel('Reconstruction vector RMSE (mm)');grid on;legend('Hide one frame','Hide two frames','Location','best');
title('Artificially hidden REAL samples; actual gap error unknown');
nexttile;name="circleSweep";a=G(G.runId==name&flag(G.oneFrameCandidate),:);first=a.firstRow(1);
P=sgRead(fullfile(folder,[char(name) '_preview.csv']));rows=max(1,first-15):min(height(P),first+15);t=P.physicalTime_s(rows)-P.physicalTime_s(first);
plot(t,1000*P.q2_raw_m(rows),'k.-','LineWidth',1.2);hold on;imputed=flag(P.imputed(rows));
plot(t(imputed),1000*P.q2_preview_m(rows(imputed)),'ro','MarkerFaceColor','r');
xlabel('Physical time from actual loss (s)');ylabel('l2 (mm)');grid on;legend('Real measurements','Offline interpolation candidate','Location','best');
title('One-frame preview, not measured ground truth');
sgtitle('Bounded NDI interpolation feasibility; no training/validation data changed');drawnow;
savefig(f,fullfile(folder,'short_gap_audit.fig'));exportgraphics(f,fullfile(folder,'short_gap_audit.png'),'Resolution',300);
exportgraphics(f,fullfile(folder,'short_gap_audit.pdf'),'ContentType','vector');close(f);
copyfile(mfilename('fullpath')+".m",fullfile(folder,'plotNdiShortGapAudit.m'));
end
function x=flag(x)
if islogical(x)||isnumeric(x),x=logical(x);else,x=strcmpi(string(x),'true')|string(x)=="1";end
end
