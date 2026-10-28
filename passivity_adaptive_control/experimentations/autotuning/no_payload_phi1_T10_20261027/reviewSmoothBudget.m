function summary=reviewSmoothBudget()
%REVIEWSMOOTHBUDGET Read completed ledger and plot accuracy/smoothness tradeoff.
folder=fileparts(mfilename('fullpath'));pac=fileparts(fileparts(fileparts(folder)));
addpath(pac,fullfile(pac,'autotuning'));root=fullfile(folder,'campaign_smooth_v2');
state=pacTuneManualCampaign('status',root);
assert(strcmp(state.stopReason,'evaluation_budget_exhausted')&&state.remainingEvaluations==0);
dest=fullfile(root,'completed_search_review');assert(~isfolder(dest),'Choose a fresh review directory.');mkdir(dest);
rows=struct([]);
for k=1:height(state.observations)
    id=char(string(state.observations.trialId(k,:)));
    path=char(string(state.observations.artifactFolder(k,:)));
    a=load(fullfile(path,'evaluation.mat'),'result');r=a.result;
    candidateId=char(string(state.observations.candidateId(k,:)));
    c=load(fullfile(root,'candidates',candidateId,'candidate.mat'),'candidate');g=c.candidate.gamma;
    z=r.scoring;
    rows=[rows;struct('trialId',id,'candidateId',candidateId,'gammaK',g(1),'gammaD',g(2), ...
        'gammaH',g(3),'q2Rmse_mm',z.q2Rmse_mm,'q3Rmse_mm',z.q3Rmse_mm, ...
        'worstCoordinateRmse_mm',max(z.q2Rmse_mm,z.q3Rmse_mm), ...
        'roughness_mPerS2',z.velocityErrorRateRmse_mPerS2, ...
        'velocityRmse_mPerS',z.velocityRmse_mPerS,'feasible',r.success)]; %#ok<AGROW>
end
summary=struct2table(rows);writetable(summary,fullfile(dest,'all_observations.csv'));
fig=figure('Visible','off','Color','w','Position',[50 50 1100 700]);cl=onCleanup(@()close(fig));
hold on;rectangle('Position',[0 0 2 .025],'FaceColor',[.9 .98 .9],'EdgeColor','none');
scatter(summary.worstCoordinateRmse_mm,summary.roughness_mPerS2,65,double(summary.feasible),'filled');
xline(2,'k--','2 mm limit');yline(.025,'k--','Roughness cap');
labels={'Baseline1','Baseline2','03','04','05','06','07','08','Smooth09','Smooth10'};
for k=1:height(summary)
    offset=.00045;if k==2,offset=-.0009;end
    text(summary.worstCoordinateRmse_mm(k)+.025,summary.roughness_mPerS2(k)+offset,labels{k},'FontSize',11);
end
xlim([.9 3.5]);ylim([.012 .037]);grid on
xlabel('Larger of l2/l3 RMSE (mm)');ylabel('Fixed-analysis roughness (m/s^2)');
title('Screening observations: lower is better; green region meets these two limits');
exportgraphics(fig,fullfile(dest,'accuracy_smoothness_tradeoff.png'),'Resolution',200);
save(fullfile(dest,'review.mat'),'state','summary','-v7');disp(summary);
end
