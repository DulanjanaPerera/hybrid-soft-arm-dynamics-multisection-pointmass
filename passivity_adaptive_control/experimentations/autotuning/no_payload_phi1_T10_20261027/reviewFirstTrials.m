function results=reviewFirstTrials(outputFolder)
%REVIEWFIRSTTRIALS Offline audit of the first actual-arm framework recordings.
folder=fileparts(mfilename('fullpath'));pac=fileparts(fileparts(fileparts(folder)));
addpath(pac,fullfile(pac,'autotuning'));
if nargin<1,outputFolder=fullfile(folder,'analysis_review_v1');end
assert(~isfolder(outputFolder),'Choose a fresh output folder.');
names={'smoke_00g_phi1_T10_3cycles_r01','baseline_00g_phi1_T10_12cycles_r01', ...
    'baseline_00g_phi1_T10_12cycles_r02'};
ends=[50 140 140];starts=[40 110 110];results=cell(3,1);
for k=1:3
    cfg=pacTuneCampaignConfig();cfg.trackingWindow_s=[20 ends(k)];
    cfg.scoringWindow_s=[starts(k) ends(k)];cfg.expectedPeriod_s=10;
    trial=pacTuneReadTrial(fullfile(folder,[names{k} '.mat']),cfg);
    result=pacTuneEvaluateTrial(trial,cfg);results{k}=result;
    destination=fullfile(outputFolder,names{k});pacTuneReportTrial(trial,result,destination);
    fid=fopen(fullfile(destination,'actual_settings_and_windows.txt'),'w');cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'SOURCE %s\nSHA256 %s\nPERIOD %.9f CLASS %s\n',trial.file,trial.sha256,result.period_s,result.classification);
    for f=fieldnames(result.initialTuning)'
        fprintf(fid,'%s %s\n',f{1},mat2str(result.initialTuning.(f{1}),9));
    end
    fprintf(fid,'INITIAL estimates%s history%s\n',mat2str(trial.signals.parameterEstimate.x(1,:),7),mat2str(trial.signals.historyState.x(1,:),7));
    t=trial.signals.parameterEstimate.t;circle=t>=20&t<ends(k);
    parameters=trial.signals.parameterEstimate.x;
    fprintf(fid,'END_TRACKING estimates%s\n',mat2str(parameters(find(circle,1,'last'),:),7));
    fprintf(fid,'bounds lowerFrac%s upperFrac%s\n',mat2str(result.lowerBoundFraction,5),mat2str(result.upperBoundFraction,5));
    fprintf(fid,'ratecap%s gateMismatch%d\n',mat2str(result.rates.capFraction,5),result.gateActiveMismatch);
    fprintf(fid,'pressure minmax%s residualPeak%g replay%s\n',mat2str(result.pressureRange_bar,5),result.nominalResidualPeak_N,jsonencode(result.allocationReplay));
    fprintf(fid,'timing maxControllerStep%.7g\n',result.maximumObservedStep_s);
    for win=[20 40 50 80 110;40 50 80 110 140]
        if win(2)>ends(k),continue;end
        c=cfg;c.scoringWindow_s=win';r=pacTuneEvaluateTrial(trial,c);
        fprintf(fid,'WINDOW %.0f-%.0f %s\n',win(1),win(2),jsonencode(r.scoring));
    end
    fprintf(fid,'FAULT whole run ');status=trial.signals.status;
    fprintf(fid,'%d first %s\n',nnz(status.x==5),mat2str(status.t(find(status.x==5,1))));
    clear cleanup
    fprintf('%s: coverage %.3f%% active %.3f%% limited %.3f%% changes%d faults%d; late theta%.5f phi%.5f length[%.3f %.3f]mm\n', ...
        names{k},100*result.tracking.coverage,100*result.tracking.activeFraction,100*result.tracking.limitedFraction, ...
        height(result.tuningEvents),result.tracking.faultSamples,result.scoring.thetaMae_rad, ...
        result.scoring.phiMae_rad,result.scoring.q2Rmse_mm,result.scoring.q3Rmse_mm);
end
save(fullfile(outputFolder,'comparison.mat'),'results','-v7');
end
