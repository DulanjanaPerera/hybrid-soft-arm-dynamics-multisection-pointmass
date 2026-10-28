function results=pacTuneReplayHistorical(outputRoot)
%PACTUNEREPLAYHISTORICAL Evaluate v8/v9/v10 without running any model.
here=fileparts(mfilename('fullpath'));addpath(fileparts(here));
data=fullfile(fileparts(here),'experimentations','adaptCommonP');
assert(~isfolder(outputRoot),'pacTune:OutputExists','Choose a new output directory.');
periods=[10 2 4];ends=[320 80 140];scores=[290 74 94;320 80 102]';
results=cell(3,1);
for k=1:3
    cfg=pacTuneCampaignConfig();cfg.trackingWindow_s=[20 ends(k)];
    cfg.scoringWindow_s=scores(k,:);cfg.expectedPeriod_s=periods(k);
    name=sprintf('1rad_10s_gainTesting_v%d',k+7);
    trial=pacTuneReadTrial(fullfile(data,[name '.mat']),cfg);
    result=pacTuneEvaluateTrial(trial,cfg);results{k}=result;
    pacTuneReportTrial(trial,result,fullfile(outputRoot,name));
    fprintf('v%d period %.4g, active %.4f%%, limited %.4f%%, %s\n', ...
        k+7,result.period_s,100*result.tracking.activeFraction, ...
        100*result.tracking.limitedFraction,result.classification);
end
save(fullfile(outputRoot,'historical_results.mat'),'results','-v7');
end
