function state=pacTuneDemoManualCampaign(outputRoot)
%PACTUNEDEMOMANUALCAMPAIGN Historical replay demo, never a hardware campaign.
% Bounds below are illustrative software-test values, not recommended gains.
here=fileparts(mfilename('fullpath'));addpath(here,fileparts(here));
assert(~isfolder(outputRoot)&&~isfile(outputRoot),'pacTune:Exists','Choose a new demo directory.');
source=fullfile(fileparts(here),'experimentations','adaptCommonP','1rad_10s_gainTesting_v9.mat');
E=pacTuneCampaignConfig();E.trackingWindow_s=[20 80];E.scoringWindow_s=[74 80];E.expectedPeriod_s=2;
trial=pacTuneReadTrial(source,E);
C=pacTuneManualConfig();C.id='historical_v9_offline_demo';C.purpose='offline_demo';C.evaluator=E;
C.baselineTuning=pacTuneBaselineFromTrial(trial,20);
g=C.baselineTuning.adaptationGamma([1 4 6]);C.gammaBounds=[g/2 g*2];
C.evaluationBudget=2;C.confirmationRepeats=2;C.minimumConfirmationDuration_s=100;
root=fullfile(outputRoot,'campaign');pacTuneManualCampaign('create',root,C);
pacTuneManualCampaign('propose',root,'recorded_initial_gammas',g);
m=pacTuneTrialManifest();m.operator='offline historical replay';
m.notes='Existing recording; not newly acquired. Reset/payload conditions are unverified. No operator-run trial was executed.';
pacTuneManualCampaign('import',root,'historical_v9','recorded_initial_gammas',source,m);
state=pacTuneReportCampaign(root,fullfile(outputRoot,'summary'));
assert(~state.success&&~state.configurationReady);
fprintf('OFFLINE DEMO: %s; observations %d; success %d.\n',state.stopReason,numel(state.trialIds),state.success);
end
