function validatePacTuneManualCampaign()
%VALIDATEPACTUNEMANUALCAMPAIGN Synthetic manual-loop test; no hardware.
here=fileparts(mfilename('fullpath'));addpath(here,fileparts(here));
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
source=fullfile(fileparts(here),'experimentations','adaptCommonP','1rad_10s_gainTesting_v8.mat');
cfg=pacTuneCampaignConfig();trial=pacTuneReadTrial(source,cfg);
C=pacTuneManualConfig();C.id='synthetic_test';C.purpose='manual';
C.baselineTuning=pacTuneBaselineFromTrial(trial,20);
C.gammaBounds=[.01 2;1e-6 .1;.01 10];C.evaluationBudget=6;
C.confirmationRepeats=2;C.minimumConfirmationDuration_s=100;
C.conditions=struct('payloadMass_kg',0,'payloadAttachedBeforeTrim','yes', ...
    'calibrationId','synthetic','modelId','synthetic','resetProcedure','synthetic', ...
    'preconditioning','synthetic','phi_rad',1,'direction',1,'codeRevision','synthetic');
C.evaluator.expectedPeriod_s=10;
for f=fieldnames(C.evaluator.targets)',C.evaluator.targets.(f{1})=1e6;end
C.evaluator.targets.limitedFraction=1;
C.evaluator.trackingWindow_s=[20 40];C.evaluator.scoringWindow_s=[30 40];
C.confirmationEvaluator=C.evaluator;
C.confirmationEvaluator.trackingWindow_s=[20 320];C.confirmationEvaluator.scoringWindow_s=[290 320];
root=fullfile(folder,'campaign');state=pacTuneManualCampaign('create',root,C);
assert(state.configurationReady&&isempty(state.candidateIds));
g=C.baselineTuning.adaptationGamma([1 4 6]);
candidate=pacTuneManualCampaign('propose',root,'candidate01',g);
assert(isequal(candidate.tuning.adaptationGamma,C.baselineTuning.adaptationGamma));
expectError(@()pacTuneManualCampaign('propose',root,'outside',[20;1;20]),'pacTune:Bounds');
expectError(@()pacTuneManualCampaign('propose',root,'../escape',g),'pacTune:Id');
expectError(@()pacTuneManualCampaign('propose',root,'candidate01',g),'pacTune:Duplicate');
% Synthetic fixed-gamma fixture: copy the Dataset; original recording untouched.
S=load(source,'out');ds=S.out.logsout;ix=find(strcmp(ds.getElementNames,'adaptationGamma'));
e=ds{ix};e.Values=timeseries(repmat(C.baselineTuning.adaptationGamma',numel(e.Values.Time),1),e.Values.Time);
ds=ds.setElement(ix,e);out=struct('logsout',ds);
fixture=fullfile(folder,'synthetic.mat');marker='synthetic fixture 1';save(fixture,'out','marker');
m=pacTuneTrialManifest();m.conditions=C.conditions;m.resetVerified=true;m.operator='synthetic test';
o=pacTuneManualCampaign('import',root,'trial01','candidate01',fixture,m);
assert(o.eligible&&strcmp(o.classification,'valid_feasible'));
expectError(@()pacTuneManualCampaign('import',root,'trial01','candidate01',fixture,m),'pacTune:Duplicate');
expectError(@()pacTuneManualCampaign('import',root,'renamed','candidate01',fixture,m),'pacTune:Duplicate');
state=pacTuneManualCampaign('status',root);
assert(strcmp(state.bestScreenedCandidate,'candidate01')&&~state.success);
% Distinct synthetic files model independent repeated observations; labelled only.
m.kind='confirmation';
for k=2:3
    marker=sprintf('synthetic fixture %d',k);file=fullfile(folder,sprintf('synthetic%d.mat',k));save(file,'out','marker');
    pacTuneManualCampaign('import',root,sprintf('trial%02d',k),'candidate01',file,m);
end
state=pacTuneManualCampaign('status',root);
assert(state.success&&strcmp(state.bestConfirmedCandidate,'candidate01')&&state.remainingEvaluations==3);
% Metadata mismatch cannot win despite permissive numeric metrics.
marker='synthetic mismatch';file=fullfile(folder,'mismatch.mat');save(file,'out','marker');
m.conditions.payloadMass_kg=.1;
o=pacTuneManualCampaign('import',root,'mismatch','candidate01',file,m);
assert(~o.eligible&&strcmp(o.classification,'protocol_deviation'));
state=pacTuneManualCampaign('status',root);
assert(~state.success&&isempty(state.bestConfirmedCandidate),'Failed confirmation must withdraw confirmation status.');
assert(isempty(state.bestScreenedCandidate)&&strcmp(state.stopReason,'awaiting_candidate'));
% Logged gamma changes rejected against fixed candidate even with matching manifest.
m.conditions=C.conditions;m.kind='confirmation';
o=pacTuneManualCampaign('import',root,'exploratory','candidate01',source,m);
assert(~o.eligible&&strcmp(o.evaluatorClassification,'protocol_deviation'));
% Restart, locking and immutable configuration detection.
clear pacTuneManualCampaign
state=pacTuneManualCampaign('status',root);assert(height(state.observations)==5);
lock=java.io.File(fullfile(root,'.writer.lock'));assert(lock.createNewFile());
expectError(@()pacTuneManualCampaign('propose',root,'locked',g),'pacTune:Locked');lock.delete();
marker='synthetic budget fixture';file=fullfile(folder,'budget.mat');save(file,'out','marker');
m.resetVerified=false;pacTuneManualCampaign('import',root,'trial06','candidate01',file,m);
state=pacTuneManualCampaign('status',root);
assert(strcmp(state.stopReason,'evaluation_budget_exhausted')&&~state.success);
expectError(@()pacTuneManualCampaign('propose',root,'after_budget',g),'pacTune:Budget');
pacTuneReportCampaign(root,fullfile(folder,'snapshot'));
% Objective search can continue after a feasible screen in a NEW config only.
searchC=C;searchC.id='continue_search_test';searchC.optimizer.continueSearchAfterFeasible=true;
searchRoot=fullfile(folder,'search_campaign');pacTuneManualCampaign('create',searchRoot,searchC);
pacTuneManualCampaign('propose',searchRoot,'baseline',g);
sm=pacTuneTrialManifest();sm.conditions=searchC.conditions;sm.operator='synthetic';sm.resetVerified=true;
pacTuneManualCampaign('import',searchRoot,'screen','baseline',fixture,sm);
ss=pacTuneManualCampaign('status',searchRoot);assert(~isempty(ss.bestScreenedCandidate));
proposal=pacTuneSuggestNext(searchRoot);assert(strcmp(proposal.action,'propose'));
pacTuneQueueSuggestion(searchRoot,'next',proposal);
ss=pacTuneManualCampaign('status',searchRoot);assert(strcmp(ss.stopReason,'awaiting_observation'));
C.conditions.payloadMass_kg=.5;save(fullfile(root,'campaign.mat'),'C');
expectError(@()pacTuneManualCampaign('status',root),'pacTune:ChangedConfig');
fprintf('PASS: manual create/propose/import/status, bounds, duplicates, repeats, metadata, restart, locks, immutable config.\n');
end

function expectError(fn,id)
failed=false;
try,fn();catch err,failed=strcmp(err.identifier,id);if ~failed,rethrow(err);end,end
assert(failed,'Expected error %s.',id);
end
