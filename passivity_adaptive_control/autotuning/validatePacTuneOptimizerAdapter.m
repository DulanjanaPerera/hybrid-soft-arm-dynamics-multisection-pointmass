function validatePacTuneOptimizerAdapter()
%VALIDATEPACTUNEOPTIMIZERADAPTER Temporary test-only campaign on saved data.
here=fileparts(mfilename('fullpath'));addpath(here,fileparts(here));
data=fullfile(fileparts(here),'experimentations','adaptCommonP');
source=fullfile(data,'1rad_10s_gainTesting_v8.mat');
trial=pacTuneReadTrial(source);C=pacTuneManualConfig();C.id='adapter_test';C.purpose='manual';
C.baselineTuning=pacTuneBaselineFromTrial(trial,20);g=C.baselineTuning.adaptationGamma([1 4 6]);
C.gammaBounds=[g/2 g*2];C.evaluationBudget=6;C.confirmationRepeats=2;
C.minimumConfirmationDuration_s=20;
C.evaluator.trackingWindow_s=[20 40];C.evaluator.scoringWindow_s=[30 40];C.evaluator.expectedPeriod_s=10;
for f=fieldnames(C.evaluator.targets)',C.evaluator.targets.(f{1})=1e6;end
C.evaluator.targets.qVectorRmse_mm=.01; % Test-only infeasibility, not a real target.
C.conditions=struct('payloadMass_kg',0,'payloadAttachedBeforeTrim','yes', ...
    'calibrationId','synthetic','modelId','synthetic','resetProcedure','synthetic', ...
    'preconditioning','synthetic','phi_rad',1,'direction',1,'codeRevision','synthetic');
C.optimizer.initialDesignSize=1;C.optimizer.poolSize=400;
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
root=fullfile(folder,'campaign');pacTuneManualCampaign('create',root,C);
p=pacTuneSuggestNext(root);assert(strcmp(p.reason,'baseline')&&isequal(p.gamma,g));
candidate=pacTuneQueueSuggestion(root,'baseline',p);assert(~isempty(candidate.optimizer));
next=pacTuneSuggestNext(root);assert(strcmp(next.action,'observe'));
failed=false;
try,pacTuneQueueSuggestion(root,'stale',p);catch err,failed=strcmp(err.identifier,'pacTune:StaleSuggestion');end
assert(failed,'Stale optimizer snapshots must not queue.');
m=pacTuneTrialManifest();m.conditions=C.conditions;m.operator='test only';m.resetVerified=true;
obs=pacTuneManualCampaign('import',root,'screen','baseline',source,m);
assert(strcmp(obs.classification,'valid_infeasible'));
q=pacTuneSuggestNext(root);assert(strcmp(q.action,'propose'));
assert(q.training.usableObjectives==1&&q.training.knownConstraints==1);
% Confirmation data must not train or change the GP suggestion.
m.kind='confirmation';
pacTuneManualCampaign('import',root,'held_out','baseline',fullfile(data,'1rad_10s_gainTesting_v9.mat'),m);
r=pacTuneSuggestNext(root);
assert(r.excluded.confirmation==1&&r.training.observations==1&&isequal(r.gamma,q.gamma));
% Read-only replay from the same ledger is deterministic.
clear pacTuneSuggestNext
s=pacTuneSuggestNext(root);assert(isequaln(r,s));
% A fault must not hide live-tuning/protocol changes from the training filter.
source10=fullfile(data,'1rad_10s_gainTesting_v10.mat');t10=pacTuneReadTrial(source10);
C2=C;C2.id='fault_mapping_test';C2.evaluator.trackingWindow_s=[20 140];
C2.evaluator.scoringWindow_s=[94 102];C2.evaluator.expectedPeriod_s=4;
C2.baselineTuning=pacTuneBaselineFromTrial(t10,20);g10=C2.baselineTuning.adaptationGamma([1 4 6]);
C2.gammaBounds=[g10/2 g10*2];root2=fullfile(folder,'fault_campaign');
pacTuneManualCampaign('create',root2,C2);p=pacTuneSuggestNext(root2);
pacTuneQueueSuggestion(root2,'fault_baseline',p);m.kind='screening';m.conditions=C2.conditions;
pacTuneManualCampaign('import',root2,'changed_gamma','fault_baseline',source10,m);
p=pacTuneSuggestNext(root2);
assert(p.training.knownConstraints==0&&p.training.usableObjectives==0&&p.excluded.protocolOrData==1);
% Explicit synthetic fixed-gamma copy tests comparable fault mapping only.
S=load(source10,'out');ds=S.out.logsout;ix=find(strcmp(ds.getElementNames,'adaptationGamma'));
e=ds{ix};e.Values=timeseries(repmat(C2.baselineTuning.adaptationGamma',numel(e.Values.Time),1),e.Values.Time);
ds=ds.setElement(ix,e);out=struct('logsout',ds);fixture=fullfile(folder,'synthetic_fault.mat');save(fixture,'out');
pacTuneManualCampaign('import',root2,'comparable_fault','fault_baseline',fixture,m);
p=pacTuneSuggestNext(root2);
assert(p.training.knownConstraints==1&&p.training.usableObjectives==0,'A comparable fault supplies a constraint label, never a partial score.');
demo=fullfile(here,'reports','manual_demo_v1','campaign');
if isfolder(demo)
    p=pacTuneSuggestNext(demo);assert(strcmp(p.action,'blocked'));
end
fprintf('PASS: adapter baseline, queue provenance, pending observations, stale snapshots, held-out confirmation, fault/protocol separation and replay.\n');
end
