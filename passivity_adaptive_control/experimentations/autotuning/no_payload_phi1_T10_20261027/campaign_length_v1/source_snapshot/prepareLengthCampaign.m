function state=prepareLengthCampaign()
%PREPARELENGTHCAMPAIGN Archive agreed manual campaign and existing baselines.
% No hardware/model execution or automatic live tuning. Refuses overwrite.
folder=fileparts(mfilename('fullpath'));pac=fileparts(fileparts(fileparts(folder)));
repo=fileparts(pac);addpath(pac,fullfile(pac,'autotuning'),folder);
root=fullfile(folder,'campaign_length_v1');
assert(~isfolder(root),'pacTune:Exists','Campaign exists; use status/suggest rather than recreating it.');
C=pacTuneManualConfig();C.id='unloaded_phi1_T10_length_v1';C.purpose='manual';
C.evaluator=agreedEvaluationConfig();C.gammaBounds=[.25 .75;.001 .008;2.5 7.5];
C.evaluationBudget=10;C.confirmationRepeats=2;C.minimumConfirmationDuration_s=120;
source=fullfile(folder,'baseline_00g_phi1_T10_12cycles_r01.mat');
trial=pacTuneReadTrial(source,C.evaluator);C.baselineTuning=pacTuneBaselineFromTrial(trial,20);
model=fullfile(repo,'PAC_autotuning_manual.slx');modelHash=pacTuneFileHash(model);
C.conditions=struct('payloadMass_kg',0,'payloadAttachedBeforeTrim','yes', ...
 'calibrationId','unloaded_full_reset_recalibration_trim_procedure_v1', ...
 'modelId',['PAC_autotuning_manual.slx SHA256 ' modelHash], ...
 'resetProcedure','User-attested full depressurization, fresh model restart and recalibration/trim between runs.', ...
 'preconditioning','Same 20 s startup after full reset; 5 s ramp-down and 5 s pressurized hold, then manual pressure kill.', ...
 'phi_rad',1,'direction',1, ...
 'codeRevision','Framework base d6020cb plus archived local length-policy sources; controller model identified separately by saved SLX SHA256.');
C.protocol=struct('startup_s',20,'cycles',12,'period_s',10,'trackingEnd_s',140, ...
 'rampDownEnd_s',145,'recordingEnd_s',150,'pressureKill','manual after final hold');
C.budgetNotes='Ten observations including the two existing baseline screenings and any confirmations; review before adding evaluations. No new hardware runs were launched.';
C.provenance=struct('baseline1',source,'baseline1Hash',trial.sha256, ...
 'metadataSource','User attestation of full reset; calibration ID denotes a procedure, not an independently measured calibration record.', ...
 'modelSnapshotNote','Current saved duplicate snapshot; exact historical model bytes are not embedded in the baseline MAT files.');
pacTuneManualCampaign('create',root,C);
snapshot=fullfile(root,'source_snapshot');mkdir(snapshot);
files=[dir(fullfile(pac,'autotuning','*.m'));dir(fullfile(folder,'agreedEvaluationConfig.m'));dir(fullfile(folder,'prepareLengthCampaign.m'))];
for k=1:numel(files),copyfile(fullfile(files(k).folder,files(k).name),fullfile(snapshot,files(k).name));end
copyfile(model,fullfile(snapshot,'PAC_autotuning_manual.slx'));
g=C.baselineTuning.adaptationGamma([1 4 6]);pacTuneManualCampaign('propose',root,'baseline',g);
m=pacTuneTrialManifest();m.conditions=C.conditions;m.resetVerified=true;
m.operator='user (full-reset procedure attested in chat)';
m.notes='Previously recorded baseline, retrospectively evaluated under newly agreed length targets. Software reset verified from initial logs; physical reset/calibration per user attestation. Saved model snapshot is current, not an embedded historical artifact.';
for k=1:2
    file=fullfile(folder,sprintf('baseline_00g_phi1_T10_12cycles_r%02d.mat',k));
    pacTuneManualCampaign('import',root,sprintf('baseline_r%02d',k),'baseline',file,m);
end
p=pacTuneSuggestNext(root);assert(strcmp(p.action,'propose'));
candidate=pacTuneQueueSuggestion(root,'candidate03',p);
save(fullfile(root,'next_candidate.mat'),'candidate','p','-v7');
state=pacTuneReportCampaign(root,fullfile(root,'initial_status'));
fprintf('NEXT CANDIDATE gammaK %.17g gammaD %.17g gammaH %.17g\n',candidate.gamma);
fprintf('Remaining observations %d; pending %s. No gains applied.\n',state.remainingEvaluations,strjoin(state.pendingCandidateIds,', '));
end
