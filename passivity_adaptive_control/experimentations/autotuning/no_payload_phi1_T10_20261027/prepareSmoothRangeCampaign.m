function state=prepareSmoothRangeCampaign()
%PREPARESMOOTHRANGECAMPAIGN New immutable smoothness-first wider-range campaign.
% Offline only. Reuses verified original-filter recordings, never applies gains.
folder=fileparts(mfilename('fullpath'));pac=fileparts(fileparts(fileparts(folder)));
addpath(pac,fullfile(pac,'autotuning'),folder);
root=fullfile(folder,'campaign_smooth_v2');
assert(~isfolder(root),'pacTune:Exists','Campaign already exists.');
old=load(fullfile(folder,'campaign_length_v1','campaign.mat'),'C');
C=pacTuneSmoothCampaignDraft();C.id='unloaded_phi1_T10_smooth_v2';C.purpose='manual';
C.baselineTuning=old.C.baselineTuning;C.conditions=old.C.conditions;
C.conditions.codeRevision='Framework d6020cb plus archived smoothness-objective/fixed-analysis sources; controller kept at original logged filter.';
C.evaluationBudget=10;C.confirmationRepeats=2;C.minimumConfirmationDuration_s=120;
C.optimizer.initialDesignSize=1; % Eight historical observations seed the surrogate.
C.optimizer.continueSearchAfterFeasible=true; % User asks to optimize, not stop at a feasible baseline.
C.protocol=old.C.protocol;
C.provenance=struct('sourceCampaign','campaign_length_v1', ...
 'reuse','Eight old-filter screenings re-evaluated under new objective/2mm constraints; not prospective validation.', ...
 'reset','Standing user attestation; software states verified; no new physical observation.', ...
 'exclusion','filter_smoke changed-filter run explicitly excluded.', ...
 'restore','Operator must restore original measured filter before new trials; future logs will be replay-checked.');
C.budgetNotes='Ten observations total: eight reused plus two new exploratory runs. Any required fresh confirmations beyond that need budget review; no automatic extension.';
names={'baseline_00g_phi1_T10_12cycles_r01','baseline_00g_phi1_T10_12cycles_r02'};
for k=3:8,names{end+1}=sprintf('candidate%02d_00g_phi1_T10_12cycles_r01',k);end
trials=cell(8,1);hashes=strings(8,1);
for k=1:8
    trials{k}=pacTuneReadTrial(fullfile(folder,[names{k} '.mat']),C.evaluator);
    assertOriginalFilter(trials{k});hashes(k)=string(trials{k}.sha256);
end
C.provenance.sourceHashes=hashes;
pacTuneManualCampaign('create',root,C);
snapshot=fullfile(root,'source_snapshot');mkdir(snapshot);
files=[dir(fullfile(pac,'autotuning','*.m'));dir(fullfile(pac,'pacMeasuredDerivativesOriginal.m')); ...
    dir(fullfile(folder,'prepareSmoothRangeCampaign.m'))];
for k=1:numel(files),copyfile(fullfile(files(k).folder,files(k).name),fullfile(snapshot,files(k).name));end
m=pacTuneTrialManifest();m.conditions=C.conditions;m.resetVerified=true;
m.operator='user (standing full-reset attestation)';
m.notes='Historical old-filter run, retrospectively rescored under new objective/2mm limits. Original filter verified by exact recorded-input replay. Physical reset/calibration per standing user attestation; not reserved confirmation.';
for k=1:8
    if k<=2,id='baseline';else,id=sprintf('candidate%02d',k);end
    if k~=2
        b=pacTuneBaselineFromTrial(trials{k},20);
        pacTuneManualCampaign('propose',root,id,b.adaptationGamma([1 4 6]));
    end
    observation=pacTuneManualCampaign('import',root,sprintf('historical_%02d',k),id,trials{k}.file,m);
    fprintf('HISTORY %d %s class%s objective%.8g\n',k,id,observation.classification,observation.objective);
end
p=pacTuneSuggestNext(root);assert(strcmp(p.action,'propose'));
candidate=pacTuneQueueSuggestion(root,'smooth_candidate09',p);
save(fullfile(root,'next_candidate.mat'),'candidate','p','-v7');
state=pacTuneReportCampaign(root,fullfile(root,'initial_status'));
fprintf('NEXT gammaK %.17g gammaD %.17g gammaH %.17g reason%s\n',candidate.gamma,p.reason);
fprintf('PRED normalized mean%.8g std%.8g heuristicFeasibility%.8g\n',p.predictedMean,p.predictedStd,p.feasibilityScore);
fprintf('Best screened %s; remaining%d. No live gains applied.\n',state.bestScreenedCandidate,state.remainingEvaluations);
end

function assertOriginalFilter(trial)
S=load(trial.file,'out');ds=S.out.logsout;t=trial.signals.parameterEstimate.t;
q=pacTuneAlignSignal(trial.signals.meas_length,t,.075);
valid=pacTuneAlignSignal(trial.signals.orientationValid,t,.075)==1;
ix=find(strcmp(ds.getElementNames,'pollStart_s'));assert(isscalar(ix));e=ds{ix}.Values;
raw=double(e.Data);if e.IsTimeFirst,x=reshape(raw,numel(e.Time),[]);else,x=reshape(raw,[],numel(e.Time))';end
stamp=pacTuneAlignSignal(struct('t',double(e.Time(:)),'x',x),t,.075);
logged=pacTuneAlignSignal(trial.signals.filtered_d_length,t,.075);v=zeros(size(logged));
clear pacMeasuredDerivativesOriginal
for k=1:numel(t),[~,d,~]=pacMeasuredDerivativesOriginal(q(k,:)',valid(k),stamp(k));v(k,:)=d';end
assert(max(abs(v-logged),[],'all')<1e-12,'pacTune:FilterMismatch','Cannot reuse a different-filter recording.');
end
