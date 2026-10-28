function suggestion=pacTuneSuggestNext(root)
%PACTUNESUGGESTNEXT Read ledger and propose/confirm/observe/stop; no mutation.
state=pacTuneManualCampaign('status',root); % Also verifies immutable config hash.
assert(~isfile(fullfile(root,'.writer.lock')),'pacTune:Locked','Wait for the campaign writer.');
s=load(fullfile(root,'campaign.mat'),'C');C=s.C;
files=dir(fullfile(root,'events','*.mat'));[~,order]=sort({files.name});files=files(order);
assert(numel(files)==state.ledgerEventCount,'pacTune:StaleSuggestion','Ledger changed during read; retry.');
suggestion=struct('action','blocked','reason',state.stopReason,'gamma',NaN(3,1), ...
    'candidateId','','ledgerEventCount',numel(files), ...
    'configHash',pacTuneFileHash(fullfile(root,'campaign.mat')));
if state.success||state.remainingEvaluations<=0
    suggestion.action='stop';return
end
if ~state.configurationReady,suggestion.configurationIssues=state.configurationIssues;return;end
if ~isempty(state.pendingCandidateIds)
    suggestion.action='observe';suggestion.candidateId=char(state.pendingCandidateIds(1));return
end
if ~isempty(state.bestScreenedCandidate)
    suggestion.action='confirm';suggestion.candidateId=state.bestScreenedCandidate;return
end
assert(isfield(C,'optimizer'),'pacTune:OptimizerConfig', ...
    'Campaign lacks an optimizer snapshot; create a new immutable campaign.');
cfg=C.optimizer;cfg.gammaBounds=C.gammaBounds;cfg.maxObservations=C.evaluationBudget;
g=C.baselineTuning.adaptationGamma;
assert(all(g(1:3)==g(1))&&g(4)==g(5)&&all(g(6:8)==g(6)), ...
    'pacTune:OptimizerBaseline','Baseline must use grouped k/d/h gamma values.');
cfg.baselineGamma=g([1 4 6]);
o=pacTuneOptimizerObservations();excluded=struct('confirmation',0,'protocolOrData',0,'unknownFailure',0);
candidates=containers.Map('KeyType','char','ValueType','any');
for k=2:numel(files)
    s=load(fullfile(files(k).folder,files(k).name),'event');e=s.event;
    if strcmp(e.type,'candidate'),candidates(e.payload.id)=e.payload;continue;end
    if ~strcmp(e.type,'observation'),continue;end
    p=e.payload;
    if strcmp(p.kind,'confirmation'),excluded.confirmation=excluded.confirmation+1;continue;end
    candidate=candidates(p.candidateId);gamma=candidate.gamma;
    r=load(fullfile(p.artifactFolder,'evaluation.mat'),'result');
    objective=NaN;usable=false;known=false;feasible=false;
    % Do not train from protocol mismatches, even when a fault was logged.
    comparable=isempty(p.manifestIssues)&&isfield(r.result,'protocolMismatch')&& ...
        isempty(r.result.protocolMismatch)&&isempty(r.result.tuningEvents);
    if comparable
        if ismember(p.classification,{'valid_feasible','valid_infeasible'})&&isfinite(p.objective)
            objective=p.objective;usable=true;known=true;feasible=strcmp(p.classification,'valid_feasible');
        elseif strcmp(p.classification,'aborted_or_faulted')
            known=isfield(r.result,'firstFault_s')&&isfinite(r.result.firstFault_s);
            if ~known,excluded.unknownFailure=excluded.unknownFailure+1;end
        else
            excluded.protocolOrData=excluded.protocolOrData+1;
        end
    else
        excluded.protocolOrData=excluded.protocolOrData+1;
    end
    duration=NaN;
    if isfield(r.result,'diagnostics')
        time=r.result.diagnostics.t;step=median(diff(time));
        duration=max(0,min(time(end)+step,C.evaluator.trackingWindow_s(2))- ...
            max(time(1),C.evaluator.trackingWindow_s(1)));
    end
    o(end+1,:)={string(p.candidateId),gamma(1),gamma(2),gamma(3),objective,usable,known,feasible,duration}; %#ok<AGROW>
end
next=pacTuneOptimizerSuggest(cfg,o);
fields=fieldnames(next);
for k=1:numel(fields),suggestion.(fields{k})=next.(fields{k});end
suggestion.excluded=excluded;
suggestion.trainingCandidateIds=o.candidateId;
% All files are append-only. A later writer is also detected atomically on queue.
assert(~isfile(fullfile(root,'.writer.lock')),'pacTune:Locked','Campaign changed during read; retry.');
end
