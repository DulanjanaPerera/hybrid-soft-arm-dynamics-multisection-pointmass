function output=pacTuneManualCampaign(action,root,varargin)
%PACTUNEMANUALCAMPAIGN Offline manual candidate/manifest/observation ledger.
% create(root,config), propose(root,id,[gK gD gH]),
% import(root,trialId,candidateId,matFile,manifest), status(root).
% No model launch, tuning write, optimizer or hardware interface.
root=char(java.io.File(root).getCanonicalPath());
if strcmp(action,'create')
    C=varargin{1};validateConfig(C);
    assert(~isfolder(root)&&~isfile(root),'pacTune:Exists','Campaign root already exists.');
    mkdir(root);mkdir(fullfile(root,'events'));mkdir(fullfile(root,'artifacts'));mkdir(fullfile(root,'candidates'));
    save(fullfile(root,'campaign.mat'),'C','-v7');
    event=struct('type','created','configHash',pacTuneFileHash(fullfile(root,'campaign.mat')), ...
        'utc',utcNow(),'schemaVersion','pac-ledger-1');
    append(root,event,1);output=status(root);return
end
assert(isfolder(root),'pacTune:MissingCampaign','Campaign does not exist.');
if strcmp(action,'status'),output=status(root);return;end
% Exclusive writer lock; a failed process leaves a visible lock for inspection.
lock=java.io.File(fullfile(root,'.writer.lock'));
assert(lock.createNewFile(),'pacTune:Locked','Campaign writer lock exists; inspect before recovery.');
cleanup=onCleanup(@()deleteLock(lock));
[C,events]=read(root);state=status(root);
assert(isfinite(C.evaluationBudget),'pacTune:BudgetUnconfigured','Set a finite evaluation budget in a new campaign.');
assert(numel(state.trialIds)<C.evaluationBudget,'pacTune:Budget','Evaluation budget exhausted.');
switch action
    case 'propose'
        id=varargin{1};identifier(id);g=double(varargin{2}(:));
        assert(numel(g)==3&&all(isfinite(g))&&all(g>0),'pacTune:Gamma','Need three positive finite gamma values.');
        assert(all(isfinite(C.gammaBounds(:)))&&all(C.gammaBounds(:,1)>0)&& ...
            all(C.gammaBounds(:,2)>=C.gammaBounds(:,1)),'pacTune:Bounds','Set valid positive gamma bounds first.');
        assert(all(g>=C.gammaBounds(:,1)&g<=C.gammaBounds(:,2)),'pacTune:Bounds','Candidate is outside configured bounds.');
        assert(~ismember(string(id),state.candidateIds),'pacTune:Duplicate','Candidate ID already exists.');
        optimizer=[];
        if numel(varargin)>=3
            optimizer=varargin{3};
            assert(strcmp(optimizer.action,'propose')&&isequal(g,optimizer.gamma(:)), ...
                'pacTune:OptimizerProposal','Proposal and gamma differ.');
            assert(optimizer.ledgerEventCount==numel(events)&& ...
                strcmp(optimizer.configHash,pacTuneFileHash(fullfile(root,'campaign.mat'))), ...
                'pacTune:StaleSuggestion','Ledger changed since suggestion; recompute before queueing.');
        end
        tuning=C.baselineTuning;validateBaseline(tuning);
        tuning.adaptationGamma(1:3)=g(1);tuning.adaptationGamma(4:5)=g(2);tuning.adaptationGamma(6:8)=g(3);
        validateBaseline(tuning);
        payload=struct('id',char(id),'gamma',g,'tuning',tuning,'ready',state.configurationReady, ...
            'configurationIssues',state.configurationIssues,'optimizer',optimizer,'utc',utcNow());
        writeCandidate(fullfile(root,'candidates',char(id)),payload);
        event=struct('type','candidate','payload',payload,'utc',utcNow());
        append(root,event,numel(events)+1);
        output=payload;
    case 'import'
        trialId=varargin{1};candidateId=varargin{2};file=varargin{3};manifest=varargin{4};
        identifier(trialId);identifier(candidateId);
        assert(~ismember(string(trialId),state.trialIds),'pacTune:Duplicate','Trial ID already exists.');
        candidates=events(cellfun(@(e)strcmp(e.type,'candidate'),events));
        chosen=find(cellfun(@(e)strcmp(e.payload.id,candidateId),candidates));
        assert(isscalar(chosen),'pacTune:Candidate','Unknown candidate ID.');
        candidate=candidates{chosen}.payload;
        hash=pacTuneFileHash(file);
        assert(~ismember(string(hash),state.sourceHashes),'pacTune:Duplicate','This MAT file was already ingested.');
        issues=manifestIssues(manifest,C);
        cfg=C.evaluator;
        if strcmp(manifest.kind,'confirmation')&&~isempty(C.confirmationEvaluator),cfg=C.confirmationEvaluator;end
        cfg.expectedTuning=candidate.tuning;
        trial=pacTuneReadTrial(file,cfg);result=pacTuneEvaluateTrial(trial,cfg);
        assert(strcmp(hash,trial.sha256),'pacTune:ChangedFile','Source file changed during import.');
        % Logged startup can support a software reset check, not prove physical reset.
        if isfield(trial.signals,'parameterEstimate')
            initial=trial.signals.parameterEstimate.x(1,:);
            if ~isequal(initial,candidate.tuning.parameterInitial(:)')
                issues(end+1,1)="Logged initial parameter estimates differ from candidate.";
            end
        end
        if isfield(trial.signals,'historyState')&&any(trial.signals.historyState.x(1,:)~=0)
            issues(end+1,1)="Logged initial history is not zero.";
        end
        if isfield(result,'direction')&&isfinite(C.conditions.direction)&&result.direction~=C.conditions.direction
            issues(end+1,1)="Logged reference direction differs from campaign.";
        end
        if isfield(trial.signals,'des_phi')&&isfinite(C.conditions.phi_rad)
            sig=trial.signals.des_phi;mask=sig.t>=cfg.trackingWindow_s(1)&sig.t<cfg.trackingWindow_s(2);
            if ~any(mask)||any(abs(sig.x(mask)-C.conditions.phi_rad)>1e-6)
                issues(end+1,1)="Logged reference phi differs from campaign.";
            end
        end
        evaluatorClass=result.classification;eligible=result.success&&state.configurationReady&&isempty(issues);
        classification=evaluatorClass;
        if ~isempty(issues)&&~ismember(evaluatorClass,{'invalid_data','aborted_or_faulted'})
            classification='protocol_deviation';
        elseif ~state.configurationReady&&strcmp(evaluatorClass,'valid_feasible')
            classification='campaign_unconfigured';
        end
        % Store evaluator and campaign acceptance separately; never relabel raw metrics.
        result.campaignAcceptance=struct('eligible',eligible,'classification',classification, ...
            'manifestIssues',issues,'configurationIssues',state.configurationIssues);
        destination=fullfile(root,'artifacts',char(trialId));
        assert(~isfolder(destination),'pacTune:OrphanArtifact', ...
            'Trial artifact already exists without ledger entry; inspect before recovery.');
        pacTuneReportTrial(trial,result,destination);
        save(fullfile(destination,'manifest.mat'),'manifest','candidate','-v7');
        payload=struct('trialId',char(trialId),'candidateId',char(candidateId),'sourceHash',hash, ...
            'sourceFile',trial.file,'artifactFolder',destination,'kind',manifest.kind, ...
            'evaluatorClassification',evaluatorClass,'classification',classification, ...
            'eligible',eligible,'objective',NaN,'duration_s',diff(cfg.trackingWindow_s), ...
            'manifestIssues',issues,'utc',utcNow());
        if isfield(result,'objective'),payload.objective=result.objective;end
        event=struct('type','observation','payload',payload,'utc',utcNow());
        append(root,event,numel(events)+1);output=payload;
    otherwise
        error('pacTune:Action','Unknown action: %s.',action);
end
end

function validateConfig(C)
assert(strcmp(C.schemaVersion,'pac-manual-1'),'pacTune:Schema','Unknown campaign schema.');
identifier(C.id);
assert(ismember(C.purpose,{'draft','manual','offline_demo'}),'pacTune:Purpose','Invalid purpose.');
assert(isequal(size(C.gammaBounds),[3 2]),'pacTune:Bounds','Gamma bounds must be 3x2.');
assert(isnan(C.evaluationBudget)||(isfinite(C.evaluationBudget)&&C.evaluationBudget>=1&&fix(C.evaluationBudget)==C.evaluationBudget));
assert(isnan(C.confirmationRepeats)||(isfinite(C.confirmationRepeats)&&C.confirmationRepeats>=2&&fix(C.confirmationRepeats)==C.confirmationRepeats));
assert(isnan(C.minimumConfirmationDuration_s)||(isfinite(C.minimumConfirmationDuration_s)&&C.minimumConfirmationDuration_s>0));
pacTuneTargetPolicy(C.evaluator);
pacTuneObjectivePolicy(C.evaluator);
if ~isempty(C.confirmationEvaluator)
    a=C.evaluator;b=C.confirmationEvaluator;
    a=rmfield(a,{'trackingWindow_s','scoringWindow_s'});b=rmfield(b,{'trackingWindow_s','scoringWindow_s'});
    assert(isequaln(a,b),'pacTune:ConfirmationConfig','Only windows may differ for confirmation.');
end
end

function validateBaseline(b)
required={'adaptationGamma','parameterScale','parameterMaxRate','parameterLower','parameterUpper', ...
 'parameterInitial','adaptationOptions','hystoryOptions','pressureOptions','supervisorOptions', ...
 'lambda_gains','damping_gains','enableAdaptation','enableCommonPressureAdaptation','freezeHysteresisWeights'};
width=[10 10 10 10 10 10 3 4 10 7 2 2 1 1 1];
for k=1:numel(required)
    assert(isfield(b,required{k})&&isequal(size(b.(required{k})),[width(k) 1]), ...
        'pacTune:Baseline','Missing/wrong-size baseline %s.',required{k});
end
pacTuningConfig(b.adaptationGamma,b.parameterScale,b.parameterMaxRate,b.parameterLower, ...
    b.parameterUpper,b.parameterInitial,b.adaptationOptions,b.hystoryOptions,b.pressureOptions,b.supervisorOptions);
assert(all(isfinite([b.lambda_gains;b.damping_gains]))&&all([b.lambda_gains;b.damping_gains]>0));
assert(all(ismember([b.enableAdaptation;b.enableCommonPressureAdaptation;b.freezeHysteresisWeights],[0 1])));
end

function issues=configurationIssues(C)
issues=strings(0,1);
if ~strcmp(C.purpose,'manual'),issues(end+1,1)="Campaign purpose is not manual.";end
try
    validateBaseline(C.baselineTuning);
catch err
    issues(end+1,1)="Baseline: "+string(err.message);
end
if any(~isfinite(C.gammaBounds(:)))||any(C.gammaBounds(:,1)<=0)||any(C.gammaBounds(:,2)<C.gammaBounds(:,1))
    issues(end+1,1)="Positive finite gamma bounds are unset/invalid.";
end
if ~isfinite(C.evaluationBudget),issues(end+1,1)="Evaluation budget is unset.";end
if ~isfinite(C.confirmationRepeats),issues(end+1,1)="Confirmation repeat count is unset.";end
if ~isfinite(C.minimumConfirmationDuration_s),issues(end+1,1)="Confirmation duration is unset.";end
confirmation=C.evaluator;
if ~isempty(C.confirmationEvaluator),confirmation=C.confirmationEvaluator;end
if isfinite(C.minimumConfirmationDuration_s)&&diff(confirmation.trackingWindow_s)<C.minimumConfirmationDuration_s
    issues(end+1,1)="Confirmation tracking window is shorter than the minimum duration.";
end
if ~pacTuneTargetPolicy(C.evaluator),issues(end+1,1)="Required performance targets are unset.";end
if ~isfinite(C.evaluator.expectedPeriod_s)||C.evaluator.expectedPeriod_s<=0,issues(end+1,1)="Expected period is unset.";end
c=C.conditions;
if ~isfinite(c.payloadMass_kg)||c.payloadMass_kg<0,issues(end+1,1)="Payload mass is unset.";end
if ~ismember(c.payloadAttachedBeforeTrim,{'yes','no'}),issues(end+1,1)="Payload attachment timing is unknown.";end
for f={'calibrationId','modelId','resetProcedure','preconditioning','codeRevision'}
    if strlength(string(c.(f{1})))==0,issues(end+1,1)="Missing condition: "+string(f{1});end %#ok<AGROW>
end
if ~isfinite(c.phi_rad)||c.phi_rad<=0||~ismember(c.direction,[-1 1]),issues(end+1,1)="Trajectory conditions are unset.";end
end

function issues=manifestIssues(m,C)
template=pacTuneTrialManifest();issues=strings(0,1);
for f=fieldnames(template)'
    assert(isfield(m,f{1}),'pacTune:Manifest','Manifest missing %s.',f{1});
end
assert(strcmp(m.schemaVersion,template.schemaVersion),'pacTune:Manifest','Unknown manifest schema.');
assert(ismember(m.kind,{'screening','confirmation'}),'pacTune:Manifest','Invalid trial kind.');
if ~isequal(m.resetVerified,true),issues(end+1,1)="Operator has not verified reset procedure.";end
if isfield(m,'conditions')&&(~isstruct(m.conditions)||~isscalar(m.conditions))
    error('pacTune:Manifest','Manifest conditions must be a scalar struct.');
end
if strlength(string(m.operator))==0,issues(end+1,1)="Operator identifier is missing.";end
for f=fieldnames(C.conditions)'
    if ~isfield(m.conditions,f{1})||~isequaln(m.conditions.(f{1}),C.conditions.(f{1}))
        issues(end+1,1)="Condition mismatch: "+string(f{1}); %#ok<AGROW>
    end
end
end

function output=status(root)
[C,events]=read(root);issues=configurationIssues(C);
output=struct('id',C.id,'purpose',C.purpose,'configurationReady',isempty(issues), ...
    'configurationIssues',issues,'candidateIds',strings(0,1),'trialIds',strings(0,1), ...
    'sourceHashes',strings(0,1),'bestScreenedCandidate','','bestConfirmedCandidate','', ...
    'success',false,'stopReason','awaiting_candidate','nextAction','Propose a bounded candidate.');
observations=struct([]);candidates=struct([]);
for k=2:numel(events)
    e=events{k};
    if strcmp(e.type,'candidate')
        output.candidateIds(end+1,1)=string(e.payload.id);
        candidates=[candidates;struct('id',e.payload.id,'gammaK',e.payload.gamma(1), ...
            'gammaD',e.payload.gamma(2),'gammaH',e.payload.gamma(3),'configurationReady',e.payload.ready)]; %#ok<AGROW>
    end
    if strcmp(e.type,'observation')
        output.trialIds(end+1,1)=string(e.payload.trialId);
        output.sourceHashes(end+1,1)=string(e.payload.sourceHash);
        observations=[observations;e.payload]; %#ok<AGROW>
    end
end
output.observations=struct2table(observations,'AsArray',true);
output.candidates=struct2table(candidates,'AsArray',true);
output.ledgerEventCount=numel(events);
best=Inf;confirmedBest=Inf;
for id=output.candidateIds'
    if isempty(observations),continue;end
    matches=strcmp({observations.candidateId},id);eligible=[observations.eligible]&matches;
    allScreening=matches&strcmp({observations.kind},'screening');
    allConfirmation=matches&strcmp({observations.kind},'confirmation');
    screening=eligible&allScreening;
    confirmation=eligible&allConfirmation& ...
        [observations.duration_s]>=C.minimumConfirmationDuration_s;
    cleanScreening=any(screening)&&all(eligible(allScreening));
    failedConfirmation=any(allConfirmation&~confirmation);
    if cleanScreening&&~failedConfirmation
        score=mean([observations(screening).objective]);
        if score<best,best=score;output.bestScreenedCandidate=char(id);end
    end
    if cleanScreening&&nnz(confirmation)>=C.confirmationRepeats&&all(confirmation(allConfirmation))
        score=mean([observations(confirmation).objective]);
        if score<confirmedBest,confirmedBest=score;output.bestConfirmedCandidate=char(id);end
    end
end
output.remainingEvaluations=C.evaluationBudget-numel(output.trialIds);
observedIds=strings(0,1);
if ~isempty(observations),observedIds=string({observations.candidateId});end
output.pendingCandidateIds=setdiff(output.candidateIds,observedIds,'stable');
if ~isempty(output.bestConfirmedCandidate)
    output.success=true;output.stopReason='confirmed_target_met';output.nextAction='Review repeated confirmation and retain the candidate.';
elseif output.remainingEvaluations<=0
    output.stopReason='evaluation_budget_exhausted';output.nextAction='Review results; a revised budget requires a new campaign.';
elseif ~output.configurationReady
    output.stopReason='configuration_incomplete';output.nextAction='Resolve configuration issues in a new immutable campaign; offline diagnostics may continue.';
elseif ~isempty(output.bestScreenedCandidate)
    if isfield(C,'optimizer')&&isfield(C.optimizer,'continueSearchAfterFeasible')&&C.optimizer.continueSearchAfterFeasible
        if ~isempty(output.pendingCandidateIds)
            output.stopReason='awaiting_observation';output.nextAction='Import the pending objective-search trial; screened candidates remain unconfirmed.';
        else
            output.stopReason='awaiting_candidate';output.nextAction='Continue bounded objective search; screened candidates remain unconfirmed.';
        end
    else
        output.stopReason='awaiting_confirmation';output.nextAction='Acquire repeated operator-run confirmation for the screened candidate.';
    end
elseif ~isempty(output.pendingCandidateIds)
    output.stopReason='awaiting_observation';output.nextAction='Import an operator-run trial for a proposed candidate.';
elseif ~isempty(output.trialIds)
    output.stopReason='awaiting_candidate';output.nextAction='Review failed attempts and propose a revised bounded candidate or campaign.';
end
end

function [C,events]=read(root)
s=load(fullfile(root,'campaign.mat'),'C');C=s.C;validateConfig(C);
files=dir(fullfile(root,'events','*.mat'));[~,order]=sort({files.name});files=files(order);
assert(~isempty(files),'pacTune:Ledger','Missing creation event.');events=cell(numel(files),1);
for k=1:numel(files)
    assert(strcmp(files(k).name,sprintf('%06d.mat',k)),'pacTune:Ledger','Ledger sequence gap.');
    s=load(fullfile(files(k).folder,files(k).name),'event');events{k}=s.event;
end
assert(strcmp(events{1}.type,'created')&&strcmp(events{1}.configHash,pacTuneFileHash(fullfile(root,'campaign.mat'))), ...
    'pacTune:ChangedConfig','Campaign configuration changed after creation.');
end

function append(root,event,index)
target=fullfile(root,'events',sprintf('%06d.mat',index));
assert(~isfile(target),'pacTune:Ledger','Refusing to replace an event.');
temporary=[tempname(fullfile(root,'events')) '.pending'];
save(temporary,'event','-v7');[ok,msg]=movefile(temporary,target);assert(ok,'pacTune:Ledger','%s',msg);
end

function identifier(id)
assert((ischar(id)||(isstring(id)&&isscalar(id)))&& ...
    ~isempty(regexp(char(id),'^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$','once')), ...
    'pacTune:Id','IDs must contain 1-64 letters, numbers, underscores or hyphens.');
end

function t=utcNow()
t=char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function deleteLock(lock)
lock.delete();
end

function writeCandidate(folder,candidate)
assert(~isfolder(folder),'pacTune:OrphanArtifact','Candidate artifact exists; inspect before recovery.');
mkdir(folder);save(fullfile(folder,'candidate.mat'),'candidate','-v7');
fid=fopen(fullfile(folder,'review.md'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Manual PAC candidate %s\n\nConfiguration complete: %d. No gains have been applied.\n\n',candidate.id,candidate.ready);
for s=candidate.configurationIssues',fprintf(fid,'- %s\n',s);end
fprintf(fid,'\nRecorded proposal (review before any manually run trial):\n\n```matlab\n');
for f=fieldnames(candidate.tuning)'
    fprintf(fid,'%s = %s;\n',f{1},mat2str(candidate.tuning.(f{1}),17));
end
fprintf(fid,'```\n\nThe helper uses hystoryOptions to match existing log spelling. Map it to the live historyOptions input. lambda_gains/damping_gains are lambda/Ks inputs.\n');
if isfield(candidate,'optimizer')&&~isempty(candidate.optimizer)
    p=candidate.optimizer;
    fprintf(fid,'\n## Optimizer proposal provenance\n\nReason: %s. Ledger events at suggestion: %d.\n\n',p.reason,p.ledgerEventCount);
    fprintf(fid,'Predicted objective mean: %.6g; latent standard deviation: %.6g. Feasibility score: %.6g.\n\n', ...
        p.predictedMean,p.predictedStd,p.feasibilityScore);
    fprintf(fid,'These model-dependent diagnostics do not guarantee physical feasibility or stability. Full optimizer configuration and training counts are in candidate.mat.\n');
end
end
