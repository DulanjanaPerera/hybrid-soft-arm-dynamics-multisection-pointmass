function result=pacTuneEvaluateTrial(trial,cfg)
%PACTUNEEVALUATETRIAL Offline tracking, smoothness and constraint evaluation.
if nargin<2,cfg=pacTuneCampaignConfig();end
[configured,requiredTargets,activeTargets]=pacTuneTargetPolicy(cfg);
lengthPolicy=strcmp(cfg.schemaVersion,'pac-tune-2');
assert(all(diff(cfg.trackingWindow_s)>0)&&all(diff(cfg.scoringWindow_s)>0));
assert(cfg.lengthNormalizer_mm>0&&cfg.minimumCoverage>0&&cfg.minimumCoverage<=1);
assert(cfg.scoringWindow_s(1)>=cfg.trackingWindow_s(1)&& ...
    cfg.scoringWindow_s(2)<=cfg.trackingWindow_s(2),'pacTune:Window','Scoring must be inside tracking.');
result=struct('schemaVersion',cfg.schemaVersion,'file',trial.file,'sha256',trial.sha256, ...
    'classification','invalid_data','success',false,'issues',strings(0,1));
result.config=cfg;
result.requiredTargets=requiredTargets;result.activeTargets=activeTargets;
if ~isempty(trial.missing)
    result.issues="Missing required signals: "+strjoin(trial.missing,', ');return
end
t=trial.signals.parameterEstimate.t;D=struct;result.extendedConstantInputs=strings(0,1);
for f=fieldnames(trial.signals)'
    signal=trial.signals.(f{1});
    if ismember(f{1},cfg.constantInputFields)&&all(all(signal.x==signal.x(1,:)))
        x=repmat(signal.x(1,:),numel(t),1);x(t<signal.t(1),:)=NaN;
        D.(f{1})=x;
        if signal.t(end)<t(end),result.extendedConstantInputs(end+1,1)=string(f{1});end
    else
        D.(f{1})=pacTuneAlignSignal(signal,t,cfg.maxAlignmentAge_s);
    end
end
N=numel(t);dt=median(diff(t));circle=t>=cfg.trackingWindow_s(1)&t<cfg.trackingWindow_s(2);
assert(any(circle),'pacTune:Window','No data in tracking window.');
% Reject gaps/nonfinite aligned required data; Inf is valid only as stop time.
alignedGood=true(N,1);
for f=fieldnames(D)'
    if ismember(f{1},trial.requiredFields)
        if lengthPolicy&&ismember(f{1},{'theta','phi','thetaObservable'}),continue;end
        x=D.(f{1});
        if strcmp(f{1},'supervisorOptions'),x=x(:,1:6);end
        alignedGood=alignedGood&all(isfinite(x),2);
    end
end
% All finite guards precede Boolean comparisons: logical(NaN) would be true.
age=t-D.measurementTime_s;
fresh=isfinite(age)&age>=-1e-9&age<=D.supervisorOptions(:,5)+1e-9;
valid=alignedGood&D.modelValid==1&D.orientationValid==1&D.configurationValid==1& ...
    D.baseOreintationValid==1&fresh&D.status==3;
q=D.meas_length(:,2:3);qd=D.sampledDes_length;
qerr=qd-q;verr=D.d_length-D.filtered_d_length;
thetaerr=atan2(sin(D.des_theta-D.theta),cos(D.des_theta-D.theta));
phierr=D.des_phi-D.phi;
qvalid=valid&all(isfinite(qerr),2);phivalid=valid&isfinite(phierr);
velocityValid=qvalid&all(isfinite(verr),2);
velocityPair=velocityValid&[false;velocityValid(1:end-1)]& ...
    [false;diff(t)<=cfg.maxAlignmentAge_s+1e-9];
velocityErrorRate=[NaN(1,2);diff(verr)./diff(t)];
thvalid=valid&D.thetaObservable==1&isfinite(thetaerr);
pair=[false;thvalid(1:end-1)]&thvalid&[false;diff(t)<=cfg.maxAlignmentAge_s+1e-9];
omega=[NaN;atan2(sin(diff(D.theta)),cos(diff(D.theta)))./diff(t)];omega(~pair)=NaN;
refRate=[NaN;diff(D.des_theta)./diff(t)];
motion=circle&isfinite(refRate)&abs(refRate)>1e-6;
result.period_s=2*pi/median(abs(refRate(motion)));
assert(isfinite(result.period_s)&&result.period_s>0,'pacTune:Reference','No circular reference motion.');
result.direction=sign(median(refRate(motion)));
result.tracking=windowMetrics(cfg.trackingWindow_s);
result.scoring=windowMetrics(cfg.scoringWindow_s);
result.objective=result.scoring.qVectorRmse_mm/cfg.lengthNormalizer_mm;
result.firstFault_s=NaN;i=find(circle&D.status==5,1);
if ~isempty(i),result.firstFault_s=t(i);end
% Logged changes include Inf equality; structural NaNs are never silently equal.
tuning={'adaptationGamma','parameterLower','parameterUpper','parameterInitial', ...
 'parameterMaxRate','parameterScale','adaptationOptions','hystoryOptions', ...
 'pressureOptions','supervisorOptions','lambda_gains','damping_gains', ...
 'enableAdaptation','enableCommonPressureAdaptation','freezeHysteresisWeights'};
eventSignal=strings(0,1);eventTime=zeros(0,1);eventValues=cell(0,1);
for f=tuning
    signal=trial.signals.(f{1});x=signal.x;
    changed=[false;any(x(2:end,:)~=x(1:end-1,:),2)];
    inWindow=signal.t>=cfg.trackingWindow_s(1)&signal.t<cfg.trackingWindow_s(2);
    for j=find(inWindow&changed)'
        eventSignal(end+1,1)=string(f{1});eventTime(end+1,1)=signal.t(j);eventValues{end+1,1}=x(j,:); %#ok<AGROW>
    end
end
result.tuningEvents=table(eventSignal,eventTime,eventValues);
result.initialTuning=struct;
for f=tuning,result.initialTuning.(f{1})=D.(f{1})(find(circle,1),:);end
result.protocolMismatch=strings(0,1);
for f=fieldnames(cfg.expectedTuning)'
    if ~isfield(D,f{1})||any(any(D.(f{1})(circle,:)~=reshape(cfg.expectedTuning.(f{1}),1,[])))
        result.protocolMismatch(end+1,1)=string(f{1});
    end
end
if isfinite(cfg.expectedPeriod_s)&&abs(result.period_s-cfg.expectedPeriod_s)>1e-3
    result.protocolMismatch(end+1,1)="period";
end
% Bounds refer to stored (pre-update) estimates. Measured rates use next output.
th=D.parameterEstimate;
result.lowerBoundFraction=mean(th(circle,:)<=D.parameterLower(circle,:)+1e-8,1);
result.upperBoundFraction=mean(th(circle,:)>=D.parameterUpper(circle,:)-1e-8,1);
eigmin=(th(:,1)+th(:,2)+2*th(:,3)-sqrt((th(:,1)-th(:,2)).^2+4*th(:,3).^2))/2;
result.minimumEstimatedKEigenvalue_NPerM=min(eigmin(circle));
result.nonpositiveKFraction=mean(eigmin(circle)<=0);
result.pressureRange_bar=[min(D.pressureCommand_bar(circle,:),[],1);max(D.pressureCommand_bar(circle,:),[],1)];
pressureRate=diff(D.pressureCommand_bar)./diff(t);
pressurePair=circle(1:end-1)&circle(2:end)&diff(t)<=cfg.maxAlignmentAge_s;
result.pressurePeakCommandRate_barPerS=max(abs(pressureRate(pressurePair,:)),[],1);
result.nominalResidualPeak_N=max(vecnorm(D.allocationResidual_N(circle,:),2,2));
result.alignedCoverage=mean(alignedGood(circle));
% Reconstruct actual update dt and first-failing gate from the wrapper logic.
good=D.configurationValid==1&D.baseOreintationValid==1&fresh;
updateDt=zeros(N,1);last=-1;
for j=1:N
    if good(j)&&D.measurementTime_s(j)>last
        if last>=0,updateDt(j)=D.measurementTime_s(j)-last;end
        last=D.measurementTime_s(j);
    end
end
reason=zeros(N,1);
if isfield(D,'s_mPerS')
    checks={D.enableAdaptation~=1,~good|D.modelValid~=1,D.status~=3,D.limited~=0, ...
        vecnorm(D.allocationResidual_N,2,2)>D.adaptationOptions(:,2), ...
        updateDt<=0|updateDt>D.adaptationOptions(:,3), ...
        vecnorm(D.s_mPerS,2,2)<=D.adaptationOptions(:,1)};
    labels=["eligible projected update";"disabled";"invalid freshness/model"; ...
        "supervisor";"limited";"nominal residual";"update interval";"deadband";"zero projected rate"];
    for k=1:numel(checks),reason(reason==0&checks{k})=k;end
    reason(reason==0&D.adaptationActive~=1)=8;
    counts=arrayfun(@(k)nnz(circle&reason==k),(0:8)');
    result.gates=table(labels,counts);
    result.gateActiveMismatch=nnz(circle&reason~=0&D.adaptationActive==1);
else
    result.gates=table;result.gateActiveMismatch=NaN;
end
rates=diff(th)./updateDt(1:end-1);
rateMask=circle(1:end-1)&D.adaptationActive(1:end-1)==1&updateDt(1:end-1)>0;
caps=D.parameterMaxRate(1:end-1,:);
result.rates=struct('activeTransitions',nnz(rateMask), ...
    'capFraction',mean(abs(rates(rateMask,:))>=caps(rateMask,:)-1e-6 & caps(rateMask,:)>0,1), ...
    'peakAbsolutePerSecond',max(abs(rates(rateMask,:)),[],1));
result.unclippedRateAvailable=all(isfield(D,{'historyState','s_mPerS'}));
if result.unclippedRateAvailable
    v=D.d_length-D.lambda_gains.*(q-qd);z=D.historyState;s=D.s_mPerS;
    gradient=[q(:,1).*s(:,1),q(:,2).*s(:,2),sum(q,2).*sum(s,2), ...
        v.*s,-z(:,1).*sum(s,2),z(:,2).*s(:,1),z(:,3).*s(:,2),s];
    gamma=D.adaptationGamma;gamma(D.freezeHysteresisWeights==1,6:8)=0;
    rawRate=-D.parameterScale.^2.*gamma.*gradient;
    eligible=circle&(reason==0|reason==8);
    result.rates.unclippedPeakAbsolutePerSecond=max(abs(rawRate(eligible,:)),[],1);
    result.rates.unclippedAboveCapFraction=mean(abs(rawRate(eligible,:))>D.parameterMaxRate(eligible,:),1);
end
result.forceRms_N=struct;
for f={'dampingForce_N','hysteresisForce_N','feedbackForce_N','inertialForce_N','staticForce_N'}
    if isfield(D,f{1}),result.forceRms_N.(f{1})=sqrt(mean(D.(f{1})(circle&valid,:).^2,1));end
end
result.parameterChangeTracking=th(find(circle,1,'last'),:)-th(find(circle,1),:);
result.maximumObservedStep_s=max(diff(t));
% Exact nominal allocation replay distinguishes feasible target vs slew mismatch.
allocation=false(N,1);slew=false(N,1);replayed=false(N,1);
if all(isfield(D,{'tauBlended_N','forceBias_N'}))
    fields={'common_bar','sTrim_bar','aTrim_bar','minimum_bar','maximum_bar', ...
        'maxDeltaS_bar','maxDeltaA_bar','slew_barPerS','activation_bar','forcePerBar_N'};
    for j=find(circle&valid&D.enableCommonPressureAdaptation==1)'
        p=D.pressureOptions(j,:);pc=pacPressureConfig();
        for k=1:10,pc.(fields{k})=p(k);end
        [target,~,~,limited]=pacAllocateVariableCommon(D.tauBlended_N(j,:)',D.forceBias_N(j,:)',pc);
        allocation(j)=limited;slew(j)=any(abs(target'-D.pressureCommand_bar(j,:))>1e-9);replayed(j)=true;
    end
end
result.allocationReplay=struct('samples',nnz(replayed),'allocationLimited',nnz(allocation), ...
 'slewTargetMismatch',nnz(slew),'loggedFlagMismatch',nnz(replayed&(D.limited~=(allocation|slew))));
% Per-cycle windows come from actual unwrapped reference speed, not filename.
rows=struct([]);nLaps=floor(diff(cfg.trackingWindow_s)/result.period_s+1e-6);
for k=1:nLaps
    w=cfg.trackingWindow_s(1)+[k-1 k]*result.period_s;
    if w(1)>t(end),break;end
    m=windowMetrics(w);m.lap=k;rows=[rows;m]; %#ok<AGROW>
end
result.laps=struct2table(rows);
% Hard outcomes precede protocol/performance labels. No configured targets = no success.
targetNames=fieldnames(cfg.targets);
result.targetsConfigured=configured;result.violations=strings(0,1);
for f=targetNames'
    value=result.scoring.(f{1});limit=cfg.targets.(f{1});
    if isfinite(limit)&&(~all(isfinite(value))||any(value>limit))
        result.violations(end+1,1)=string(f{1});
    end
end
complete=t(1)<=cfg.trackingWindow_s(1)+dt/2&&t(end)+dt>=cfg.trackingWindow_s(2)-1e-8;
coverageBad=result.alignedCoverage<cfg.minimumCoverage|| ...
    result.tracking.coverage<cfg.minimumCoverage||result.scoring.coverage<cfg.minimumCoverage;
if lengthPolicy
    coverageBad=coverageBad||result.scoring.velocityCoverage<cfg.minimumCoverage|| ...
        result.scoring.velocityPairCoverage<cfg.minimumCoverage;
    if any(ismember(activeTargets,{'thetaMae_rad','thetaRateStd_radPerS'}))
        coverageBad=coverageBad||result.scoring.thetaCoverage<cfg.minimumCoverage;
    end
    if ismember('phiMae_rad',activeTargets),coverageBad=coverageBad||result.scoring.phiCoverage<cfg.minimumCoverage;end
else
    coverageBad=coverageBad||result.scoring.thetaCoverage<cfg.minimumCoverage;
end
if isfinite(result.firstFault_s)||~complete
    result.classification='aborted_or_faulted';
elseif coverageBad
    result.classification='invalid_data';
elseif ~isempty(eventTime)||~isempty(result.protocolMismatch)
    result.classification='protocol_deviation';
elseif ~isempty(result.violations)
    result.classification='valid_infeasible';
elseif ~configured
    result.classification='unconfigured_targets';
else
    result.classification='valid_feasible';result.success=true;
end
if ~configured,result.issues(end+1,1)="Targets are unset; historical diagnostics cannot certify success.";end
if ~complete,result.issues(end+1,1)="Recording ends before planned tracking completes.";end
if coverageBad
    result.issues(end+1,1)="Coverage below configured minimum of "+string(cfg.minimumCoverage)+".";
end
if ~isempty(eventTime),result.issues(end+1,1)="Live tuning changes detected; exploratory trial.";end
if ~isempty(result.protocolMismatch),result.issues(end+1,1)="Expected configuration mismatch: "+strjoin(result.protocolMismatch,', ');end
result.diagnostics=table(t,valid,thvalid,omega,refRate,eigmin,allocation,slew,replayed,reason,updateDt);

    function m=windowMetrics(w)
        ix=t>=w(1)-1e-8&t<w(2)-1e-8;expected=max(1,round(diff(w)/dt));
        iq=ix&qvalid;it=ix&thvalid;ip=ix&phivalid;iv=ix&velocityValid;ir=ix&velocityPair;
        m=struct('start_s',w(1),'end_s',w(2),'samples',nnz(ix), ...
            'coverage',nnz(iq)/expected,'thetaCoverage',nnz(it)/expected, ...
            'phiCoverage',nnz(ip)/expected,'velocityCoverage',nnz(iv)/expected, ...
            'velocityPairCoverage',nnz(ir)/expected, ...
            'q2Rmse_mm',1000*sqrt(mean(qerr(iq,1).^2)), ...
            'q3Rmse_mm',1000*sqrt(mean(qerr(iq,2).^2)), ...
            'qVectorRmse_mm',1000*sqrt(mean(sum(qerr(iq,:).^2,2))), ...
            'thetaMae_rad',mean(abs(thetaerr(it))),'phiMae_rad',mean(abs(phierr(ip))), ...
            'phiMax_rad',maxOrNaN(D.phi(ip)), ...
            'velocityRmse_mPerS',sqrt(mean(sum(verr(iv,:).^2,2))), ...
            'velocityErrorRateRmse_mPerS2',sqrt(mean(sum(velocityErrorRate(ir,:).^2,2))), ...
            'thetaRateStd_radPerS',std(omega(ix&pair),'omitnan'), ...
            'nearStopFraction',mean(abs(omega(ix&pair))<.2*abs(refRate(ix&pair))), ...
            'slowOrReverseFraction',mean(sign(refRate(ix&pair)).*omega(ix&pair)<.2*abs(refRate(ix&pair))), ...
            'activeFraction',mean(D.adaptationActive(ix)==1), ...
            'limitedFraction',mean(D.limited(ix)==1),'faultSamples',nnz(ix&D.status==5));
    end
end

function x=maxOrNaN(x)
if isempty(x),x=NaN;else,x=max(x);end
end
