function prepareMotionPredictionPlots(root,folder)
%PREPAREMOTIONPREDICTIONPLOTS Export/audit Nov3; frozen-model plotting only.
assert(~isfolder(folder));mkdir(folder);
addpath(root,fullfile(root,'RL_data_collection'),fullfile(root,'residual_dynamics_learning'), ...
    fullfile(root,'residual_dynamics_learning','force_model'),fullfile(root,'residual_dynamics_learning','force_model','streaming'));
P=load(fullfile(root,'residual_dynamics_learning','datasets','two_day_v2','preparation.mat'),'cfg');cfg=P.cfg;
ids={'circleSweep','SineWave','helix','StepWave','combinePressure'};sources=table;quality=table;gaps=table;trials=table;
for j=1:numel(ids)
    id=ids{j};raw=fullfile(root,'RL_data_collection','20261103',[id '.mat']);hash=sgHash(raw);
    output=fullfile(folder,'exports',id);
    opts=struct('writeNative',false,'outputDirectory',output,'metadata', ...
        struct('cohort','20261103','purpose','frozen model visualization and gap audit; no retraining','rawPreserved',true));
    exportRlTrainingData(raw,opts);
    O=sgRead(fullfile(output,'observations.csv'));S=sequence(O,cfg,id);
    save(fullfile(folder,[id '_prepared.mat']),'S','O','cfg','-v7.3');
    sources=[sources;table(string(raw),string(hash),'VariableNames',{'source','sha256'})]; %#ok<AGROW>
    assert(strcmp(hash,sgHash(raw)));
    pressureAge=O{:,{'p1_command_bar_last_log_age_sim_s','p2_command_bar_last_log_age_sim_s','p3_command_bar_last_log_age_sim_s'}};
    quality=[quality;table(string(id),height(O),mean(O.orientation_valid>.5),mean(O.fresh_sensor_frame>.5), ...
        nnz(S.valid),nnz(S.mature),max(S.segment),median(S.dt(S.dt>0)),nnz(S.dt>cfg.maxPhysicalGap_s), ...
        nnz(any(pressureAge>cfg.maxCommandAge_s,2)), ...
        'VariableNames',{'runId','samples','orientationValidFraction','freshFrameFraction','validTrainingConventionRows', ...
        'matureRows','segments','medianPhysicalDt_s','longPhysicalIntervals','staleCommandRows'})]; %#ok<AGROW>
    % Native NDI usability, independent of trial/command eligibility.
    good=O.orientation_valid>.5&O.base_orientation_valid>.5&O.fresh_sensor_frame>.5&all(isfinite(S.q),2);
    starts=find(~good&[true;good(1:end-1)]);ends=find(~good&[good(2:end);true]);
    for g=1:numel(starts)
        a=starts(g);b=ends(g);bracketed=a>1&&b<height(O)&&good(a-1)&&good(b+1);
        duration=NaN;sameTrial=false;
        if bracketed,duration=S.physicalTime(b+1)-S.physicalTime(a-1);sameTrial=S.trial(a-1)==S.trial(b+1)&&S.channel(a-1)==S.channel(b+1);end
        tipMissing=NaN;reason=NaN;
        if ismember('ndi_tip_sensor_status',O.Properties.VariableNames),tipMissing=nnz(O.ndi_tip_sensor_status(a:b)==2);end
        if ismember('ndi_tip_rejection_reason',O.Properties.VariableNames),reason=mode(O.ndi_tip_rejection_reason(a:b));end
        gaps=[gaps;table(string(id),a,b,S.t(a),S.t(b),b-a+1,bracketed,duration,sameTrial,tipMissing,reason, ...
            'VariableNames',{'runId','firstRow','lastRow','startSim_s','endSim_s','missingSamples','bracketed', ...
            'bracketingPhysicalSpan_s','sameChannelTrial','tipMissingStatusRows','dominantTipRejectionReason'})]; %#ok<AGROW>
    end
    for tr=unique(S.trial(S.trial>0))'
        take=S.trial==tr;phase=take&ismember(S.phase,[5 6 7]);cycles=unique(O.cycleId(phase));
        starts=[];for c=cycles',k=find(phase&O.cycleId==c,1);starts(end+1)=S.t(k);end %#ok<AGROW>
        period=NaN;if numel(starts)>1,period=median(diff(starts));end
        maxPhi=NaN;if ismember('phi_reference_rad',O.Properties.VariableNames),maxPhi=max(O.phi_reference_rad(take));end
        trials=[trials;table(string(id),tr,S.channel(find(take,1)),min(S.t(take)),max(S.t(take)), ...
            numel(cycles),period,maxPhi,mean(good(take)), ...
            'VariableNames',{'runId','trial','channel','startSim_s','endSim_s','loggedCycles','cycleSpacingSim_s','maximumRequestedPhi_rad','nativeUsableFraction'})]; %#ok<AGROW>
    end
    fprintf('PREPARED %s: %d rows, %d mature, %d native NDI gap episodes.\n',id,height(O),nnz(S.mature),numel(ends));
end
writetable(sources,fullfile(folder,'raw_sources.csv'));writetable(quality,fullfile(folder,'quality.csv'));
writetable(gaps,fullfile(folder,'ndi_gaps.csv'));writetable(trials,fullfile(folder,'recorded_trials.csv'));
save(fullfile(folder,'preparation.mat'),'cfg','ids');disp(cfg);disp(quality);
end
function S=sequence(O,cfg,id)
S=struct('t',O.time_sim_s,'physicalTime',O.sensor_poll_start_host_s,'runId',id,'cohort','day4_20261103');
S.dt=[NaN;diff(S.physicalTime)];S.q=O{:,{'length_inferred_m_2','length_inferred_m_3'}};
S.p=O{:,{'p1_command_bar','p2_command_bar','p3_command_bar'}};S.trial=O.trialId;S.phase=O.phaseId;
S.channel=zeros(height(O),1);if ismember('channelId',O.Properties.VariableNames),S.channel=O.channelId;end
S.split=repmat("visualization_only",height(O),1);
ages=O{:,{'p1_command_bar_last_log_age_sim_s','p2_command_bar_last_log_age_sim_s','p3_command_bar_last_log_age_sim_s'}};
S.valid=O.dynamics_observation_valid==1&O.base_orientation_valid==1& ...
    all(isfinite([S.q S.p]),2)&all(ages>=-1e-9&ages<=cfg.maxCommandAge_s,2)& ...
    all(S.p>=-1e-9&S.p<=3+1e-9,2)&isfinite(S.physicalTime)&isfinite(S.trial)&S.trial>0;
N=height(O);S.segment=zeros(N,1);S.mature=false(N,1);S.v=NaN(N,2);seg=0;first=0;
for k=1:N
    if ~S.valid(k),continue;end
    continuous=k>1&&S.valid(k-1)&&S.trial(k)==S.trial(k-1)&&S.channel(k)==S.channel(k-1)&& ...
        S.dt(k)>0&&S.dt(k)<=cfg.maxPhysicalGap_s&&S.t(k)-S.t(k-1)<=cfg.maxSimulationGap_s;
    if ~continuous,seg=seg+1;first=k;S.v(k,:)=0;
    else,dt=S.dt(k);raw=(S.q(k,:)-S.q(k-1,:))/dt;S.v(k,:)=S.v(k-1,:)+dt/(cfg.velocityTau_s+dt)*(raw-S.v(k-1,:));end
    S.segment(k)=seg;S.mature(k)=S.physicalTime(k)-S.physicalTime(first)>=cfg.warmup_s;
end
end
