function summary=prepareResidualData(root,outFolder)
%PREPARERESIDUALDATA Offline two-cohort preparation. No training/hardware.
assert(~isfolder(outFolder),'Choose a fresh output directory.');
addpath(root,fullfile(root,'RL_data_collection'));
audit=fullfile(root,'residual_dynamics_learning','audits','interfaces_20261101');
manifest=readtable(fullfile(audit,'cohort_provenance.csv'),'TextType','string');
cfg=struct('historyFrames',10,'previewFrames',5,'maxPhysicalGap_s',.10, ...
    'maxSimulationGap_s',.075,'maxCommandAge_s',.075,'velocityTau_s',.20, ...
    'accelerationTau_s',.30,'warmup_s',.60,'schema','residual-data-1', ...
    'boundMode',0); % Historical two_day_v2 contract; do not silently relabel.
mkdir(outFolder);save(fullfile(outFolder,'configuration.mat'),'cfg','manifest');
summary=table;trainStats=zeros(0,7);inverseStats=zeros(0,2);
for r=1:height(manifest)
    source=char(manifest.source(r));
    assert(strcmp(filehash(source),manifest.sha256(r)),'Raw source changed since audit.');
    id=char(manifest.runId(r));folder=fullfile(outFolder,id);mkdir(folder);
    options=struct('writeNative',false,'outputDirectory',fullfile(folder,'export'), ...
        'maxTransitionHostInterval_s',cfg.maxPhysicalGap_s, ...
        'metadata',struct('cohort',char(manifest.cohort(r)),'runId',id));
    exportRlTrainingData(source,options);
    O=readtable(fullfile(folder,'export','observations.csv'),'VariableNamingRule','preserve');
    S=struct;S.t=O.time_sim_s;
    assert(all(diff(S.t)>0),'Non-increasing sensor simulation clock.');
    assert(ismember('sensor_poll_start_host_s',O.Properties.VariableNames),'No physical clock.');
    S.physicalTime=O.sensor_poll_start_host_s;S.dt=[NaN;diff(S.physicalTime)];
    S.q=O{:,{'length_inferred_m_2','length_inferred_m_3'}};
    S.p=O{:,{'p1_command_bar','p2_command_bar','p3_command_bar'}};
    S.trial=O.trialId;S.channel=zeros(height(O),1);
    if ismember('channelId',O.Properties.VariableNames),S.channel=O.channelId;end
    S.phase=O.phaseId;S.split=repmat("excluded",height(O),1);
    positive=isfinite(S.trial)&S.trial>0;
    if manifest.cohort(r)=="day2_20261027"
        S.split(positive)="external_test";
    else
        S.split(positive)="train";
        % Whole trial held out per channel/mode stratum; never random rows.
        for channel=unique(S.channel(positive))'
            ids=unique(S.trial(positive&S.channel==channel));
            assert(numel(ids)>=2,'Need at least two trials in each development stratum.');
            held=ids(ceil(numel(ids)/2));
            S.split(positive&S.channel==channel&S.trial==held)="validation";
        end
    end
    ages=O{:,{'p1_command_bar_last_log_age_sim_s','p2_command_bar_last_log_age_sim_s','p3_command_bar_last_log_age_sim_s'}};
    S.valid=O.dynamics_observation_valid==1&O.base_orientation_valid==1& ...
        all(isfinite([S.q S.p]),2)&all(ages>=-1e-9&ages<=cfg.maxCommandAge_s,2)& ...
        all(S.p>=-1e-9&S.p<=3+1e-9,2)&isfinite(S.physicalTime)&positive;
    N=height(O);S.segment=zeros(N,1);S.mature=false(N,1);
    S.v=NaN(N,2);S.a=NaN(N,2);seg=0;start=0;
    for k=1:N
        if ~S.valid(k),continue;end
        continuous=k>1&&S.valid(k-1)&&S.trial(k)==S.trial(k-1)&& ...
            S.split(k)==S.split(k-1)&&S.channel(k)==S.channel(k-1)&& ...
            S.dt(k)>0&&S.dt(k)<=cfg.maxPhysicalGap_s&& ...
            S.t(k)-S.t(k-1)<=cfg.maxSimulationGap_s;
        if ~continuous
            seg=seg+1;start=k;S.v(k,:)=0;S.a(k,:)=0;
        else
            dt=S.dt(k);raw=(S.q(k,:)-S.q(k-1,:))/dt;
            S.v(k,:)=S.v(k-1,:)+dt/(cfg.velocityTau_s+dt)*(raw-S.v(k-1,:));
            ra=(S.v(k,:)-S.v(k-1,:))/dt;
            S.a(k,:)=S.a(k-1,:)+dt/(cfg.accelerationTau_s+dt)*(ra-S.a(k-1,:));
        end
        S.segment(k)=seg;S.mature(k)=S.physicalTime(k)-S.physicalTime(start)>=cfg.warmup_s;
    end
    S.nominalAcceleration=NaN(N,2);S.inverseDifferentialNominal=NaN(N,2);
    S.inversePressureNominal=NaN(N,3);S.inverseNominalLimited=false(N,1);
    for k=find(S.mature)'
        [a,pn,dn,lim]=residualNominal(S.q(k,:)',S.v(k,:)',S.a(k,:)',S.p(k,:)',cfg.boundMode);
        S.nominalAcceleration(k,:)=a';S.inverseDifferentialNominal(k,:)=dn';
        S.inversePressureNominal(k,:)=pn';S.inverseNominalLimited(k)=lim;
    end
    S.p0=mean(S.p,2);S.d=[S.p(:,2)-S.p(:,1),S.p(:,3)-S.p(:,1)];
    S.inverseResidualCandidate=S.d-S.inverseDifferentialNominal;
    S.accelerationResidualCandidate=S.a-S.nominalAcceleration;
    % Window index definitions, not large duplicated tensors. Every row in a
    % window must be mature, same segment and split; no gap bridging.
    forward=[];inverse=[];
    for k=cfg.historyFrames+1:N-cfg.previewFrames
        a=k-cfg.historyFrames;b=k+cfg.previewFrames;
        if S.segment(k)==0||~all(S.mature(a:b))||any(S.segment(a:b)~=S.segment(k)),continue;end
        forward(end+1,1)=k; %#ok<AGROW>
        if all(isfinite(S.inverseResidualCandidate(k,:))),inverse(end+1,1)=k;end %#ok<AGROW>
    end
    S.forwardAnchors=forward;S.inverseAnchors=inverse;
    S.runId=id;S.cohort=char(manifest.cohort(r));S.sha256=char(manifest.sha256(r));
    train=S.mature&S.split=="train";
    trainStats=[trainStats;S.q(train,:) S.v(train,:) S.p(train,:)]; %#ok<AGROW>
    take=inverse(S.split(inverse)=="train");
    inverseStats=[inverseStats;S.inverseResidualCandidate(take,:)]; %#ok<AGROW>
    save(fullfile(folder,'sequences.mat'),'S','cfg');
    rows=table(S.t,S.physicalTime,S.trial,S.channel,S.phase,S.split,S.valid, ...
        S.segment,S.mature,'VariableNames',{'time_sim_s','time_physical_s','trial', ...
        'channel','phase','split','valid','segment','mature'});
    writetable(rows,fullfile(folder,'row_manifest.csv'));
    anchors=unique([forward;inverse]);
    index=table(anchors,S.split(anchors),ismember(anchors,forward),ismember(anchors,inverse), ...
        'VariableNames',{'anchorRow','split','forward','inverse'});
    writetable(index,fullfile(folder,'window_index.csv'));
    row=table(string(id),manifest.cohort(r),N,nnz(S.valid),seg,numel(forward),numel(inverse), ...
        nnz(S.split(forward)=="train"),nnz(S.split(forward)=="validation"), ...
        nnz(S.split(forward)=="external_test"),'VariableNames', ...
        {'runId','cohort','rows','validRows','segments','forwardWindows','inverseWindows', ...
        'trainForward','validationForward','externalForward'});
    summary=[summary;row]; %#ok<AGROW>
    fprintf('Prepared %s: %d forward / %d inverse windows\n',id,numel(forward),numel(inverse));
end
assert(~isempty(trainStats)&&~isempty(inverseStats));
normalizer=struct('columns',{{'q2','q3','dq2','dq3','p1','p2','p3'}}, ...
    'mean',mean(trainStats,1),'std',std(trainStats,0,1), ...
    'inverseResidualMean',mean(inverseStats,1),'inverseResidualStd',std(inverseStats,0,1), ...
    'fitScope','day1 mature TRAIN rows only; inverse target scales on TRAIN anchors only');
normalizer.std(normalizer.std<1e-12)=1;
normalizer.inverseResidualStd(normalizer.inverseResidualStd<1e-12)=1;
save(fullfile(outFolder,'normalizer.mat'),'normalizer');
writetable(summary,fullfile(outFolder,'summary.csv'));
save(fullfile(outFolder,'preparation.mat'),'summary','cfg','manifest');
end
function h=filehash(path)
f=fopen(path,'rb');assert(f>=0);cleanup=onCleanup(@()fclose(f)); %#ok<NASGU>
bytes=fread(f,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(bytes,'int8'));h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
