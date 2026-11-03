function D=sgExternalData(root,cfg)
%SGEXTERNALDATA Evaluation only; no normalization or fitting.
base=fullfile(root,'residual_dynamics_learning','validations','fresh_20261102_v2');
ids={'StepWave','SineWave4','combinePressure','circleSweep','helix'};
D=struct('sequences',{{}},'segments',struct([]),'index',table,'sources',table);
for run=1:numel(ids)
    path=fullfile(base,[ids{run} '_sequences.mat']);a=load(path);S=a.S;
    assert(strcmp(S.cohort,'day3_20261102')&&all(S.split=="fresh_test"));D.sequences{run}=S;
    D.sources=[D.sources;table(string(path),string(sgHash(path)),'VariableNames',{'source','sha256'})];
    raw=fullfile(root,'RL_data_collection','20261102',[ids{run} '.mat']);
    D.sources=[D.sources;table(string(raw),string(sgHash(raw)),'VariableNames',{'source','sha256'})];
    anchors=[];sids=[];
    for sid=unique(S.segment(S.mature))'
        rows=find(S.mature&S.segment==sid);first=rows(1);last=rows(end);
        k=find(S.physicalTime>=S.physicalTime(first)+cfg.commonBurnIn_s,1);
        if isempty(k)||k>=last,continue;end
        take=k:last-1;take=take(S.physicalTime(last)-S.physicalTime(take)>=1);
        if isempty(take),continue;end
        ends=zeros(size(take));aux=ends;
        for j=1:numel(take)
            ends(j)=find(S.physicalTime<=S.physicalTime(take(j))+.5+1e-10,1,'last');
            aux(j)=find(S.physicalTime<=S.physicalTime(take(j))+1+1e-10,1,'last');
        end
        seg=struct('run',run,'sourceSegment',sid,'firstMature',first,'last',last, ...
            'starts',take,'ends',ends,'auxEnds',aux,'partition','fresh_test');
        D.segments=[D.segments;seg];anchors=[anchors take];sids=[sids repmat(numel(D.segments),size(take))]; %#ok<AGROW>
    end
    choose=unique(round(linspace(1,numel(anchors),min(12,numel(anchors)))));
    D.index=[D.index;table(repmat(run,numel(choose),1),anchors(choose)',sids(choose)', ...
        'VariableNames',{'run','anchor','segmentId'})];
end
end
