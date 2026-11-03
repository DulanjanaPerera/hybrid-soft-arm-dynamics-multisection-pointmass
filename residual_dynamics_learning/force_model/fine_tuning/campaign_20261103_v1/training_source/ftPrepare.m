function ftPrepare(root,folder,cfg)
%FTPREPARE Preserve old splits; deterministic whole-trial holdouts in Nov3.
assert(~isfolder(folder));mkdir(folder);base=fullfile(root,'residual_dynamics_learning');
old=load(fullfile(base,'force_model','streaming','data_20261103_v1','data.mat'));
sequences=old.sequences;segments=old.segments;index=old.index;normalizer=old.normalizer;
manifest=old.manifest;runGroup=ones(numel(sequences),1);addedManifest=table;trialSplit=table;
ids={'circleSweep','SineWave','helix','StepWave','combinePressure'};
prepared=fullfile(base,'figures','motion_predictions_20261103_v1');
raw=sgRead(fullfile(prepared,'raw_sources.csv'),'TextType','string');
assert(all(contains(raw.source,filesep+"20261103"+filesep)));
for j=1:height(raw),assert(strcmp(sgHash(raw.source(j)),raw.sha256(j)));end
for j=1:numel(ids)
    file=fullfile(prepared,[ids{j} '_prepared.mat']);a=load(file,'S');S=a.S;
    assert(strcmp(S.cohort,cfg.newCohort));S.runId=[cfg.newCohort '__' ids{j}];
    S.split=repmat("excluded",numel(S.t),1);positive=isfinite(S.trial)&S.trial>0;S.split(positive)="train";
    for ch=unique(S.channel(positive))'
        tr=unique(S.trial(positive&S.channel==ch));assert(numel(tr)>=2);
        held=tr(ceil(numel(tr)/2));S.split(positive&S.channel==ch&S.trial==held)="selection";
        for t=tr'
            take=S.channel==ch&S.trial==t;part=S.split(find(take,1));
            trialSplit=[trialSplit;table(string(S.runId),ch,t,part,nnz(take),nnz(S.mature&take), ...
                'VariableNames',{'runId','channel','trial','partition','rows','matureRows'})]; %#ok<AGROW>
        end
    end
    for ch=unique(S.channel(positive))'
        train=unique(S.trial(S.channel==ch&S.split=="train"));held=unique(S.trial(S.channel==ch&S.split=="selection"));
        assert(isempty(intersect(train,held)));
    end
    sequences{end+1}=S;run=numel(sequences);runGroup(run,1)=2;allAnchors=[];allIds=[];
    for sid=unique(S.segment(S.mature))'
        rows=find(S.mature&S.segment==sid);first=rows(1);last=rows(end);
        k=find(S.physicalTime>=S.physicalTime(first)+cfg.commonBurnIn_s,1);
        if isempty(k)||k>=last||S.segment(k)~=sid,continue;end
        starts=[];ends=[];auxEnds=[];
        while S.physicalTime(last)-S.physicalTime(k)>=cfg.maxHorizon_s
            e=find(S.physicalTime<=S.physicalTime(k)+cfg.chunk_s+1e-10,1,'last');
            aend=find(S.physicalTime<=S.physicalTime(k)+cfg.maxHorizon_s+1e-10,1,'last');
            if e<=k||aend>last,break;end
            starts(end+1)=k;ends(end+1)=e;auxEnds(end+1)=aend;k=e; %#ok<AGROW>
        end
        if isempty(starts),continue;end
        part=S.split(starts(1));assert(part=="train"||part=="selection");
        assert(all(S.split(rows)==part));
        seg=struct('run',run,'sourceSegment',sid,'firstMature',first,'last',last,'starts',starts, ...
            'ends',ends,'auxEnds',auxEnds,'partition',char(part));
        segments=[segments;seg];allAnchors=[allAnchors starts];allIds=[allIds repmat(numel(segments),size(starts))]; %#ok<AGROW>
    end
    eligible=S.split(allAnchors)=="selection";aa=allAnchors(eligible);ss=allIds(eligible);
    take=unique(round(linspace(1,numel(aa),min(cfg.selectionAnchorsPerRun,numel(aa)))));
    assert(~isempty(take),'No eligible Nov3 holdout for %s.',ids{j});
    index=[index;table(repmat(run,numel(take),1),aa(take)',ss(take)', ...
        'VariableNames',{'run','anchor','segmentId'})]; %#ok<AGROW>
    source=fullfile(root,'RL_data_collection','20261103',[ids{j} '.mat']);
    addedManifest=[addedManifest;table(string(cfg.newCohort),string(S.runId),string(source),string(sgHash(source)), ...
        string(file),string(sgHash(file)),'VariableNames',{'cohort','runId','source','sha256','preparedSource','preparedSHA256'})]; %#ok<AGROW>
end
assert(isequaln(sequences(1:10),old.sequences)&&isequaln(index(1:height(old.index),:),old.index));
assert(isequaln(normalizer,old.normalizer));assert(all(ismember(manifest.cohort,{'day1_root','day2_20261027'})));
coverage=table;
for run=1:numel(sequences)
    for part={'train','selection'}
        ss=find([segments.run]==run&strcmp({segments.partition},part{1}));chunks=0;
        for sid=ss,chunks=chunks+numel(segments(sid).starts);end
        coverage=[coverage;table(string(sequences{run}.runId),runGroup(run),string(part{1}),numel(ss),chunks, ...
            'VariableNames',{'runId','group','partition','segments','chunks'})]; %#ok<AGROW>
    end
end
save(fullfile(folder,'data.mat'),'sequences','segments','index','normalizer','manifest','addedManifest','runGroup','cfg','-v7.3');
writetable(trialSplit,fullfile(folder,'nov3_trial_split.csv'));writetable(coverage,fullfile(folder,'coverage.csv'));
writetable(index,fullfile(folder,'selection_anchors.csv'));writetable(manifest,fullfile(folder,'old_raw_sources.csv'));
writetable(addedManifest,fullfile(folder,'nov3_sources.csv'));
disp(trialSplit(trialSplit.partition=="selection",:));disp(coverage(coverage.group==2,:));
fprintf('FINE-TUNE DATA: old sequences/splits/normalizer preserved; %d old + %d Nov3 selection anchors.\n',height(old.index),height(index)-height(old.index));
end
