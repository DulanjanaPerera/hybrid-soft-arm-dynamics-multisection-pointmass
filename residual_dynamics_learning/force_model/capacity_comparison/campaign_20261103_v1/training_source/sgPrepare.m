function sgPrepare(dataFolder,folder,cfg)
%SGPREPARE Common anchors, ordered chunks and physical time caps. OLD days only.
assert(~isfolder(folder));mkdir(folder);D=load(fullfile(dataFolder,'training_data.mat'));
assert(all(ismember(D.manifest.cohort,{'day1_root','day2_20261027'})));
sequences=D.sequences;normalizer=D.normalizer;manifest=D.manifest;
segments=struct([]);index=table;coverage=table;
for run=1:numel(sequences)
    S=sequences{run};allAnchors=[];segIds=[];before=numel(segments);
    for segment=unique(S.segment(S.mature))'
        rows=find(S.segment==segment&S.mature);if isempty(rows),continue;end
        first=rows(1);last=rows(end);
        k=find(S.physicalTime>=S.physicalTime(first)+cfg.commonBurnIn_s,1);
        if isempty(k)||k>=last||S.segment(k)~=segment,continue;end
        starts=[];ends=[];auxEnds=[];
        while S.physicalTime(last)-S.physicalTime(k)>=cfg.maxHorizon_s
            endPrimary=find(S.physicalTime<=S.physicalTime(k)+cfg.chunk_s+1e-10,1,'last');
            endAux=find(S.physicalTime<=S.physicalTime(k)+cfg.maxHorizon_s+1e-10,1,'last');
            if endPrimary<=k||endAux>last,break;end
            starts(end+1)=k;ends(end+1)=endPrimary;auxEnds(end+1)=endAux; %#ok<AGROW>
            k=endPrimary;
        end
        if isempty(starts),continue;end
        part=S.split(starts(1));assert(part=="train"||part=="selection");
        seg=struct('run',run,'sourceSegment',segment,'firstMature',first,'last',last, ...
            'starts',starts,'ends',ends,'auxEnds',auxEnds,'partition',char(part));
        segments=[segments;seg];sid=numel(segments); %#ok<AGROW>
        allAnchors=[allAnchors starts];segIds=[segIds repmat(sid,size(starts))]; %#ok<AGROW>
    end
    for part={'train','selection'}
        part=part{1};take=find(arrayfun(@(s)strcmp(s.partition,part),segments(before+1:end)))+before;
        total=0;for s=take',total=total+numel(segments(s).starts);end
        coverage=[coverage;table(string(S.runId),string(part),numel(take),total, ...
            'VariableNames',{'runId','partition','segments','chunks'})]; %#ok<AGROW>
    end
    if isempty(allAnchors),continue;end
    keep=S.split(allAnchors)=="selection";a=allAnchors(keep);ids=segIds(keep);
    take=unique(round(linspace(1,numel(a),min(cfg.selectionAnchorsPerRun,numel(a)))));
    index=[index;table(repmat(run,numel(take),1),a(take)',ids(take)', ...
        'VariableNames',{'run','anchor','segmentId'})]; %#ok<AGROW>
end
assert(~isempty(index));
save(fullfile(folder,'data.mat'),'sequences','segments','normalizer','manifest','index','cfg','-v7.3');
writetable(index,fullfile(folder,'selection_anchors.csv'));writetable(coverage,fullfile(folder,'coverage.csv'));
writetable(manifest,fullfile(folder,'old_raw_sources.csv'));
disp(coverage);fprintf('Streaming preparation: %d segments, %d common selection anchors.\n',numel(segments),height(index));
end
