function prepareForceResidualData(root,folder,cfg)
%PREPAREFORCERESIDUALDATA Old two days only; whole-trial internal holdouts.
assert(~isfolder(folder));mkdir(folder);
old=fullfile(root,'residual_dynamics_learning','datasets','two_day_v2');
P=load(fullfile(old,'preparation.mat'));manifest=P.manifest;
assert(height(manifest)==10&&all(ismember(manifest.cohort,{'day1_root','day2_20261027'})));
sequences=cell(height(manifest),1);normalizerRows=[];index=table;coverage=table;
for r=1:height(manifest)
    id=char(manifest.runId(r));source=char(manifest.source(r));
    assert(strcmp(filehash(source),manifest.sha256(r)),'Old raw recording changed.');
    a=load(fullfile(old,id,'sequences.mat'));S=struct;
    for f={'t','physicalTime','dt','q','v','p','trial','channel','phase','valid','segment','mature'}
        S.(f{1})=a.S.(f{1});
    end
    S.runId=id;S.split=repmat("excluded",numel(S.t),1);positive=S.trial>0;
    S.split(positive)="train";
    for channel=unique(S.channel(positive))'
        trials=unique(S.trial(positive&S.channel==channel));assert(numel(trials)>=2);
        held=trials(ceil(numel(trials)/2));S.split(positive&S.channel==channel&S.trial==held)="selection";
    end
    use=S.mature&S.split=="train";
    normalizerRows=[normalizerRows;S.q(use,:) S.v(use,:) S.p(use,:) S.dt(use)]; %#ok<AGROW>
    counts=zeros(1,2);
    for partition={'train','selection'}
        part=partition{1};K=cfg.trainSteps;cap=cfg.maxTrainWindowsPerRun;
        if strcmp(part,'selection'),K=cfg.selectionSteps;cap=cfg.maxSelectionWindowsPerRun;end
        anchors=[];
        for k=cfg.historyFrames+1:numel(S.t)-K
            if S.split(k)~=string(part),continue;end
            rows=k-cfg.historyFrames:k+K;
            if all(S.mature(rows))&&S.segment(k)>0&&all(S.segment(rows)==S.segment(k))&&all(S.split(rows)==S.split(k))
                anchors(end+1)=k; %#ok<AGROW>
            end
        end
        if isempty(anchors),continue;end
        chosen=anchors(unique(round(linspace(1,numel(anchors),min(cap,numel(anchors))))));
        index=[index;table(repmat(r,numel(chosen),1),chosen(:),repmat(string(part),numel(chosen),1), ...
            S.trial(chosen),repmat(string(id),numel(chosen),1),'VariableNames', ...
            {'run','anchor','partition','trial','runId'})]; %#ok<AGROW>
        counts(1+strcmp(part,'selection'))=numel(chosen);
    end
    sequences{r}=S;
    coverage=[coverage;table(string(id),nnz(S.mature&S.split=="train"),counts(1),counts(2), ...
        'VariableNames',{'runId','matureTrainingRows','trainingWindows','selectionWindows'})]; %#ok<AGROW>
end
assert(all(isfinite(normalizerRows),'all'));
normalizer=struct('mean',mean(normalizerRows,1)','std',std(normalizerRows,0,1)', ...
    'columns',{{'q2','q3','v2','v3','p1','p2','p3','previousPhysicalDt'}}, ...
    'scope','October19+October27 mature TRAIN rows only; November2 never loaded');
normalizer.std(normalizer.std<1e-10)=1;
save(fullfile(folder,'training_data.mat'),'sequences','index','normalizer','cfg','manifest','-v7.3');
writetable(index,fullfile(folder,'window_manifest.csv'));writetable(coverage,fullfile(folder,'coverage.csv'));
writetable(manifest,fullfile(folder,'raw_sources.csv'));
disp(coverage);fprintf('Prepared force model: %d train / %d internal-selection windows; old days only.\n',nnz(index.partition=="train"),nnz(index.partition=="selection"));
end
function h=filehash(path)
f=fopen(path,'rb');assert(f>=0);c=onCleanup(@()fclose(f)); %#ok<NASGU>
bytes=fread(f,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
