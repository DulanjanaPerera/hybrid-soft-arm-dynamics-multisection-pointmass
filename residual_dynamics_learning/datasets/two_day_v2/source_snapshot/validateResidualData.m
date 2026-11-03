function validateResidualData(folder)
%VALIDATERESIDUALDATA Real prepared-data isolation and label-leak checks.
d=load(fullfile(folder,'preparation.mat'));n=load(fullfile(folder,'normalizer.mat'));
train=zeros(0,7);total=0;target=zeros(0,2);
for i=1:height(d.summary)
    p=load(fullfile(folder,char(d.summary.runId(i)),'sequences.mat'));S=p.S;cfg=p.cfg;
    if strcmp(S.cohort,'day2_20261027')
        assert(~any(S.split=="train"|S.split=="validation"));
    else
        assert(~any(S.split=="external_test"));
    end
    for trial=unique(S.trial(S.trial>0))'
        assert(numel(unique(S.split(S.trial==trial)))==1,'A trial crosses partitions.');
    end
    use=S.mature&S.split=="train";
    train=[train;S.q(use,:) S.v(use,:) S.p(use,:)]; %#ok<AGROW>
    a=S.inverseAnchors(S.split(S.inverseAnchors)=="train");
    target=[target;S.inverseResidualCandidate(a,:)]; %#ok<AGROW>
    for ids={S.forwardAnchors,S.inverseAnchors}
        for k=ids{1}'
            rows=k-cfg.historyFrames:k+cfg.previewFrames;
            assert(all(S.mature(rows))&&S.segment(k)>0&& ...
                all(S.segment(rows)==S.segment(k))&&all(S.split(rows)==S.split(k)));
        end
    end
    if ~isempty(S.inverseAnchors)
        k=S.inverseAnchors(1);a=residualBatch(S,k,'inverse',cfg,n.normalizer);
        changed=S;changed.p(k,:)=changed.p(k,:)+100; % Current target command.
        b=residualBatch(changed,k,'inverse',cfg,n.normalizer);
        assert(isequaln(a.Xsequence,b.Xsequence),'Current target pressure leaked into inverse history.');
        f=residualBatch(S,k,'forward',cfg,n.normalizer);changed=S;
        changed.q(k+1:k+cfg.previewFrames,:)=1e6;
        g=residualBatch(changed,k,'forward',cfg,n.normalizer);
        assert(isequaln(f.Xsequence,g.Xsequence)&&all(isnan(f.preview(:))), ...
            'Future measured motion leaked into forward input.');
    end
    total=total+numel(S.forwardAnchors);
end
assert(max(abs(mean(train)-n.normalizer.mean))<1e-12);
assert(max(abs(mean(target)-n.normalizer.inverseResidualMean))<1e-12);
fprintf('PASS: %d windows; whole-trial/day isolation, gap guards, train-only scaling, inverse pressure and forward future-target leak checks.\n',total);
end
