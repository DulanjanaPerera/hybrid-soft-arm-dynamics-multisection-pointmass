function B=residualBatch(S,anchors,task,cfg,normalizer)
%RESIDUALBATCH Materialize explicit leakage-checked windows for MLP/GRU.
% Xsequence: features x history x examples. No test statistics fitted here.
H=cfg.historyFrames;F=cfg.previewFrames;anchors=anchors(:);
B.Xsequence=zeros(8,H,numel(anchors));B.preview=zeros(2,F,numel(anchors));
B.Y=zeros(numel(anchors),2);B.p0=zeros(numel(anchors),1);B.dtNext=B.p0;
for i=1:numel(anchors)
    k=anchors(i);rows=k-H+1:k;allRows=k-H:k+F;
    assert(all(S.mature(allRows))&&all(S.segment(allRows)==S.segment(k))&& ...
        all(S.split(allRows)==S.split(k)),'Window crosses invalid gap or split.');
    if strcmp(task,'inverse')
        assert(ismember(k,S.inverseAnchors));pressureRows=rows-1;
        B.Y(i,:)=S.inverseResidualCandidate(k,:);
        B.preview(:,:,i)=S.q(k+1:k+F,:)';B.p0(i)=S.p0(k);
    else
        assert(strcmp(task,'forward')&&ismember(k,S.forwardAnchors));pressureRows=rows;
        B.Y(i,:)=S.q(k+1,:);
        B.preview(:,:,i)=NaN; % Unknown future response is never forward input.
    end
    values=[S.q(rows,:) S.v(rows,:) S.p(pressureRows,:)];
    B.Xsequence(:,:,i)=[(values-normalizer.mean)./normalizer.std S.dt(rows)]';
    B.dtNext(i)=S.physicalTime(k+1)-S.physicalTime(k);
end
B.Xflat=reshape(B.Xsequence,8*H,[])';B.anchors=anchors;
B.split=S.split(anchors);B.runId=S.runId;B.task=task;
% Inverse preview represents an OFFLINE motion requirement. Replace it with
% planned desired motion for deployment; don't access future measurements.
end
