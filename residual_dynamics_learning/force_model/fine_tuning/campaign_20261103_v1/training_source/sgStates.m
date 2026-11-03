function H=sgStates(params,D,index,cfg)
%SGSTATES Replay measured inputs once, capture PRE-current states at anchors.
H=sgZero(cfg,height(index));
for sid=unique(index.segmentId)'
    seg=D.segments(sid);S=D.sequences{seg.run};ix=find(index.segmentId==sid);
    [anchors,order]=sort(index.anchor(ix));ix=ix(order);
    h=sgWarm(params,S,seg,seg.starts(1),cfg,D.normalizer);j=1;
    for k=seg.starts(1):anchors(end)
        while j<=numel(anchors)&&anchors(j)==k
            for l=1:numel(h),H{l}(:,ix(j))=h{l};end
            j=j+1;
        end
        x=[S.q(k,:) S.v(k,:) S.p(k,:) S.dt(k)]';
        h=sgCell(params,(x-D.normalizer.mean)./D.normalizer.std,h,cfg);
    end
end
end
