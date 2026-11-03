function h=sgWarm(params,S,seg,anchor,cfg,n)
%SGWARM Observed prefix only; excludes the anchor's current input.
h=sgZero(cfg,1);
first=find(S.physicalTime>=S.physicalTime(seg.starts(1))-cfg.burnIn_s,1);
first=max(first,seg.firstMature);
assert(first<seg.starts(1)&&anchor>=seg.starts(1));
for k=first:anchor-1
    x=[S.q(k,:) S.v(k,:) S.p(k,:) S.dt(k)]';
    h=sgCell(params,(x-n.mean)./n.std,h,cfg);
end
end
