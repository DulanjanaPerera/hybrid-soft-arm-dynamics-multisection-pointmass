function B=sgBatch(D,jobs,cfg)
%SGBATCH Variable-length ordered chunks, safe padding, masks, no gap crossing.
N=size(jobs,1);lastP=zeros(N,1);lastA=lastP;first=lastP;
for i=1:N
    seg=D.segments(jobs(i,1));S=D.sequences{seg.run};block=jobs(i,2);
    k=seg.starts(block);first(i)=k;lastP(i)=seg.ends(block);
    lastA(i)=find(S.physicalTime<=S.physicalTime(k)+cfg.horizon_s+1e-10,1,'last');
    assert(lastA(i)<=seg.auxEnds(block)&&lastA(i)>k);
end
P=max(lastP-first);A=max(lastA-first);
B=struct('x',zeros(8,N,P),'measuredQ',zeros(2,N,P),'measuredV',zeros(2,N,P), ...
    'pressure',zeros(3,N,P),'dt',.05*ones(1,N,P),'truth',zeros(2,N,P),'mask',zeros(1,N,P), ...
    'q0',zeros(2,N),'v0',zeros(2,N),'auxP',zeros(3,N,A),'auxDt',.05*ones(1,N,A), ...
    'auxPrevDt',.05*ones(1,N,A),'auxTruth',zeros(2,N,A),'auxMask',zeros(1,N,A));
for i=1:N
    seg=D.segments(jobs(i,1));S=D.sequences{seg.run};k=first(i);np=lastP(i)-k;na=lastA(i)-k;
    rows=k:lastA(i);assert(all(S.mature(rows))&&all(S.segment(rows)==seg.sourceSegment));
    assert(all(S.split(rows)==string(seg.partition)));
    B.x(:,i,1:np)=[S.q(k:k+np-1,:) S.v(k:k+np-1,:) S.p(k:k+np-1,:) S.dt(k:k+np-1)]';
    B.measuredQ(:,i,1:np)=S.q(k:k+np-1,:)';B.measuredV(:,i,1:np)=S.v(k:k+np-1,:)';
    B.pressure(:,i,1:np)=S.p(k:k+np-1,:)';B.dt(:,i,1:np)=S.dt(k+1:k+np)';
    B.truth(:,i,1:np)=S.q(k+1:k+np,:)';B.mask(:,i,1:np)=1;
    B.q0(:,i)=S.q(k,:)';B.v0(:,i)=S.v(k,:)';
    B.auxP(:,i,1:na)=S.p(k:k+na-1,:)';B.auxDt(:,i,1:na)=S.dt(k+1:k+na)';
    B.auxPrevDt(:,i,1:na)=S.dt(k:k+na-1)';B.auxTruth(:,i,1:na)=S.q(k+1:k+na,:)';B.auxMask(:,i,1:na)=1;
end
end
