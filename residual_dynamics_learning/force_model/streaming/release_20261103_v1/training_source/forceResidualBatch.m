function B=forceResidualBatch(sequences,index,cfg,K)
%FORCERESIDUALBATCH Raw causal history and separate future length targets.
N=height(index);H=cfg.historyFrames;
B=struct('history',zeros(8,N,H),'q0',zeros(2,N),'v0',zeros(2,N), ...
    'pressure',zeros(3,N,K),'dt',zeros(1,N,K),'previousDt',zeros(1,N,K), ...
    'truth',zeros(2,N,K));
for i=1:N
    S=sequences{index.run(i)};k=index.anchor(i);rows=k-H:k+K;
    assert(all(S.mature(rows))&&all(S.segment(rows)==S.segment(k))&& ...
        all(S.split(rows)==S.split(k)),'Window crosses a gap or partition.');
    history=k-H+1:k;B.history(:,i,:)=[S.q(history,:) S.v(history,:) S.p(history,:) S.dt(history)]';
    B.q0(:,i)=S.q(k,:)';B.v0(:,i)=S.v(k,:)';
    B.pressure(:,i,:)=S.p(k:k+K-1,:)';B.dt(:,i,:)=S.dt(k+1:k+K)';
    B.previousDt(:,i,:)=S.dt(k:k+K-1)';B.truth(:,i,:)=S.q(k+1:k+K,:)';
end
end
