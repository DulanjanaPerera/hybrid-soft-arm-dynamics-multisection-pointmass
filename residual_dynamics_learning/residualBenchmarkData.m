function D=residualBenchmarkData(folder,task,partition,stride)
%RESIDUALBENCHMARKDATA Explicit partition gate; no test data in selection.
d=load(fullfile(folder,'preparation.mat'));n=load(fullfile(folder,'normalizer.mat'));
if isfield(d.cfg,'boundMode')
    assert(d.cfg.boundMode==0,'Legacy screening requires coordinate-bounds-v1 data/context.');
end
D=struct('Seq',[],'Context',[],'Y',[],'Base',[],'Truth',[],'p0',[], ...
    'Previous',[],'Current',[],'Velocity',[],'dt',[],'Run',strings(0,1),'Anchor',[]);
for i=1:height(d.summary)
    external=strcmp(d.summary.cohort{i},'day2_20261027');
    if external~=strcmp(partition,'external_test'),continue;end
    id=d.summary.runId{i};a=load(fullfile(folder,id,'sequences.mat'));S=a.S;
    if strcmp(task,'inverse'),k=S.inverseAnchors;else,k=S.forwardAnchors;end
    k=k(S.split(k)==string(partition));k=k(1:stride:end);
    B=residualBatch(S,k,task,a.cfg,n.normalizer);N=numel(k);
    if strcmp(task,'inverse')
        ctx=[reshape(B.preview,[],N)' B.p0 S.inverseDifferentialNominal(k,:)];
        base=S.inverseDifferentialNominal(k,:);truth=S.p(k,:);Y=B.Y;
        prev=S.d(k-1,:);
    else
        base=zeros(N,2);
        for j=1:N
            z=residualFrozenStep(S.q(k(j),:)',S.v(k(j),:)',S.p(k(j),:)',B.dtNext(j));
            base(j,:)=z(1:2)';
        end
        ctx=[B.dtNext base];truth=B.Y;Y=1000*(truth-base);prev=S.q(k,:);
    end
    D.Seq=cat(3,D.Seq,single(B.Xsequence));D.Context=[D.Context;ctx];
    D.Y=[D.Y;Y];D.Base=[D.Base;base];D.Truth=[D.Truth;truth];D.Previous=[D.Previous;prev];
    D.p0=[D.p0;B.p0];D.Current=[D.Current;S.q(k,:)];D.Velocity=[D.Velocity;S.v(k,:)];
    D.dt=[D.dt;B.dtNext];D.Run=[D.Run;repmat(string(id),N,1)];D.Anchor=[D.Anchor;k];
end
assert(~isempty(D.Y));
end
