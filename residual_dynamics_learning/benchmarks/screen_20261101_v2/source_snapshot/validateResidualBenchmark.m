function validateResidualBenchmark(dataset)
%VALIDATERESIDUALBENCHMARK Feature and nominal integration integrity.
validateResidualNominal();
q=[.002;-.003];v=[.001;-.001];p=[.2;.7;1];
a=residualNominal(q,v,[0;0],p);h=1e-7;z=residualFrozenStep(q,v,p,h);
assert(norm((z-[q;v])/h-[v;a])<1e-3,'Frozen flow has wrong local dynamics.');
d=load(fullfile(dataset,'preparation.mat'));n=load(fullfile(dataset,'normalizer.mat'));
s=load(fullfile(dataset,d.summary.runId{1},'sequences.mat'));S=s.S;k=S.inverseAnchors(50);
b=residualBatch(S,k,'inverse',d.cfg,n.normalizer);
assert(all(isfinite(b.Xsequence(:)))&&all(isfinite(b.preview(:))));
changed=S;changed.p(k,:)=100;
c=residualBatch(changed,k,'inverse',d.cfg,n.normalizer);assert(isequal(b.Xsequence,c.Xsequence));
fprintf('PASS: frozen nominal local derivative, finite inverse preview and no current-pressure history leak.\n');
end
