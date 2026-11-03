function verifyResidualRecordedMaps(folder)
%VERIFYRESIDUALRECORDEDMAPS Verify model force/pressure conventions in logs.
% Not a predictive test against real measured motion or historical byte proof.
d=load(fullfile(folder,'configuration.mat'));M=d.manifest;report=table;
alpha=pi*.0065^2*1e5;A=[-1 1 0;-1 0 1];
for k=1:height(M)
    raw=load(char(M.source(k)),'out');ds=raw.out.logsout;
    [t,p]=sig(ds,'pressureUsed_bar');[tt,tau]=sig(ds,'tau_N');assert(isequal(t,tt));
    good=all(isfinite([p tau]),2);e=max(abs(p(good,:)*(alpha*A)'-tau(good,:)),[],'all');
    [t,p]=sig(ds,'pressureDesired_bar');[tt,tau]=sig(ds,'tauAvailable_N');assert(isequal(t,tt));
    ok=all(isfinite([p tau]),2);ei=max(abs(p(ok,:)*(alpha*A)'-tau(ok,:)),[],'all');
    e8=max(abs(max(p(ok,:)-.8,0)*(alpha*A)'-tau(ok,:)),[],'all');
    activation=0;best=ei;if e8<ei,activation=.8;best=e8;end
    report=[report;table(M.runId(k),e,ei,e8,activation,best,nnz(good),nnz(ok),'VariableNames', ...
        {'runId','forwardEffectivePressureMapError_N','inverseZeroActivationError_N', ...
        'inversePoint8ActivationError_N','identifiedInverseActivation_bar', ...
        'bestInverseMapError_N','finiteForward','finiteInverse'})]; %#ok<AGROW>
end
writetable(report,fullfile(folder,'recorded_map_verification.csv'));
assert(all(report.forwardEffectivePressureMapError_N<1e-6&report.bestInverseMapError_N<1e-6));
fprintf('PASS: recorded effective-pressure maps verified; historical inverse activation convention recorded per run (not assumed identical).\n');
end
function [t,x]=sig(ds,name)
i=find(strcmp(ds.getElementNames,name));assert(numel(i)==1);
v=ds{i}.Values;t=double(v.Time(:));x=double(v.Data);
if v.IsTimeFirst,x=reshape(x,numel(t),[]);else,x=reshape(x,[],numel(t))';end
end
