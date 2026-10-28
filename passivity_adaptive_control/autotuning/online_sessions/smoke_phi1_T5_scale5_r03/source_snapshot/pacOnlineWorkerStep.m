function s=pacOnlineWorkerStep(s,tx,c,rt)
%PACONLINEWORKERSTEP Pure ask/tell worker, independent of network and model.
if isempty(s)
    s=struct('observations',pacTuneOptimizerObservations(),'lastCompleted',0, ...
        'requestId',1,'gamma',rt.baselineGamma,'action',1, ...
        'passes',0,'reason','baseline','events',{{}});
end
if numel(tx)<16||tx(1)~=rt.session||tx(4)~=1||tx(3)==0,return;end
id=tx(3);
if id<=s.lastCompleted,return;end
if id~=s.requestId||tx(2)~=id+1||any(abs(tx(14:16)-s.gamma)>1e-10)
    s.action=4;s.reason='candidate_or_gain_mismatch';return;
end
m=tx(7:10)';coverage=tx(11);limited=tx(12);
if any(~isfinite(m))||any(m<0)||~isfinite(coverage)|| ...
        coverage<c.minimumCoverage||coverage>1||~isfinite(limited)||limited~=0
    s.action=4;s.reason='invalid_or_limited_batch';return;
end
pass=all(m<=c.targets);
s.observations(end+1,:)={string(id),s.gamma(1),s.gamma(2),s.gamma(3), ...
    m(4)/c.targets(4),true,true,pass,rt.period_s*rt.circles};
s.events{end+1}=struct('id',id,'gamma',s.gamma,'metrics',m,'coverage',coverage, ...
    'time',tx(5),'protocol','preserved-state');
s.lastCompleted=id;s.requestId=id+1;
if pass,s.passes=s.passes+1;else,s.passes=0;end
if s.passes>=c.confirmations
    s.action=2;s.reason='operational_targets_repeated';
elseif height(s.observations)>=c.optimizer.maxObservations
    s.action=3;s.reason='budget_exhausted';
elseif pass
    s.reason='repeat_same_gains';
else
    oc=c.optimizer;
    oc.gammaBounds(:,1)=max(oc.gammaBounds(:,1),s.gamma/rt.maxGammaFactor);
    oc.gammaBounds(:,2)=min(oc.gammaBounds(:,2),s.gamma*rt.maxGammaFactor);
    oc.baselineGamma=s.gamma;
    o=s.observations;g=o{:,{'gammaK','gammaD','gammaH'}};
    inside=all(g>=oc.gammaBounds(:,1)'&g<=oc.gammaBounds(:,2)',2);
    p=pacTuneOptimizerSuggest(oc,o(inside,:));
    if strcmp(p.action,'propose'),s.gamma=p.gamma;s.reason=p.reason;
    else,s.action=4;s.reason=p.reason;end
end
end
