function suggestion=pacTuneOptimizerSuggest(cfg,observations,pending)
%PACTUNEOPTIMIZERSUGGEST Pure deterministic bounded GP ask/tell proposal.
% Does not evaluate an objective, save state, run hardware or change global RNG.
if nargin<3,pending=zeros(0,3);end
validate(cfg,observations,pending);
lo=log(cfg.gammaBounds(:,1))';span=log(cfg.gammaBounds(:,2))'-lo;
span(span==0)=1;
fixed=cfg.gammaBounds(:,1)==cfg.gammaBounds(:,2);
normalize=@(g)(log(g)-lo)./span;
denormalize=@(x)min(max(exp(lo+x.*span),cfg.gammaBounds(:,1)'),cfg.gammaBounds(:,2)');
G=observations{:,{'gammaK','gammaD','gammaH'}};
X=normalize(G);P=normalize(pending);seen=[X;P];
suggestion=struct('schemaVersion',cfg.schemaVersion,'action','propose','reason','', ...
    'gamma',NaN(3,1),'predictedMean',NaN,'predictedStd',NaN, ...
    'feasibilityScore',NaN,'acquisition',NaN,'config',cfg);
suggestion.training=struct('observations',height(observations), ...
    'usableObjectives',nnz(observations.usableObjective), ...
    'knownConstraints',nnz(observations.constraintKnown), ...
    'excludedObjectives',nnz(~observations.usableObjective),'uniqueObjectiveLocations',0, ...
    'maximumReplicates',0,'knownDurationObservations',nnz(isfinite(observations.elapsed_s)), ...
    'totalRecordedDuration_s',sum(observations.elapsed_s,'omitnan'));
if height(observations)>0&&~any(isfinite(observations.elapsed_s))
    suggestion.training.totalRecordedDuration_s=NaN;
end
suggestion.replicates=table;
if any(observations.usableObjective)
    use=observations.usableObjective;[g,~,groups]=unique(G(use,:),'rows');
    y=observations.objective(use);n=accumarray(groups,1);
    means=accumarray(groups,y,[],@mean);sd=accumarray(groups,y,[],@(a)std(a,0));
    suggestion.replicates=table(g(:,1),g(:,2),g(:,3),n,means,sd,sd./sqrt(n), ...
        'VariableNames',{'gammaK','gammaD','gammaH','count','meanObjective','observedStd','standardError'});
    suggestion.training.uniqueObjectiveLocations=size(g,1);
    suggestion.training.maximumReplicates=max(n);
end
if height(observations)>=cfg.maxObservations
    suggestion.action='stop';suggestion.reason='observation_budget_exhausted';return
end
% A local stream makes replay independent of the caller's global random stream.
stream=RandStream('mt19937ar','Seed',cfg.seed);
design=lhs(stream,max(1,cfg.initialDesignSize-1),3);design(:,fixed)=0;
design=[normalize(cfg.baselineGamma(:)');design];
for k=1:cfg.initialDesignSize
    if nearest(design(k,:),seen)>cfg.minimumDistance
        suggestion.gamma=denormalize(design(k,:))';
        if k==1
            suggestion.gamma=cfg.baselineGamma(:);suggestion.reason='baseline';
        else
            suggestion.reason='initial_log_space_design';
        end
        return
    end
end
pool=lhs(stream,cfg.poolSize,3);pool(:,fixed)=0;
pool=[pool;zeros(1,3);double(~fixed)'];
distance=nearest(pool,seen);pool=pool(distance>cfg.minimumDistance,:);
distance=distance(distance>cfg.minimumDistance);
if isempty(pool),suggestion.action='stop';suggestion.reason='search_space_exhausted';return;end
known=observations.constraintKnown;
% Beta(1,1) kernel-weighted frequency: a heuristic score, not a safety probability.
feasibility=.5*ones(size(pool,1),1);
if any(known)
    w=kernel(pool,X(known,:),cfg.feasibilityLengthScale);
    feasibility=(1+w*double(observations.feasible(known)))./(2+sum(w,2));
end
use=observations.usableObjective;
mu=NaN(size(pool,1),1);sigma=mu;
if any(use)
    y=observations.objective(use);center=mean(y);scale=max(std(y),.1);
    [locations,~,groups]=unique(X(use,:),'rows');
    normalized=(y-center)/scale;
    means=accumarray(groups,normalized,[],@mean);
    counts=accumarray(groups,1);
    variance=accumarray(groups,normalized,[],@(a)var(a,0));
    noise=max(variance,cfg.noiseFloor^2)./counts;
    K=kernel(locations,locations,cfg.lengthScale)+diag(noise+cfg.jitter);
    L=chol(K,'lower');cross=kernel(pool,locations,cfg.lengthScale);
    m=cross*(L'\(L\means));v=L\cross';
    sd=sqrt(max(1-sum(v.^2,1)',0));
    mu=center+scale*m;sigma=scale*sd;
    feasibleUse=use&known&observations.feasible;
    if any(feasibleUse)
        incumbent=(min(observations.objective(feasibleUse))-center)/scale;
        improvement=incumbent-m-cfg.improvementOffset;
        z=improvement./max(sd,1e-12);
        ei=improvement.*(.5*erfc(-z/sqrt(2)))+sd.*exp(-.5*z.^2)/sqrt(2*pi);
        acq=max(ei,0).*feasibility;
        suggestion.reason='feasibility_weighted_expected_improvement';
    else
        acq=feasibility.*(1+sd);
        suggestion.reason='seek_feasible_region';
    end
    if max(acq)<eps
        acq=feasibility.*(sd+.1*distance);suggestion.reason='numerical_ei_exploration_fallback';
    end
else
    acq=feasibility.*(1+distance);
    suggestion.reason='explore_without_usable_objective';
end
[~,j]=max(acq);
suggestion.gamma=denormalize(pool(j,:))';
suggestion.predictedMean=mu(j);suggestion.predictedStd=sigma(j);
suggestion.feasibilityScore=feasibility(j);suggestion.acquisition=acq(j);
end

function validate(c,o,p)
assert(strcmp(c.schemaVersion,'pac-gp-1'),'pacTune:OptimizerSchema','Unsupported optimizer.');
assert(isequal(size(c.gammaBounds),[3 2])&&all(isfinite(c.gammaBounds(:)))&& ...
    all(c.gammaBounds(:,1)>0)&&all(c.gammaBounds(:,2)>=c.gammaBounds(:,1)), ...
    'pacTune:OptimizerBounds','Positive finite 3x2 gamma bounds required.');
assert(numel(c.baselineGamma)==3&&all(isfinite(c.baselineGamma(:)))&& ...
    all(c.baselineGamma(:)>=c.gammaBounds(:,1)&c.baselineGamma(:)<=c.gammaBounds(:,2)), ...
    'pacTune:OptimizerBounds','Baseline must be inside bounds.');
for f={'seed','initialDesignSize','poolSize','maxObservations'}
    x=c.(f{1});assert(isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x,'pacTune:OptimizerConfig','Invalid %s.',f{1});
end
for f={'lengthScale','noiseFloor','jitter','feasibilityLengthScale','minimumDistance'}
    x=c.(f{1});assert(isscalar(x)&&isfinite(x)&&x>0,'pacTune:OptimizerConfig','Invalid %s.',f{1});
end
assert(isscalar(c.improvementOffset)&&isfinite(c.improvementOffset)&&c.improvementOffset>=0);
assert(istable(o),'pacTune:OptimizerData','Observations must be a table.');
template=pacTuneOptimizerObservations();
assert(all(ismember(template.Properties.VariableNames,o.Properties.VariableNames)), ...
    'pacTune:OptimizerData','Missing observation fields.');
g=o{:,{'gammaK','gammaD','gammaH'}};
assert(all(isfinite(g(:)))&&all(all(g>=c.gammaBounds(:,1)'&g<=c.gammaBounds(:,2)')), ...
    'pacTune:OptimizerData','Observation gamma is outside bounds.');
assert(all(isfinite(o.objective(o.usableObjective))), ...
    'pacTune:OptimizerData','Usable objectives must be finite; missing values are not penalties.');
for f={'usableObjective','constraintKnown','feasible'}
    assert(islogical(o.(f{1}))&&iscolumn(o.(f{1})),'pacTune:OptimizerData','Flags must be logical columns.');
end
assert(all(isnan(o.elapsed_s)|(isfinite(o.elapsed_s)&o.elapsed_s>=0)), ...
    'pacTune:OptimizerData','Durations must be nonnegative or unknown (NaN).');
assert(size(p,2)==3&&all(isfinite(p(:)))&&all(all(p>=c.gammaBounds(:,1)'&p<=c.gammaBounds(:,2)')), ...
    'pacTune:OptimizerData','Pending gamma is outside bounds.');
end

function X=lhs(stream,n,d)
X=zeros(n,d);
for j=1:d,X(:,j)=(randperm(stream,n)'-rand(stream,n,1))/n;end
end

function K=kernel(a,b,ell)
distance=max(sum(a.^2,2)+sum(b.^2,2)'-2*a*b',0);
K=exp(-.5*distance/ell^2);
end

function d=nearest(a,b)
if isempty(b),d=Inf(size(a,1),1);else,d=sqrt(min(max(sum(a.^2,2)+sum(b.^2,2)'-2*a*b',0),[],2));end
end
