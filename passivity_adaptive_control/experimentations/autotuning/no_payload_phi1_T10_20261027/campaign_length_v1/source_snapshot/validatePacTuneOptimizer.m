function validatePacTuneOptimizer()
%VALIDATEPACTUNEOPTIMIZER Offline bounded/search/replay tests, no hardware.
cfg=pacTuneOptimizerConfig();cfg.gammaBounds=[.01 2;1e-6 .1;.05 10];
cfg.baselineGamma=[.4;.002;1];cfg.poolSize=400;cfg.maxObservations=20;
o=pacTuneOptimizerObservations();before=rng;
p=pacTuneOptimizerSuggest(cfg,o);assert(strcmp(p.reason,'baseline'));
assert(max(abs(log(p.gamma./cfg.baselineGamma)))<1e-12);
assert(isequal(rng,before),'Optimizer must not modify global RNG.');
pending=p.gamma';q=pacTuneOptimizerSuggest(cfg,o,pending);
assert(norm(log(q.gamma./p.gamma))>.01,'Pending baseline must not be reproposed.');
for k=1:9
    p=pacTuneOptimizerSuggest(cfg,o);g=p.gamma;
    assert(all(g>=cfg.gammaBounds(:,1)&g<=cfg.gammaBounds(:,2)));
    x=(log(g)-log(cfg.gammaBounds(:,1)))./log(cfg.gammaBounds(:,2)./cfg.gammaBounds(:,1));
    o(end+1,:)={string(k),g(1),g(2),g(3),sum((x-[.3;.4;.3]).^2),true,true,sum(x)<1.6,8}; %#ok<AGROW>
end
% Preserve repeated raw measurements; GP groups their mean and observed variance.
repeat=o(1,:);repeat.candidateId="repeat";repeat.objective=repeat.objective+.05;o=[o;repeat];
fault=o(2,:);fault.candidateId="fault";fault.objective=NaN;fault.usableObjective=false;fault.feasible=false;o=[o;fault];
missing=fault;missing.candidateId="missing";missing.constraintKnown=false;o=[o;missing];
p=pacTuneOptimizerSuggest(cfg,o);
assert(isfinite(p.predictedMean)&&isfinite(p.predictedStd)&&p.predictedStd>=0);
assert(p.feasibilityScore>0&&p.feasibilityScore<1);
assert(p.training.maximumReplicates>=2&&p.training.excludedObjectives==2);
assert(p.training.knownConstraints==height(o)-1);
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
save(fullfile(folder,'restart.mat'),'cfg','o','p');
s=load(fullfile(folder,'restart.mat'));q=pacTuneOptimizerSuggest(s.cfg,s.o);
assert(isequaln(p,q),'Reloaded configuration/observations must reproduce proposal.');
cfg.maxObservations=height(o);q=pacTuneOptimizerSuggest(cfg,o);assert(strcmp(q.action,'stop'));
bad=cfg;bad.gammaBounds(1,1)=0;expect(@()pacTuneOptimizerSuggest(bad,o),'pacTune:OptimizerBounds');
bad=cfg;bad.baselineGamma(1)=100;expect(@()pacTuneOptimizerSuggest(bad,o),'pacTune:OptimizerBounds');
badData=o;badData.objective(1)=NaN;
expect(@()pacTuneOptimizerSuggest(cfg,badData),'pacTune:OptimizerData');
fprintf('PASS: optimizer bounds/log baseline, pending avoidance, RNG, GP uncertainty, repeats, failures, restart and budget.\n');
end

function expect(fn,id)
failed=false;
try,fn();catch err,failed=strcmp(err.identifier,id);if ~failed,rethrow(err);end,end
assert(failed,'Expected %s.',id);
end
