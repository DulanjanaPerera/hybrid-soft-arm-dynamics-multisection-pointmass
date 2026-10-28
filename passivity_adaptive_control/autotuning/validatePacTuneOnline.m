function validatePacTuneOnline()
%VALIDATEPACTUNEONLINE Synthetic only. Does not load or run a model.
c=pacTuneOnlineConfig();c.optimizer.gammaBounds=[.1 1;.0005 .008;1 10];
c.optimizer.baselineGamma=[.5;.002;5];c.optimizer.poolSize=40;
c.optimizer.maxObservations=3;c.settle_s=1;c.circles=3;c.scoreCircles=3;
u=struct('t',0,'healthy',true,'holdReady',true,'appliedId',0, ...
    'appliedGamma',zeros(3,1),'scoreReady',false,'scoreId',0, ...
    'metrics',[1 1 .003 .02],'coverage',1,'limitedFraction',0);
[s,a]=pacTuneOnlineStep([],u,c);assert(a.apply&&a.id==1);
u.appliedId=1;u.appliedGamma=a.gamma;u.t=.1;
[s,a]=pacTuneOnlineStep(s,u,c);assert(strcmp(s.phase,'settle'));
u.t=1.1;[s,a]=pacTuneOnlineStep(s,u,c);assert(strcmp(a.mode,'run'));
u.t=16.1;[s,a]=pacTuneOnlineStep(s,u,c);assert(strcmp(s.phase,'score'));
u.scoreReady=true;u.scoreId=1;u.t=16.2;
[s,a]=pacTuneOnlineStep(s,u,c);assert(height(s.observations)==1);
assert(strcmp(s.phase,'hold')&&s.consecutivePasses==1);
u.t=16.3;[s,a]=pacTuneOnlineStep(s,u,c);assert(a.id==2&&a.apply);
u.appliedId=2;u.t=16.4;[s,a]=pacTuneOnlineStep(s,u,c);
u.t=17.4;[s,a]=pacTuneOnlineStep(s,u,c);
u.t=32.4;[s,a]=pacTuneOnlineStep(s,u,c);
u.scoreId=2;u.t=32.5;[s,a]=pacTuneOnlineStep(s,u,c);
assert(strcmp(a.mode,'stop')&&strcmp(s.reason,'targets_repeated_operationally'));
% Stale scores never advance the ledger; invalid coverage stops the session.
q=s;q.phase='score';q.since=u.t;q.id=3;
u.scoreId=2;[r,~]=pacTuneOnlineStep(q,u,c);assert(height(r.observations)==2);
u.scoreId=3;u.coverage=.5;[r,a]=pacTuneOnlineStep(q,u,c);
assert(strcmp(a.mode,'stop')&&strcmp(r.reason,'invalid_score'));
u.coverage=1;u.limitedFraction=.01;[r,~]=pacTuneOnlineStep(q,u,c);
assert(strcmp(r.reason,'pressure_limiting'));
u.limitedFraction=0;u.healthy=false;[r,~]=pacTuneOnlineStep(q,u,c);
assert(strcmp(r.reason,'unhealthy_or_clock_reversal'));
% Infeasible accuracy asks BO for a bounded local candidate, not a reset.
u.healthy=true;u.metrics=[3 3 .003 .02];q.observations=q.observations([],:);
q.consecutivePasses=0;[r,a]=pacTuneOnlineStep(q,u,c);
assert(strcmp(r.phase,'hold')&&height(r.observations)==1);
assert(all(r.gamma>=q.gamma/c.maxGammaFactor & r.gamma<=q.gamma*c.maxGammaFactor));
assert(a.preserveAdaptiveState);
% Global budget and apply acknowledgement watchdog.
c.optimizer.maxObservations=1;[r,~]=pacTuneOnlineStep(q,u,c);
assert(strcmp(r.reason,'budget_exhausted'));
u.t=0;u.scoreReady=false;[q,~]=pacTuneOnlineStep([],u,c);
u.appliedId=-1;u.t=11;[r,~]=pacTuneOnlineStep(q,u,c);
assert(strcmp(r.reason,'apply_timeout'));
fprintf('PASS: offline online supervisor handshake, repeats, local BO, stale scores, constraints, watchdog and budget.\n');
end
