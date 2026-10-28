function validatePacOnlineRuntime()
%VALIDATEPACONLINERUNTIME Synthetic plant/reference; no arm model or devices.
folder=tempname;cfile=pacOnlinePrepareSession(folder,[.0515;.0005999;1.49999], ...
    [.02575 .103;.00029995 .0011998;.749995 2.99998],5);
cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
d=load(cfile);c=d.rt;hc=d.c;
c.start_s=.2;c.stationaryDwell_s=.1;c.ramp_s=.3;c.spinRamp_s=.2;c.settle_s=.1;
c.circles=3;c.scoreCircles=2;c.period_s=1;hc.confirmations=2;
hc.targets=[100 100 1 10]; % Toy convergence fixture, not hardware limits.
s=pacOnlineInitialState(c);worker=[];rx=zeros(7,1);g=c.baselineGamma;
tx=zeros(16,1);samples=zeros(2000,8);count=0;seq=0;lastQ=zeros(2,1);
raw=zeros(2000,5);
for t=0:.05:25
    % Perfect, gently delayed synthetic following (q from previous reference).
    q=[.013*cos(s.theta)/2-sqrt(3)*.013*sin(s.theta)/2; ...
       .013*cos(s.theta)/2+sqrt(3)*.013*sin(s.theta)/2]*s.phi;
    dq=(q-lastQ)/.05;lastQ=q;
    actual=[g(1);g(1);g(1);g(2);g(2);g(3);g(3);g(3);1;1];
    u=[t;q;dq;t;1;3;0;1;1;1];
    worker=pacOnlineWorkerStep(worker,tx,hc,c);seq=seq+1;
    rx=[c.session;seq;worker.requestId;worker.action;worker.gamma];
    [s,y,tx]=pacOnlineRuntimeStep(s,u,actual,c.settings,rx,c);g=y(3:5);
    count=count+1;samples(count,:)=[t y(1) y(2) y(8) y(9) y(11) y(6) y(7)];
    qd=[.013*cos(y(1))/2-sqrt(3)*.013*sin(y(1))/2; ...
        .013*cos(y(1))/2+sqrt(3)*.013*sin(y(1))/2]*y(2);
    raw(count,:)=[t q' qd'];
    if s.phase==10,break;end
end
assert(worker.lastCompleted>=1,'No completed batch.');
assert(s.reason==22,'Synthetic target repeats did not complete normally.');
assert(all(isfinite(worker.observations.objective)));
assert(s.coverage>=.95&&s.coverage<=1);
assert(all(samples(1:count,3)>=-1e-10));
% Independent offline reconstruction verifies the online metric accumulator.
raw=raw(1:count,:);tt=raw(:,1);valid=true(count,1);
vq=pacTuneFixedVelocity(raw(:,2:3),tt,valid,.2,.075);
vd=pacTuneFixedVelocity(raw(:,4:5),tt,valid,.1,.075);
run=samples(1:count,4)==5&samples(1:count,5)==1;
start=min(tt(run));score=run&tt>=start+(c.circles-c.scoreCircles)*c.period_s-1e-9;
qe=raw(score,2:3)-raw(score,4:5);ve=vq(score,:)-vd(score,:);
rr=diff(ve)./diff(tt(score));
expected=[1000*sqrt(mean(qe.^2)) sqrt(mean(sum(ve.^2,2))) sqrt(mean(sum(rr.^2,2)))];
assert(max(abs(expected-worker.events{1}.metrics))<1e-10,'Online/offline metrics differ.');
% Restarted/stale heartbeats must not keep motion alive.
s=pacOnlineInitialState(c);s.phase=5;s.since=0;s.lastT=.95;
s.lastHeartbeat=0;s.gamma=c.baselineGamma;s.theta=1;s.phi=1;s.vtheta=2*pi;
u=[1;0;0;0;0;1;1;3;0;1;1;1];
actual=[s.gamma(1)*ones(3,1);s.gamma(2)*ones(2,1);s.gamma(3)*ones(3,1);1;1];
c.heartbeatTimeout_s=.1;
[s,y]=pacOnlineRuntimeStep(s,u,actual,c.settings,zeros(7,1),c);
assert(s.phase==9&&s.reason==5&&y(7)==0);
for t=1.05:.05:1.4,u(1)=t;u(6)=t;[s,y]=pacOnlineRuntimeStep(s,u,actual,c.settings,zeros(7,1),c);end
assert(s.phase==10&&~y(6));
% Out-of-bounds proposals, changed tuning, pressure limiting and clock failure.
c.heartbeatTimeout_s=3;
s=pacOnlineInitialState(c);s.phase=1;s.stationarySince=0;
u(1)=1;u(6)=1;rx=[c.session;1;1;1;10;10;10];
[r,~]=pacOnlineRuntimeStep(s,u,actual,c.settings,rx,c);assert(r.reason==10);
s.phase=5;rx=[c.session;2;1;1;c.baselineGamma];bad=c.settings;bad(6)=1;
[r,~]=pacOnlineRuntimeStep(s,u,actual,bad,rx,c);assert(r.reason==6);
u(9)=1;[r,~]=pacOnlineRuntimeStep(s,u,actual,c.settings,rx,c);assert(r.reason==8);
s.lastT=2;[r,y]=pacOnlineRuntimeStep(s,u,actual,c.settings,rx,c);assert(r.reason==1&&~y(6));
% BO infeasible ask/tell remains inside a multiplicative local region.
w=pacOnlineWorkerStep([],zeros(16,1),hc,c);
p=[c.session;2;1;1;100;0;3;3;.003;.02;1;0;40;c.baselineGamma];
w=pacOnlineWorkerStep(w,p,hc,c);assert(w.lastCompleted==1&&w.action==1);
assert(all(w.gamma>=c.baselineGamma/1.5)&&all(w.gamma<=c.baselineGamma*1.5));
again=pacOnlineWorkerStep(w,p,hc,c);assert(height(again.observations)==1);
fprintf('PASS: synthetic continuous runtime, scored batches, preserved gains, worker BO, duplicate results, heartbeat, configuration, bounds, limits and clock stop.\n');
end
