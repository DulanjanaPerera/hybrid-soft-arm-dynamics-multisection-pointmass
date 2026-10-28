function validatePacOnlineRejection(folder)
%VALIDATEPACONLINEREJECTION Recorded completed-motion rejection, no hardware.
w=load(fullfile(folder,'worker_log.mat'));d=load(fullfile(folder,'pacOnlineSession.mat'));
p=w.packets;k=find(p(:,4)==1&p(:,3)==1,1);assert(~isempty(k));
tx=p(k,1:16)';s=pacOnlineWorkerStep([],zeros(16,1),d.c,d.rt);
s=pacOnlineWorkerStep(s,tx,d.c,d.rt);
assert(s.action==4&&strcmp(s.reason,'insufficient_coverage'));
assert(s.requestId==2&&height(s.observations)==0&&numel(s.attempts)==1);
again=pacOnlineWorkerStep(s,tx,d.c,d.rt);assert(isequaln(s,again));
% Runtime receives the error addressed to request2 while still settling.
r=pacOnlineInitialState(d.rt);r.phase=1;r.id=1;r.requestId=2;r.completedId=1;
r.since=51.6;r.lastT=51.6;r.stationarySince=-1;
u=[51.65;0;0;.001;0;51.65;1;3;0;1;1;1];
g=d.rt.baselineGamma;actual=[g(1)*ones(3,1);g(2)*ones(2,1);g(3)*ones(3,1);1;1];
rx=[d.rt.session;1;s.requestId;s.action;s.gamma];
[r,~,~]=pacOnlineRuntimeStep(r,u,actual,d.rt.settings,rx,d.rt);
assert(r.phase==9&&r.reason==24,'Terminal reply should not wait30s or for stationarity.');
fprintf('PASS: recorded coverage rejection retains attempt, replies to request2, remains terminal and immediately starts controlled shutdown.\n');
end
