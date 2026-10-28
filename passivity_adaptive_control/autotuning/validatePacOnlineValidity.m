function validatePacOnlineValidity(folder)
%VALIDATEPACONLINEVALIDITY Recorded-prefix regression, no hardware/model.
% Replay stops shortly after the old abort; not counterfactual tracking proof.
w=load(fullfile(folder,'worker_log.mat'));d=load(fullfile(folder,'pacOnlineSession.mat'));
c=d.rt;s=pacOnlineInitialState(c);p=w.packets;count=0;
for k=1:size(p,1)
    if p(k,5)>34.2+1e-9,break;end
    u=p(k,17:28)';g=p(k,31:40)';settings=p(k,41:121)';
    rx=[c.session;k;1;1;c.baselineGamma];
    [s,~,~]=pacOnlineRuntimeStep(s,u,g,settings,rx,c);count=count+1;
    assert(s.reason==0,'Recorded-prefix regression aborted at %.3f, reason %d.',u(1),s.reason);
end
assert(s.phase==5&&count>600);
% A sustained raw dropout still stops; the corrected exact grace boundary is
% accepted, but one more sample is rejected. No grace interval increase.
u=p(k-1,17:28)';g=p(k-1,31:40)';settings=p(k-1,41:121)';
start=34.25;u(7)=0;
for j=0:2
    u(1)=start+j*.05;u(6)=u(1);rx=[c.session;1000+j;1;1;c.baselineGamma];
    [s,~,~]=pacOnlineRuntimeStep(s,u,g,settings,rx,c);
    assert(s.reason==0,'Premature failure at exact grace boundary.');
end
u(1)=34.4;u(6)=u(1);rx(2)=1003;
[s,~,~]=pacOnlineRuntimeStep(s,u,g,settings,rx,c);assert(s.reason==4&&s.phase==9);
% A single invalid hold frame is not mislabeled as physical movement.
s=pacOnlineInitialState(c);s.phase=2;s.since=0;s.lastT=.95;s.gamma=c.baselineGamma;
u=[1;0;0;0;0;1;0;4;0;1;0;1];rx=[c.session;1;1;1;c.baselineGamma];
[s,~,~]=pacOnlineRuntimeStep(s,u,g,settings,rx,c);assert(s.reason==0&&s.phase==2);
u(1)=1.05;u(6)=1.05;u(7)=1;u(8)=3;u(11)=1;u(4)=.01;rx(2)=2;
[s,~,~]=pacOnlineRuntimeStep(s,u,g,settings,rx,c);assert(s.reason==12);
fprintf('PASS: recorded-prefix isolated-dropout regression, sustained-dropout stop, grace boundary and hold-frame/movement distinction.\n');
end
