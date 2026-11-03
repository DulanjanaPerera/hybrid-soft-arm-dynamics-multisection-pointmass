function prediction=greyboxSimulate(w,theta,refined,boundMode)
%GREYBOXSIMULATE Existing recursive mechanics, optional effective input lag.
% theta=[k1 k2 k3 d2 d3 b2 b3 T delay]. No measured future q in propagation.
% Default0 preserves greybox_20261102_v2; mode1 is a NEW corrected baseline.
if nargin<3,refined=false;end
if nargin<4,boundMode=0;end
q0=w.x0(:);start=w.times(1);
time=w.commandTime(:)+theta(9);pressure=w.command;
time=[w.commandTime(1);time];pressure=[pressure(1,:);pressure];
[time,ix]=unique(time,'last');pressure=pressure(ix,:);
T=theta(8);z=pressure;
if T>0
    for i=2:numel(time)
        z(i,:)=pressure(i-1,:)+(z(i-1,:)-pressure(i-1,:))*exp(-(time(i)-time(i-1))/T);
    end
end
lookup=griddedInterpolant(time,(1:numel(time))','previous','nearest');
params=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81], ...
    'K',[theta(1)+theta(3),theta(3);theta(3),theta(2)+theta(3)], ...
    'D',diag(theta(4:5)),'mu',2000,'lKbounds',[-.035;.035;1e6],'tau',[0;0], ...
    'boundMode',boundMode);
rtol=1e-4;atol=1e-7;maxStep=.02;
if refined,rtol=1e-5;atol=1e-8;maxStep=.01;end
opts=odeset('RelTol',rtol,'AbsTol',atol,'MaxStep',maxStep,'Events',@domain);
[tt,x]=ode15s(@rhs,w.times-start,q0,opts);
assert(numel(tt)==numel(w.times)&&all(isfinite(x(:))),'greybox:domain','Failed or domain-limited rollout.');
prediction=x(:,1:2);
    function dx=rhs(t,x)
        physical=t+start;j=round(lookup(physical));p=pressure(j,:);
        if T>0,p=p+(z(j,:)-p)*exp(-max(0,physical-time(j))/T);end
        params.tau=pi*.0065^2*1e5*[-1 1 0;-1 0 1]*p'+theta(6:7)';
        dx=armS_single_dynamics(t,x,params);
    end
    function [value,stop,direction]=domain(~,x)
        value=.045-max(abs(x(1:2)));stop=1;direction=-1;
    end
end
