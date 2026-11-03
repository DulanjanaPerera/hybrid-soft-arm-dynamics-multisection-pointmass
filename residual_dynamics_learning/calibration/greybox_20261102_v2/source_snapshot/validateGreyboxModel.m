function validateGreyboxModel()
%VALIDATEGREYBOXMODEL Independent RHS and delayed-input causal checks.
w=struct('times',(0:.05:.3)','x0',zeros(4,1),'truth',zeros(7,2), ...
    'commandTime',(0:.01:.5)','command',repmat([.3 .5 .2],51,1));
t=[650 650 650 40 40 0 0 0 0];q=greyboxSimulate(w,t,true);
p=struct('N',1,'L',.17411,'r',.013,'mi',.1,'g',[0;0;-9.81], ...
    'K',650*[2 1;1 2],'D',40*eye(2),'mu',2000,'lKbounds',[-.035;.035;1e6], ...
    'tau',pi*.0065^2*1e5*[-1 1 0;-1 0 1]*[.3;.5;.2]);
[~,x]=ode15s(@(time,y)armS_single_dynamics(time,y,p),w.times,w.x0, ...
    odeset('RelTol',1e-5,'AbsTol',1e-8,'MaxStep',.01));
assert(max(abs(q-x(:,1:2)),[],'all')<1e-7);
changed=w;changed.truth(:)=1e9;
assert(isequal(q,greyboxSimulate(changed,t,true)),'Future measured output leaked into simulator.');
w.command(:)=0;zero=greyboxSimulate(w,t,true);
w.command(w.commandTime>=.1,2)=1;t(8:9)=[.2 .1];delayed=greyboxSimulate(w,t,true);
assert(max(abs(delayed(w.times<.19,:)-zero(w.times<.19,:)),[],'all')<1e-6, ...
    'Response preceded delayed pressure command.');
fprintf('PASS: zero-lag original nonlinear RHS, no future-state leakage, delayed-input causality.\n');
end
