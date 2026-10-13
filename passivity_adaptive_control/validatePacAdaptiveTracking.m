function report=validatePacAdaptiveTracking()
%VALIDATEPACADAPTIVETRACKING Ideal-force continuous-time simulation, no valves.
% Tests a chosen representable mismatch, NOT the recorded physical plant.
cfg=pacAdaptiveConfig();st=pacAdaptiveInitialState();
truth=[1000;900;1100;55;65;1;.7;1.2;.3;-.2];
y0=[zeros(4,1);st.theta;zeros(3,1)];
options=odeset('RelTol',1e-5,'AbsTol',1e-8,'MaxStep',.05);
[t,yFixed]=ode15s(@(t,y)rhs(t,y,false),0:.05:30,y0,options);
[~,yAdapt]=ode15s(@(t,y)rhs(t,y,true),0:.05:30,y0,options);
qd=zeros(numel(t),2);
for i=1:numel(t),d=ref(t(i));qd(i,:)=d(1:2)';end
ix=t>20;fixed=sqrt(mean((yFixed(ix,1:2)-qd(ix,:)).^2,1));
adapt=sqrt(mean((yAdapt(ix,1:2)-qd(ix,:)).^2,1));
assert(all(isfinite(yAdapt(:))));
assert(all(yAdapt(:,5:14)>=cfg.lower'-1e-5,'all') && ...
 all(yAdapt(:,5:14)<=cfg.upper'+1e-5,'all'));
report=struct('fixedRMSE_m',fixed,'adaptiveRMSE_m',adapt, ...
 'finalEstimates',yAdapt(end,5:14),'scope', ...
 'Ideal continuous force, lambda=200/Ks=20; no sampled sensor/pressure plant');
disp(report);
 function dy=rhs(time,y,adaptEnabled)
  state.theta=y(5:14);state.z=y(15:17);state.lastSample=time;
  d=ref(time);[p,~,~]=pacNominalParameters([0;0;-9.81]);
  [tau,s,Y]=pacAdaptiveForce(y(1:4),d,p.g,[200;200],[20;20],state,true);
  q=y(1:2);dq=y(3:4);[M,C,G]=armS_single_core(q,dq,p);
  K=[truth(1)+truth(3),truth(3);truth(3),truth(2)+truth(3)];
  D=diag(truth(4:5));T=[-1 -1;1 0;0 1];
  H=T.'*diag(state.z)*truth(6:8)+truth(9:10);
  ddq=M\(tau-C*dq-G-K*q-D*dq-H);
  [next,~]=pacAdaptiveUpdate(state,Y,s,.001,adaptEnabled);
  dz=(pacAdaptiveHistory(state.z,dq,.001)-state.z)/.001;
  dy=[dq;ddq;(next.theta-state.theta)/.001;dz];
 end
end
function d=ref(t)
w=.4;q=.002*[sin(w*t);cos(w*t)];dq=.002*w*[cos(w*t);-sin(w*t)];
ddq=-w*w*q;d=[q;dq;ddq];
end
