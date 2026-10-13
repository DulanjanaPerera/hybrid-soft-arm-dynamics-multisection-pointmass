function report=validatePacAdaptive()
%VALIDATEPACADAPTIVE Offline checks only; no model load or hardware.
folder=fileparts(mfilename('fullpath'));addpath(folder,fileparts(folder));
cfg=pacAdaptiveConfig();st=pacAdaptiveInitialState();
m=[.001;-.002;.0001;-.0002];des=[.002;-.001;0;0;0;0];g=[0;0;-9.81];
[u,s,Y,fi,fc,fd,fs,ff,fh,fb,K,D,ok]= ...
 pacAdaptiveForce(m,des,g,[200;200],[20;20],st,true);
assert(ok && norm(u-fi-fc-fd-fs-ff-fh-fb)<1e-10);
assert(min(eig(K))>0 && min(eig(D))>0 && norm(K-K.')==0);
q=m(1:2);v=-200*(q-des(1:2));
assert(norm(Y*st.theta-K*q-D*v)<1e-10,'Regressor mismatch.');
[~,~,~,~,~,~,~,~,~,~,~,~,ok]=pacAdaptiveForce(m,des,g,[0;200],[20;20],st,true);
assert(~ok,'Zero gain accepted.');
st0=st;[st,active]=pacAdaptiveUpdate(st,Y,s,.05,false);
assert(isequal(st,st0)&&~active,'Disabled adaptation moved estimates.');
[st,active]=pacAdaptiveUpdate(st,Y,s,0,true);assert(~active && isequal(st,st0));
[st,active]=pacAdaptiveUpdate(st,Y,s,.2,true);assert(~active && isequal(st,st0));
[st,active]=pacAdaptiveUpdate(st,Y,s,.05,true);
assert(active && s.'*Y*(st.theta-st0.theta)<=1e-12,'Wrong descent sign.');
for i=1:2000
 [st,~]=pacAdaptiveUpdate(st,ones(2,10),[1;1],.05,true);
end
assert(all(st.theta>=cfg.lower & st.theta<=cfg.upper));
[~,~,~,~,~,~,~,~,~,~,K,D]=pacAdaptiveForce(m,des,g,[200;200],[20;20],st,true);
assert(min(eig(K))>0&&min(eig(D))>0);
z=zeros(3,1);
for i=1:200,z=pacAdaptiveHistory(z,[.01;-.005],.05);end
assert(all(isfinite(z))&&max(abs(z))<=cfg.historyBound);
zBefore=z;for i=1:200,z=pacAdaptiveHistory(z,[-.01;.005],.05);end
assert(norm(z-zBefore)>.1,'History ignored reversal.');
assert(isequal(z,pacAdaptiveHistory(z,[NaN;0],.05)));
% Projection/rate limitation preserves descent sign but is not exact
% continuous-time Lyapunov cancellation. Test ideal unprojected identity.
st=pacAdaptiveInitialState();[p,~,~]=pacNominalParameters(g);
truth=st.theta+[200;100;50;5;10;1;2;1;.2;-.1];
[M,C,~,dM]=armS_single_core(m(1:2),m(3:4),p);
Delta=st.theta-truth;Dtrue=diag(truth(4:5));
sdot=M\(Y*Delta-(C+Dtrue+diag([20;20]))*s);
tdot=-(cfg.scale.^2).*cfg.gamma.*(Y.'*s);
mdot=dM(:,:,1)*m(3)+dM(:,:,2)*m(4);
vdot=s.'*M*sdot+.5*s.'*mdot*s+Delta.'*(tdot./(cfg.scale.^2.*cfg.gamma));
assert(abs(vdot+s.'*(Dtrue+diag([20;20]))*s)<1e-8,'Ideal energy sign failed.');
% Combined wrapper startup, active updates, fresh-only update, invalid hold.
clear pacAdaptiveCommonPressureBlock
for t=0:.05:20.1
 out=cell(1,26);[out{:}]=pacAdaptiveCommonPressureBlock(zeros(4,1),zeros(6,1), ...
  g,true,true,t,t,true,[200;200],[20;20],true);
end
assert(out{13}==3 && out{15});
des=[.00005;.00005;0;0;0;0];
[out{:}]=pacAdaptiveCommonPressureBlock(zeros(4,1),des,g,true,true,20.15,20.15,true,[200;200],[20;20],true);
assert(out{24},'Engaged adaptation failed.');
[out{:}]=pacAdaptiveCommonPressureBlock(zeros(4,1),des,g,true,true,20.15,20.15,true,[200;200],[20;20],true);
assert(~out{24},'Duplicate timestamp adapted.');
[out{:}]=pacAdaptiveCommonPressureBlock(zeros(4,1),des,g,false,true,20.15,20.2,true,[200;200],[20;20],true);
assert(out{13}==4 && ~out{24});
[out{:}]=pacAdaptiveCommonPressureBlock(zeros(4,1),des,g,false,true,20.15,20.4,true,[200;200],[20;20],true);
assert(out{13}==5 && ~out{24});
report=struct('regressorAndForceBalance',true,'projectionAndSPD',true, ...
 'historyReversal',true,'idealEnergyIdentity',true,'wrapperFreshnessAndFaults',true, ...
 'scope','Offline algebra/supervisor checks; sampled physical stability not established');
disp(report);fprintf('PASS: adaptive PAC offline checks.\n');
end
