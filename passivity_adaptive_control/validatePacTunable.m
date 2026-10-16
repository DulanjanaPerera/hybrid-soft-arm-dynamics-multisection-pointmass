function validatePacTunable()
%VALIDATEPACTUNABLE Offline force feasibility, bounds and switch integration.
[g,sc,mr,lo,hi,ini,op,h,pv,sv]=pacTunableDefaults();
[~,c]=pacTuningConfig(g,sc,mr,lo,hi,ini,op,h,pv,sv);
% Force that is infeasible at fixed mean is exact with a shifted mean.
req=[25;25];bias=zeros(2,1);
[p,eq,res,lim]=pacAllocateVariableCommon(req,bias,c);
assert(~lim && norm(res)<1e-9 && norm(eq-req)<1e-9);
assert(abs(mean(p)-c.common_bar)>0.01);
% Sweep both force directions, including infeasible requests.
for a=-40:4:40
 for b=-40:4:40
  req=[a;b];[p,eq,res,lim]=pacAllocateVariableCommon(req,bias,c);
  assert(all(p>=max(c.minimum_bar,c.activation_bar)-1e-9) && all(p<=c.maximum_bar+1e-9));
  assert(abs((p(2)+p(3))/2-p(1)-c.sTrim_bar)<=c.maxDeltaS_bar+1e-9);
  assert(abs((p(3)-p(2))/2-c.aTrim_bar)<=c.maxDeltaA_bar+1e-9);
  assert(norm(eq-req-res)<1e-9 && lim==(norm(res)>1e-8));
  % Requested vector is projected onto convex force polygon: variational
  % inequality against feasible pressure-box vertices (within correction caps).
  for k=0:7
   v=max(c.minimum_bar,c.activation_bar)+(c.maximum_bar-max(c.minimum_bar,c.activation_bar))*[bitget(k,1);bitget(k,2);bitget(k,3)];
   if abs((v(2)+v(3))/2-v(1)-c.sTrim_bar)<=c.maxDeltaS_bar && abs((v(3)-v(2))/2-c.aTrim_bar)<=c.maxDeltaA_bar
    f=c.forcePerBar_N*[v(2)-v(1);v(3)-v(1)];assert((req-eq)'*(f-eq)<1e-6);
   end
  end
 end
end
% Common allocator switch at runtime preserves the existing supervisor slew.
clear pacAdaptiveTunableBlock pacAdaptiveCommonPressureBlock
prev=zeros(3,1);shifted=false;ini(1)=700;
for t=0:.05:24
 on=t>=22;
 desired=zeros(6,1);if on,desired(1:2)=.014;end
 [p,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,theta,~,~,~,~,common]= ...
 pacAdaptiveTunableBlock(zeros(4,1),desired,[0;0;-9.81],true,true,t,t,true, ...
 [40;40],[20;20],false,on,g,sc,mr,lo,hi,ini,op,h,pv,sv);
 assert(all(abs(p-prev)<=c.slew_barPerS*.05+1e-9));prev=p;
 assert(abs(common-mean(p))<1e-10);
 assert(theta(1)==700);
 if on && abs(common-c.common_bar)>.01,shifted=true;end
 if ~on
  old=pacAdaptiveCommonPressureBlock(zeros(4,1),zeros(6,1),[0;0;-9.81],true,true,t,t,true,[40;40],[20;20],false);
  assert(norm(old-p)<1e-9);
 end
end
assert(shifted,'Common-pressure switch must affect a demanding force request.');
clear pacAdaptiveTunableBlock pacAdaptiveCommonPressureBlock
fprintf('PASS: variable common allocation and tunable wrapper offline checks.\n');
end
