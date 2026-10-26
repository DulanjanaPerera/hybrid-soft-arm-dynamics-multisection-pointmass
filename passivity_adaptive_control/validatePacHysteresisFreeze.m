function validatePacHysteresisFreeze()
%VALIDATEPACHYSTERESISFREEZE Offline synthetic wrapper inputs; no models/I/O.
[g,sc,mr,lo,hi,ini,op,h,pv,sv]=pacTunableDefaults();
ini(6:8)=[1;2;3];
% Engage supervisor before supplying small nonzero position/velocity errors.
sv(1)=0;sv(2)=0;sv(3)=.05;
for mode=1:3
 clear pacAdaptiveTunableBlock
 frozenStart=[];zStart=[];thetaStart=[];lastTheta=[];
 changedOther=false;changedHistory=false;resumed=false;activeFrozen=false;
 for k=0:200
  t=k*.05;measured=zeros(4,1);desired=zeros(6,1);
  if t>=2
   measured=[.001;-.0007;.0005;-.0003];
  end
  freeze=(mode~=1 && t>=3 && t<7);enable=(mode~=3);
  [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,fh,~,theta,z,active]= ...
   pacAdaptiveTunableBlock(measured,desired,[0;0;-9.81],true,true,t,t,true, ...
   [18;16.5],[23;10],enable,true,g,sc,mr,lo,hi,ini,op,h,pv,sv,freeze);
  if freeze
   if isempty(frozenStart),frozenStart=theta(6:8);zStart=z;thetaStart=theta;end
   assert(isequal(theta(6:8),frozenStart),'Frozen history weights changed.');
   other=[1:5 9:10];changedOther=changedOther || any(theta(other)~=thetaStart(other));
   changedHistory=changedHistory || any(z~=zStart);
   activeFrozen=activeFrozen || active;
   assert(norm(fh-[-theta(6)*z(1)+theta(7)*z(2); ...
    -theta(6)*z(1)+theta(8)*z(3)])<1e-12,'History force must remain active.');
  end
  if mode==2 && t>=7.1,resumed=resumed || any(theta(6:8)~=frozenStart);end
  if mode==3,assert(isequal(theta,ini) && ~active);end
  lastTheta=theta;
 end
 if mode==1,assert(any(lastTheta(6:8)~=ini(6:8)));end
 if mode==2,assert(changedOther && changedHistory && activeFrozen && resumed);end
 if mode==3,assert(changedHistory);end
end
clear pacAdaptiveTunableBlock
fprintf('PASS: hysteresis-only freeze, evolving history, other adaptation, resume and global disable.\n');
end
