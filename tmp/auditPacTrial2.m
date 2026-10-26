root=fileparts(fileparts(mfilename('fullpath')));
S=load(fullfile(root,'passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_tr2.mat'));ds=S.out.logsout;
[t,th]=read(ds,'parameterEstimate');D=struct;names=ds.getElementNames;
for k=1:numel(names)
 name=names{k};if isempty(name),continue;end
 try,[tt,x]=read(ds,name);if numel(tt)>1,D.(matlab.lang.makeValidName(name))=interp1(tt,x,t,'previous','extrap');else,D.(matlab.lang.makeValidName(name))=repmat(x,numel(t),1);end;catch,end
end
diary(fullfile(root,'tmp/trial2_audit.txt'));fprintf('N%d time %g %g\n',numel(t),t(1),t(end));
for name={'lambda_gains','damping_gains','enableAdaptation','enableCommonPressureAdaptation','adaptationGamma','parameterMaxRate','parameterScale'}
 if isfield(D,name{1}),x=D.(name{1});fprintf('%s min %s max %s\n',name{1},mat2str(min(x),6),mat2str(max(x),6));else,fprintf('MISSING %s\n',name{1});end
end
fprintf('PARAM initial %s final %s min %s max %s\n',mat2str(th(1,:),6),mat2str(th(end,:),6),mat2str(min(th),6),mat2str(max(th),6));
z=all(D.dampingEstimate_NsPerM==0,2);starts=find(z & ~[false;z(1:end-1)]);ends=find(z & ~[z(2:end);false]);
fprintf('ZERO samples%d events%d circles%d invalidOverlap%d hold%d fault%d\n',nnz(z),numel(starts),nnz(z & t>=20 & t<320),nnz(z & ~D.modelValid),nnz(z & D.status==4),nnz(D.status==5));
fprintf('EVENTS start end count\n');disp([t(starts) t(ends) ends-starts+1]);
for name={'configurationValid','orientationValid','sensorValid','baseOreintationValid','valid','thetaObservable'}
 if isfield(D,name{1}),x=D.(name{1});fprintf('INVALID %s overall %s atZeros %s\n',name{1},mat2str(sum(~x)),mat2str(sum(~x(z,:))));end
end
age=t-D.measurementTime_s;fprintf('AGE zero minmax %s good %s futureZeros%d\n',mat2str([min(age(z)) max(age(z))],6),mat2str([min(age(~z)) max(age(~z))],6),nnz(age(z)<-1e-9));
for name={'inertialForce_N','coriolisForce_N','dampingForce_N','staticForce_N','feedbackForce_N','hysteresisForce_N','adaptiveBias_N','adaptiveForce_N','tauRequested_N','tauEquivalent_N'}
 x=D.(name{1});fprintf('FORCE %s allZeroAtDrops%d of%d finite%d zeroPeak%s\n',name{1},nnz(all(x(z,:)==0,2)),nnz(z),all(isfinite(x(:))),mat2str(max(abs(x(z,:))),5));
end
fprintf('PRESSURE atDrops min%s max%s peakChange%s\n',mat2str(min(D.pressureCommand_bar(z,:)),5),mat2str(max(D.pressureCommand_bar(z,:)),5),mat2str(max(abs(D.pressureCommand_bar(starts,:)-D.pressureCommand_bar(starts-1,:))),5));
fprintf('RATE peaks%s\n',mat2str(max(abs(diff(th)./diff(t))),6));
e=atan2(sin(D.des_theta-D.theta),cos(D.des_theta-D.theta));ep=D.des_phi-D.phi;el=D.sampledDes_length-D.meas_length(:,2:3);
fprintf('LAP thetaMAE phiMAE lengthRMSE hMean zRMS hForceRMS kMean adaptPct limitPct\n');
for lap=1:30
 ix=t>=20+(lap-1)*10 & t<20+lap*10;good=ix & D.modelValid & D.orientationValid & D.thetaObservable;
 fprintf('%d %.5f %.5f %s %s %s %s %s %.2f %.2f\n',lap,mean(abs(e(good))),mean(abs(ep(good))),mat2str(1000*sqrt(mean(el(good,:).^2)),5),mat2str(mean(th(ix,6:8)),5),mat2str(sqrt(mean(D.historyState(ix,:).^2)),5),mat2str(sqrt(mean(D.hysteresisForce_N(good,:).^2)),5),mat2str(mean(th(ix,1:3)),5),100*mean(D.adaptationActive(ix)),100*mean(D.limited(ix)));
end
% Gradient directions and feature overlap on valid control samples.
q=D.meas_length(:,2:3);zh=D.historyState;ss=D.s_mPerS;
gradH=[-zh(:,1).*sum(ss,2),zh(:,2).*ss(:,1),zh(:,3).*ss(:,2)];
gradK=[q(:,1).*ss(:,1),q(:,2).*ss(:,2),sum(q,2).*sum(ss,2)];
for range=[20 120 220;120 220 320]
 ix=t>=range(1) & t<range(2) & D.adaptationActive;
 fprintf('WINDOW%d-%d meanMinusYTs H%s K%s\n',range(1),range(2),mat2str(-mean(gradH(ix,:)),6),mat2str(-mean(gradK(ix,:)),6));
end
% Analyze compensation force components separately, using current nominal K/D.
v=D.d_length-D.lambda_gains.*(q-D.sampledDes_length);
fk=[(th(:,1)+th(:,3)-1300).*q(:,1)+(th(:,3)-650).*q(:,2),(th(:,3)-650).*q(:,1)+(th(:,2)+th(:,3)-1300).*q(:,2)];
fd=(th(:,4:5)-40).*v;
for lap=[2 10 20 30]
 ix=t>=20+(lap-1)*10 & t<20+lap*10 & D.modelValid;
 fprintf('FORCECOMP lap%d stiffnessRMS%s dampingRMS%s historyRMS%s totalAdaptiveRMS%s\n',lap,mat2str(sqrt(mean(fk(ix,:).^2)),5),mat2str(sqrt(mean(fd(ix,:).^2)),5),mat2str(sqrt(mean(D.hysteresisForce_N(ix,:).^2)),5),mat2str(sqrt(mean(D.adaptiveForce_N(ix,:).^2)),5));
end
fprintf('simout27 first/last\n');disp([reshape(S.out.simout.Data(:,:,1),[],1),reshape(S.out.simout.Data(:,:,end),[],1)]);
fprintf('DROP frame evidence\n');
for j=find(strcmp(names,'frame'))'
 vv=ds{j}.Values;xx=reshape(double(vv.Data),[],numel(vv.Time))';prev=xx(starts-1,:);cur=xx(starts,:);
 fprintf('frame%d nonfinite%s repeated%s mismatch%d\n',j,mat2str(sum(~isfinite(cur))),mat2str(sum(cur==prev)),nnz(cur(:,1)~=cur(:,2)));
 disp([t(starts(1:6)) prev(1:6,:) cur(1:6,:)]);
end
fprintf('DROP raw1-8 first6 rows\n');raw=zeros(numel(t),8);for j=1:8,raw(:,j)=D.(sprintf('raw%d',j));end;disp([t(starts(1:6)) raw(starts(1:6),:)]);
ix=t>=20 & t<320 & D.modelValid;
fprintf('PRESSURE range%s %s common%s residualPeak%g\n',mat2str(min(D.pressureCommand_bar(ix,:)),5),mat2str(max(D.pressureCommand_bar(ix,:)),5),mat2str([min(D.commonPressure_bar(ix)) max(D.commonPressure_bar(ix))],5),max(vecnorm(D.allocationResidual_N(ix,:),2,2)));
fprintf('PARAM change after drop max%s HISTORY change max%s\n',mat2str(max(abs(th(starts+1,:)-th(starts,:))),5),mat2str(max(abs(D.historyState(starts+1,:)-D.historyState(starts,:))),5));
% History feature correlation with muscle velocity, to test lag/damping overlap.
vel=D.filtered_d_length;if size(vel,2)==2,vel=[-sum(vel,2) vel];end
if size(vel,2)==3
 fprintf('HISTORY velocity correlation%s\n',mat2str(diag(corr(zh(ix,:),vel(ix,:)))',5));
end
fig=figure('Visible','off','Position',[100 100 1100 850]);tiledlayout(4,1);
nexttile;plot(t,th(:,1:3));ylabel('k (N/m)');legend('k1','k2','k3');grid on
nexttile;plot(t,th(:,6:8));ylabel('h weights (N)');legend('h1','h2','h3');grid on
nexttile;hf=D.hysteresisForce_N;hf(~D.modelValid,:)=NaN;plot(t,hf);ylabel('History force (N)');grid on
nexttile;plot(t,D.pressureCommand_bar);hold on;scatter(t(starts),repmat(.7,numel(starts),1),12,'r','filled');ylabel('Pressure command (bar)');xlabel('Time (s); red marks: invalid sensor frames');grid on
exportgraphics(fig,fullfile(root,'tmp/trial2_diagnostics.png'),'Resolution',150);close(fig);
diary off
function [t,x]=read(ds,name)
names=ds.getElementNames;v=ds{find(strcmp(names,name),1)}.Values;t=v.Time(:);raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
