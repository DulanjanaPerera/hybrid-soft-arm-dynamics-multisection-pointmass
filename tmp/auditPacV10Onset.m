root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'passivity_adaptive_control'));
folder=fullfile(root,'passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_gainTesting_v10_analysis');
S=load(fullfile(fileparts(folder),'1rad_10s_gainTesting_v10.mat'));ds=S.out.logsout;
[t,th]=read(ds,'parameterEstimate');D=struct;
names={'parameterLower','modelValid','status','limited','adaptationActive','allocationResidual_N','pressureCommand_bar','pressureOptions','measurementTime_s','sensorValid','configurationValid','baseOreintationValid','enablePressure','kill_switch','adaptationGamma','tauRequested_N','tauBlended_N','forceBias_N','filtered_d_length','d_length','filtered_dd_length','meas_length','sampled des_length','theta','phi','des_theta','des_phi','orientationValid','thetaObservable','dampingForce_N','hysteresisForce_N','staticForce_N','feedbackForce_N','inertialForce_N','adaptiveForce_N','historyState'};
for k=1:numel(names),[tt,x]=read(ds,names{k});D.(matlab.lang.makeValidName(names{k}))=interp1(tt,x,t,'previous','extrap');end
fid=fopen(fullfile(folder,'onset_audit.txt'),'w');fclose(fid);diary(fullfile(folder,'onset_audit.txt'));
ev=(th(:,1)+th(:,2)+2*th(:,3)-sqrt((th(:,1)-th(:,2)).^2+4*th(:,3).^2))/2;
qerror=D.meas_length(:,2:3)-D.sampledDes_length;vError=D.filtered_d_length-D.d_length;
thetaErr=atan2(sin(D.des_theta-D.theta),cos(D.des_theta-D.theta));
thetaRate=[NaN;atan2(sin(diff(D.theta)),cos(diff(D.theta)))./diff(t)];
good=D.modelValid&D.orientationValid&D.thetaObservable&D.status==3;
pair=good&[false;good(1:end-1)];thetaRate(~pair)=NaN;
fprintf('REFERENCE median angular rate%g period%g\n',median(diff(D.des_theta(t>=20))./.05),2*pi/median(diff(D.des_theta(t>=20))./.05));
for stamp=[20 30 40 45 50 52 54 55 55.2 56 60 65 70 75 78.7 85 93.5 100]
 i=find(t>=stamp,1);fprintf('TIME%.2f spring%s d%s h%s minEig%.3f status%d limit%d phi%.4f thetaErr%.4f\n',t(i),mat2str(th(i,1:3),5),mat2str(th(i,4:5),5),mat2str(th(i,6:8),5),ev(i),D.status(i),D.limited(i),D.phi(i),thetaErr(i));
end
for j=1:3,ii=find(t>=20&th(:,j)<=D.parameterLower(:,j)+1e-9);fprintf('BOUND k%d count%d first%s\n',j,numel(ii),mat2str(t(ii(1:min(1,end)))));end
for win=[36 44 52 60 68 76 84 94;44 52 60 68 76 78.7 92 102]
 ix=t>=win(1)&t<win(2)&good;allix=t>=win(1)&t<win(2);ir=t>=win(1)&t<win(2)&pair;
 fprintf('WIN%.2f-%.2f validControlPct%.2f thetaMAE%.5f phiMAE%.5f phiRange%s lengthRMSEmm%s velocityRMSEmm/s%s thetaRateSD%.5f thetaSlowPct%.2f limitPct%.2f dMean%s minEigRange%s\n',win(1),win(2),100*mean(good(allix)),mean(abs(thetaErr(ix))),mean(abs(D.des_phi(ix)-D.phi(ix))),mat2str([min(D.phi(ix)) max(D.phi(ix))],5),mat2str(1000*sqrt(mean(qerror(ix,:).^2)),5),mat2str(1000*sqrt(mean(vError(ix,:).^2)),5),std(thetaRate(ir)),100*mean(thetaRate(ir)<.2*pi/2),100*mean(D.limited(allix)),mat2str(mean(th(ix,4:5)),5),mat2str([min(ev(ix)) max(ev(ix))],5));
 for name={'dampingForce_N','hysteresisForce_N','staticForce_N','feedbackForce_N','inertialForce_N'}
  x=D.(name{1});fprintf('FORCE%s RMS%s peak%s\n',name{1},mat2str(sqrt(mean(x(ix,:).^2)),5),mat2str(max(abs(x(ix,:))),5));
 end
end
fault=find(D.status==5,1);fprintf('FIRSTFAULT%.3f\n',t(fault));ii=max(1,fault-8):min(numel(t),fault+4);
fprintf('FAULT rows time status sensorValid1/2 configValid baseValid age enable kill phi\n');disp([t(ii) D.status(ii) D.sensorValid(ii,:) D.configurationValid(ii) D.baseOreintationValid(ii) t(ii)-D.measurementTime_s(ii) D.enablePressure(ii) D.kill_switch(ii) D.phi(ii)]);
% Replay allocation before fault to distinguish force infeasibility/slew.
p=D.pressureOptions(1,:);cfg=pacPressureConfig();fields={'common_bar','sTrim_bar','aTrim_bar','minimum_bar','maximum_bar','maxDeltaS_bar','maxDeltaA_bar','slew_barPerS','activation_bar','forcePerBar_N'};for j=1:10,cfg.(fields{j})=p(j);end
al=false(size(t));sl=al;
for i=find(good&t>=20)'
 [pp,~,~,ll]=pacAllocateVariableCommon(D.tauBlended_N(i,:)',D.forceBias_N(i,:)',cfg);al(i)=ll;sl(i)=any(abs(pp'-D.pressureCommand_bar(i,:))>1e-9);
end
for win=[20 44 52 60 78.7;44 52 60 78.7 t(fault)]
 ix=t>=win(1)&t<win(2)&good;fprintf('ALLOC%.2f-%.2f eligible%d allocation%d slew%d loggedLimited%d\n',win(1),win(2),nnz(ix),nnz(ix&al),nnz(ix&sl),nnz(ix&D.limited));
end
fprintf('CONTIGUOUS event50-58 invalid%s limited%s\n',mat2str(t(t>=50&t<=58&~D.modelValid)'),mat2str(t(t>=50&t<=58&D.limited)'));
fprintf('freeze/enable changes are in main audit; no pressure stop rule at zero error.\n');
fig=figure('Visible','off','Color','w','Position',[60 60 1400 1200]);tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
nexttile;plot(t,D.phi,'LineWidth',1.2);hold on;plot(t,D.des_phi,'k--');ylabel('\phi (rad)');
nexttile;plot(t,1000*vError,'LineWidth',1);ylabel('Velocity error (mm/s)');
nexttile;plot(t,th(:,4:5),'LineWidth',1.3);ylabel('d estimates (N s/m)');legend('d_2','d_3','Location','northwest');
nexttile;plot(t,ev,'LineWidth',1.2);yline(0,'k:');ylabel('min eig(K) (N/m)');
nexttile;plot(t,D.pressureCommand_bar,'LineWidth',1.2);ylabel('Pressure (bar)');xlabel('Time (s)');
for a=findall(fig,'Type','axes')',set(a,'FontName','Arial','FontSize',15,'XGrid','on','YGrid','on');xlim(a,[20 t(end)]);xline(a,55,'k--','HandleVisibility','off');xline(a,78.7,'k:','HandleVisibility','off');xline(a,t(fault),'r--','HandleVisibility','off');end
exportgraphics(fig,fullfile(folder,'onset_timeline.png'),'Resolution',200);savefig(fig,fullfile(folder,'onset_timeline.fig'));close(fig);diary off
function [t,x]=read(ds,name)
ix=find(strcmp(ds.getElementNames,name),1);v=ds{ix}.Values;t=double(v.Time(:));raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
