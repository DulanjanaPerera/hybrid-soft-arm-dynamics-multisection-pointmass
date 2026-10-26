root=fileparts(fileparts(mfilename('fullpath')));
report=fullfile(root,'tmp/trial3_freeze_audit.txt');
fid=fopen(report,'w');assert(fid>=0);fclose(fid);diary(report);
R=cell(2,1);
for run=2:3
 S=load(fullfile(root,sprintf('passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_tr%d.mat',run)));ds=S.out.logsout;
 [t,th]=read(ds,'parameterEstimate');D=struct;
 names={'modelValid','status','limited','adaptationActive','historyState','hysteresisForce_N','adaptiveForce_N','adaptiveBias_N','commonPressure_bar','pressureCommand_bar','allocationResidual_N','lambda_gains','damping_gains','theta','phi','des_theta','des_phi','orientationValid','thetaObservable','sampled des_length','meas_length','dampingEstimate_NsPerM','sensorValid','measurementTime_s'};
 for k=1:numel(names)
  [tt,x]=read(ds,names{k});D.(matlab.lang.makeValidName(names{k}))=interp1(tt,x,t,'previous','extrap');
 end
 fprintf('\nTRIAL%d samples%d time%g-%g\n',run,numel(t),t(1),t(end));
 fprintf('lambda %s %s Ks %s %s\n',mat2str(min(D.lambda_gains)),mat2str(max(D.lambda_gains)),mat2str(min(D.damping_gains)),mat2str(max(D.damping_gains)));
 if run==3
  [tt,x]=read(ds,'freezeHysteresisWeights');ix=[1;find(diff(x)~=0)+1];fprintf('FREEZE transitions [time value]\n');disp([tt(ix) x(ix)]);
  fr=interp1(tt,x,t,'previous','extrap')>0;first=find(fr,1);fprintf('FIRST frozen sample %.8f h%s postFreezeRange%s\n',t(first),mat2str(th(first,6:8),8),mat2str(max(th(fr,6:8))-min(th(fr,6:8)),8));
  fprintf('POSTFREEZE historyRange%s otherParamRange%s activePct%.3f\n',mat2str(max(D.historyState(fr,:))-min(D.historyState(fr,:)),6),mat2str(max(th(fr,[1:5 9:10]))-min(th(fr,[1:5 9:10])),6),100*mean(D.adaptationActive(fr)));
 end
 z=all(D.dampingEstimate_NsPerM==0,2);starts=find(z & ~[false;z(1:end-1)]);
 fprintf('DROP samples%d events%d circles%d invalid%d hold%d tipInvalid%d fault%d\n',nnz(z),numel(starts),nnz(z&t>=20&t<320),nnz(z&~D.modelValid),nnz(z&D.status==4),nnz(z&~D.sensorValid(:,2)),nnz(D.status==5));
 fprintf('DROP times%s\n',mat2str(t(starts)',6));
 e=atan2(sin(D.des_theta-D.theta),cos(D.des_theta-D.theta));ep=D.des_phi-D.phi;el=D.sampledDes_length-D.meas_length(:,2:3);
 metrics=zeros(30,8);
 fprintf('lap thetaMAE phiMAE lengthRMSE hMean historyForceRMS springMean adaptPct limitedPct\n');
 for lap=1:30
  ix=t>=20+(lap-1)*10 & t<20+lap*10;good=ix & D.modelValid & D.orientationValid & D.thetaObservable;
  metrics(lap,:)=[mean(abs(e(good))),mean(abs(ep(good))),1000*sqrt(mean(el(good,:).^2)),sqrt(mean(D.hysteresisForce_N(good,:).^2)),sqrt(mean(D.adaptiveForce_N(good,:).^2))];
  fprintf('%d %.5f %.5f %s %s %s %s %.2f %.2f\n',lap,metrics(lap,1),metrics(lap,2),mat2str(metrics(lap,3:4),5),mat2str(mean(th(ix,6:8)),6),mat2str(metrics(lap,5:6),5),mat2str(mean(th(ix,1:3)),6),100*mean(D.adaptationActive(ix)),100*mean(D.limited(ix)));
 end
 for win=[100 120 220;120 220 320]
  ix=t>=win(1)&t<win(2)&D.modelValid&D.orientationValid&D.thetaObservable;
  fprintf('WIN%d-%d theta%.6f phi%.6f lengths%s historyRMS%s adaptiveRMS%s\n',win(1),win(2),mean(abs(e(ix))),mean(abs(ep(ix))),mat2str(1000*sqrt(mean(el(ix,:).^2)),6),mat2str(sqrt(mean(D.hysteresisForce_N(ix,:).^2)),6),mat2str(sqrt(mean(D.adaptiveForce_N(ix,:).^2)),6));
 end
 fprintf('FINAL param%s\n',mat2str(th(end,:),7));
 ix=t>=20&t<320;fprintf('LIMIT%d common%s pressures%s %s residualPeak%g\n',nnz(D.limited(ix)),mat2str([min(D.commonPressure_bar(ix)) max(D.commonPressure_bar(ix))],6),mat2str(min(D.pressureCommand_bar(ix,:)),6),mat2str(max(D.pressureCommand_bar(ix,:)),6),max(vecnorm(D.allocationResidual_N(ix&D.modelValid,:),2,2)));
 R{run-1}=struct('t',t,'th',th,'D',D,'metrics',metrics);
end
fprintf('COMPARISON trial3 minus trial2 finalLap%s\n',mat2str(R{2}.metrics(30,1:4)-R{1}.metrics(30,1:4),6));
fig=figure('Visible','off','Position',[100 100 1100 850]);tiledlayout(4,1);
nexttile;plot(R{2}.t,R{2}.th(:,6:8));xline(t(first),'k--');ylim([0 1.7]);ylabel('Trial3 h weights (N)');legend('h1','h2','h3');grid on
nexttile;plot(R{2}.t,R{2}.th(:,1:3));xline(t(first),'k--');ylabel('Trial3 k (N/m)');legend('k1','k2','k3');grid on
nexttile;plot(1:30,R{1}.metrics(:,1),'o-',1:30,R{2}.metrics(:,1),'o-');xline(10.5,'k--');ylabel('Theta MAE (rad)');legend('Trial2 adapting','Trial3 freeze');grid on
nexttile;plot(1:30,R{1}.metrics(:,2),'o-',1:30,R{2}.metrics(:,2),'o-');xline(10.5,'k--');ylabel('Phi MAE (rad)');xlabel('Circle number');grid on
exportgraphics(fig,fullfile(root,'tmp/trial3_freeze_comparison.png'),'Resolution',150);close(fig);diary off
function [t,x]=read(ds,name)
names=ds.getElementNames;v=ds{find(strcmp(names,name),1)}.Values;t=v.Time(:);raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
