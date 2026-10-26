root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'passivity_adaptive_control'));
file=fullfile(root,'passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_gainTesting_v3.mat');
folder=fullfile(fileparts(file),'1rad_10s_gainTesting_v3_analysis');if ~isfolder(folder),mkdir(folder);end
fid=fopen(fullfile(folder,'audit.txt'),'w');fclose(fid);diary(fullfile(folder,'audit.txt'));
S=load(file);ds=S.out.logsout;[t,th]=read(ds,'parameterEstimate');N=numel(t);D=struct;
names={'enableAdaptation','enableCommonPressureAdaptation','freezeHysteresisWeights','adaptationGamma','adaptationOptions','hystoryOptions','parameterInitial','parameterLower','parameterMaxRate','parameterScale','parameterUpper','pressureOptions','supervisorOptions','modelValid','status','limited','adaptationActive','parameterAtBound','s_mPerS','measurementTime_s','configurationValid','baseOreintationValid','sensorValid','allocationResidual_N','historyState','lambda_gains','damping_gains','des_theta','des_phi','theta','phi','orientationValid','thetaObservable','sampled des_length','meas_length','pressureCommand_bar','commonPressure_bar','hysteresisForce_N'};
for k=1:numel(names)
 [tt,x]=read(ds,names{k});if numel(tt)==1,x=repmat(x,N,1);else,x=interp1(tt,x,t,'previous','extrap');end
 D.(matlab.lang.makeValidName(names{k}))=x;
end
fprintf('SOURCE %s N%d TIME%g-%g nominal dt%s\n',file,N,t(1),t(end),mat2str(unique(round(diff(t),8))'));
tuning={'lambda_gains','damping_gains','adaptationGamma','parameterMaxRate','parameterLower','parameterUpper','parameterScale','adaptationOptions','hystoryOptions','pressureOptions','supervisorOptions','freezeHysteresisWeights','enableAdaptation','enableCommonPressureAdaptation','parameterInitial'};
events=[];
for k=1:numel(tuning)
 x=D.(tuning{k});same=x(2:end,:)==x(1:end-1,:) | (isnan(x(2:end,:))&isnan(x(1:end-1,:)));
 ii=[1;find(any(~same,2))+1];events=[events;t(ii(2:end))]; %#ok<AGROW>
 fprintf('%s changes%d [time values]\n',tuning{k},numel(ii)-1);disp([t(ii) x(ii,:)]);
end
events=unique(events);sn=vecnorm(D.s_mPerS,2,2);rn=vecnorm(D.allocationResidual_N,2,2);
age=t-D.measurementTime_s;
good=D.configurationValid & D.baseOreintationValid & isfinite(D.measurementTime_s) & age>=-1e-9 & age<=D.supervisorOptions(:,5)+1e-9;
dt=zeros(N,1);last=-1;
for i=1:N
 if good(i)&&D.measurementTime_s(i)>last
  if last>=0,dt(i)=D.measurementTime_s(i)-last;end
  last=D.measurementTime_s(i);
 end
end
% Disjoint first-failing gate classification, not overlapping counts.
reason=zeros(N,1);checks={~D.enableAdaptation,~good | ~D.modelValid,D.status~=3,D.limited,~isfinite(rn)|rn>D.adaptationOptions(:,2),dt<=0|dt>D.adaptationOptions(:,3),sn<=D.adaptationOptions(:,1)};
labels={'enabled update','manual disable','invalid/freshness/model','supervisor not control','allocation/slew limiting','nominal residual tolerance','no fresh/acceptable update interval','s deadband','zero projected rate/unexplained'};
for k=1:numel(checks),ix=reason==0&checks{k};reason(ix)=k;end
reason(reason==0&~D.adaptationActive)=8;
mismatch=nnz(reason~=0 & D.adaptationActive);fprintf('ACTIVE vs gates mismatches%d\n',mismatch);
circle=t>=20 & t<320;
for mask={true(N,1),circle}
 ix=mask{1};fprintf('GATES window%g-%g N%d inactive%d activePct%.4f\n',min(t(ix)),max(t(ix)),nnz(ix),nnz(ix&~D.adaptationActive),100*mean(D.adaptationActive(ix)));
 for k=0:8,fprintf('%d %s %d\n',k,labels{k+1},nnz(ix&reason==k));end
end
fprintf('INITIAL%s FINAL%s MIN%s MAX%s\n',mat2str(th(1,:),7),mat2str(th(end,:),7),mat2str(min(th),7),mat2str(max(th),7));
eigmin=(th(:,1)+th(:,2)+2*th(:,3)-sqrt((th(:,1)-th(:,2)).^2+4*th(:,3).^2))/2;
fprintf('K minimum eigenvalue min%g nonpositiveCircle%d\n',min(eigmin),nnz(circle&eigmin<=0));
neg=circle&eigmin<=0;ii=find(neg&~[false;neg(1:end-1)]);jj=find(neg&~[neg(2:end);false]);
fprintf('K nonpositive events [start end count]\n');disp([t(ii) t(jj) jj-ii+1]);
fprintf('K circle minEigen%g validNonpositive%d\n',min(eigmin(circle)),nnz(neg&D.modelValid));
fprintf('REFERENCE phi range%s theta range%s\n',mat2str([min(D.des_phi) max(D.des_phi)]),mat2str([min(D.des_theta) max(D.des_theta)]));
z=~D.modelValid;fprintf('INVALID samples%d events%d tipInvalid%d status4%d faults%d limitedCircle%d\n',nnz(z),nnz(z&~[false;z(1:end-1)]),nnz(z&~D.sensorValid(:,2)),nnz(D.status==4),nnz(D.status==5),nnz(circle&D.limited));
% Rates are based on next-sample stored estimates and actual accepted dt.
rates=diff(th)./dt(1:end-1);validRate=D.adaptationActive(1:end-1)&dt(1:end-1)>0;
caps=D.parameterMaxRate(1:end-1,:);atcap=abs(rates)>=caps-1e-6 & caps>0;
fprintf('PEAK acceptedRates%s\n',mat2str(max(abs(rates(validRate,:))),6));
fprintf('RATECAP percent of active samples%s\n',mat2str(100*mean(atcap(validRate,:)),5));
eth=atan2(sin(D.des_theta-D.theta),cos(D.des_theta-D.theta));eph=D.des_phi-D.phi;el=D.sampledDes_length-D.meas_length(:,2:3);
fprintf('LAPS lap thetaMAE phiMAE lengthsRMSE activePct deadbandPct dtGatePct k2mean k2range hmean\n');
metrics=[];
for lap=1:floor((min(t(end),320)-20)/10)
 ix=t>=20+(lap-1)*10&t<20+lap*10;ok=ix&D.orientationValid&D.thetaObservable&D.modelValid;
 m=[lap mean(abs(eth(ok))) mean(abs(eph(ok))) 1000*sqrt(mean(el(ok,:).^2)) 100*mean(D.adaptationActive(ix)) 100*mean(reason(ix)==7)];metrics=[metrics;m]; %#ok<AGROW>
 fprintf('%d %.6f %.6f %s %.2f %.2f %.2f %.4f %s %s\n',lap,m(2),m(3),mat2str(m(4:5),5),m(6),m(7),100*mean(reason(ix)==6),mean(th(ix,2)),mat2str([min(th(ix,2)) max(th(ix,2))],5),mat2str(mean(th(ix,6:8)),5));
end
writematrix(metrics,fullfile(folder,'lap_metrics.csv'));
fprintf('DEADBAND off sn range%s threshold range%s\n',mat2str([min(sn(reason==7)) max(sn(reason==7))],7),mat2str([min(D.adaptationOptions(:,1)) max(D.adaptationOptions(:,1))]));
% For each combined setting interval, compare complete circle laps only.
edges=unique([20;events(events>20&events<320);320]);
for k=1:numel(edges)-1
 a=edges(k);b=edges(k+1);ix=t>=a&t<b;ok=ix&D.modelValid&D.orientationValid&D.thetaObservable;
 complete=metrics(20+(metrics(:,1)-1)*10>=a & 20+metrics(:,1)*10<=b,1);
 fprintf('INTERVAL %.3f-%.3f theta%.6f phi%.6f activePct%.3f limiting%d invalid%d dtGate%d completeLaps%s\n',a,b,mean(abs(eth(ok))),mean(abs(eph(ok))),100*mean(D.adaptationActive(ix)),nnz(ix&reason==4),nnz(ix&reason==2),nnz(ix&reason==6),mat2str(complete'));
end
fig=figure('Visible','off','Color','w','Position',[60 60 1400 1050]);tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
nexttile;semilogy(t,max(sn,1e-5),'LineWidth',1);hold on;semilogy(t,D.adaptationOptions(:,1),'k--','LineWidth',1.6);ylabel('||s|| (m/s)');legend('Filtered error','Deadband','Location','northeast');
nexttile;stairs(t,reason,'LineWidth',1);ylabel('Update gate');yticks(0:8);yticklabels({'Update','Disabled','Invalid','Supervisor','Limited','Residual','dt','Deadband','Other'});
nexttile;plot(t,th(:,2),'LineWidth',1.4);hold on;plot(t,D.parameterLower(:,2),'k--');ylabel('k_2 (N/m)');
nexttile;plot(metrics(:,1),metrics(:,2:3),'o-','LineWidth',1.4);ylabel('Angle MAE (rad)');xlabel('Circle number');legend('\theta','\phi','Location','northeast');
ax=findall(fig,'Type','axes');for a=ax',set(a,'FontName','Arial','FontSize',16,'XGrid','on','YGrid','on','GridAlpha',.15);end
exportgraphics(fig,fullfile(folder,'adaptation_gate_diagnostics.png'),'Resolution',200);savefig(fig,fullfile(folder,'adaptation_gate_diagnostics.fig'));close(fig);
h=plotPacAdaptationPresentation(file,folder);close(h);
diary off
function [t,x]=read(ds,name)
ix=find(strcmp(ds.getElementNames,name),1);assert(~isempty(ix),'Missing %s',name);v=ds{ix}.Values;t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
