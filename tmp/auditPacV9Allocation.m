root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'passivity_adaptive_control'));
folder=fullfile(root,'passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_gainTesting_v9_analysis');
S=load(fullfile(fileparts(folder),'1rad_10s_gainTesting_v9.mat'));ds=S.out.logsout;
[t,th]=read(ds,'parameterEstimate');D=struct;
for name={'pressureOptions','status','modelValid','limited','adaptationActive','tauBlended_N','forceBias_N','pressureCommand_bar','allocationResidual_N','parameterLower','parameterUpper','parameterMaxRate','des_theta','des_phi','sampled des_length','meas_length','theta','phi','orientationValid','thetaObservable'}
 [tt,x]=read(ds,name{1});D.(matlab.lang.makeValidName(name{1}))=interp1(tt,x,t,'previous','extrap');
end
p=D.pressureOptions(1,:);cfg=pacPressureConfig();fields={'common_bar','sTrim_bar','aTrim_bar','minimum_bar','maximum_bar','maxDeltaS_bar','maxDeltaA_bar','slew_barPerS','activation_bar','forcePerBar_N'};
for j=1:10,cfg.(fields{j})=p(j);end;cfg.enableCommonAdaptation=true;
target=NaN(numel(t),3);ar=NaN(numel(t),2);al=false(size(t));slew=false(size(t));
circle=t>=20&t<80;valid=circle&D.modelValid&D.status==3;
for i=find(valid)'
 [pp,~,rr,ll]=pacAllocateVariableCommon(D.tauBlended_N(i,:)',D.forceBias_N(i,:)',cfg);
 target(i,:)=pp';ar(i,:)=rr';al(i)=ll;slew(i)=any(abs(pp'-D.pressureCommand_bar(i,:))>1e-9);
end
fid=fopen(fullfile(folder,'allocation_and_bounds.txt'),'w');fclose(fid);diary(fullfile(folder,'allocation_and_bounds.txt'));
fprintf('REFERENCE median angular speed%g period%g phi%s\n',median(diff(D.des_theta(circle))./diff(t(circle))),2*pi/median(diff(D.des_theta(circle))./diff(t(circle))),mat2str([min(D.des_phi(circle)) max(D.des_phi(circle))]));
fprintf('CIRCLE total%d validControl%d loggedLimited%d replayAllocationLimited%d slewTargetMismatch%d both%d ORmismatch%d\n',nnz(circle),nnz(valid),nnz(circle&D.limited),nnz(al),nnz(slew),nnz(al&slew),nnz(valid&(D.limited~=(al|slew))));
fprintf('allocationResidual peak%g commandResidual peak%g N\n',max(vecnorm(ar(valid,:),2,2)),max(vecnorm(D.allocationResidual_N(valid,:),2,2)));
fprintf('PRESSURE min%s max%s targetMin%s targetMax%s\n',mat2str(min(D.pressureCommand_bar(circle,:)),6),mat2str(max(D.pressureCommand_bar(circle,:)),6),mat2str(min(target(valid,:)),6),mat2str(max(target(valid,:)),6));
rates=diff(D.pressureCommand_bar)./diff(t);ix=circle(1:end-1)&circle(2:end);
fprintf('Pressure peakRate%s exact2barPerSecPct%s\n',mat2str(max(abs(rates(ix,:))),6),mat2str(100*mean(abs(abs(rates(ix,:))-2)<1e-7),6));
labels={'k1','k2','k3','d2','d3','h1','h2','h3','b2','b3'};
for j=1:10
 lo=th(:,j)<=D.parameterLower(:,j)+1e-9;hi=th(:,j)>=D.parameterUpper(:,j)-1e-9;
 fprintf('BOUND %s countLower%d countUpper%d circleMin%g circleMax%g\n',labels{j},nnz(circle&lo),nnz(circle&hi),min(th(circle,j)),max(th(circle,j)));
end
rr=diff(th)./diff(t);ix=circle(1:end-1)&D.adaptationActive(1:end-1);
fprintf('Activecircle rateCapcounts%s peaks%s\n',mat2str(sum(abs(rr(ix,:))>=D.parameterMaxRate(ix,:)-1e-6)),mat2str(max(abs(rr(ix,:))),6));
fprintf('AT endTracking parameters%s\n',mat2str(th(find(circle,1,'last'),:),6));
for win=[20 40 60 74;40 60 80 80]
 ix=t>=win(1)&t<win(2);fprintf('WIN%d-%d activePct%g limitedPct%g thetaMAE%g phiMAE%g\n',win(1),win(2),100*mean(D.adaptationActive(ix)),100*mean(D.limited(ix)),mean(abs(atan2(sin(D.des_theta(ix)-D.theta(ix)),cos(D.des_theta(ix)-D.theta(ix))))),mean(abs(D.des_phi(ix)-D.phi(ix))));
end
fig=figure('Visible','off','Color','w','Position',[60 60 1400 1100]);tiledlayout(4,1,'TileSpacing','compact','Padding','compact');ix=t>=74&t<80;
nexttile;plot(t(ix),1000*D.sampledDes_length(ix,2),'k--',t(ix),1000*D.meas_length(ix,3),'LineWidth',1.5);ylabel('\Delta l_3 (mm)');legend('Desired','Measured','Location','northeast');
for j=1:2,nexttile;plot(t(ix),target(ix,j),'k--',t(ix),D.pressureCommand_bar(ix,j),'LineWidth',1.5);ylabel(sprintf('P_%d (bar)',j));if j==1,legend('Allocator target','Command after slew','Location','northeast');end;end
nexttile;stairs(t(ix),D.adaptationActive(ix),'LineWidth',1.5);ylabel('Adaptation active');xlabel('Time (s)');ylim([-.1 1.2]);
for a=findall(fig,'Type','axes')',set(a,'FontName','Arial','FontSize',16,'XGrid','on','YGrid','on');xlim(a,[74 80]);end
exportgraphics(fig,fullfile(folder,'two_second_slew_diagnostics.png'),'Resolution',200);savefig(fig,fullfile(folder,'two_second_slew_diagnostics.fig'));close(fig);diary off
function [t,x]=read(ds,name)
ix=find(strcmp(ds.getElementNames,name),1);v=ds{ix}.Values;t=double(v.Time(:));raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
