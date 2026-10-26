root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_gainTesting_v5_analysis');
fid=fopen(fullfile(folder,'v4_v5_roughness.txt'),'w');fclose(fid);diary(fullfile(folder,'v4_v5_roughness.txt'));
R=cell(2,1);
for run=4:5
 S=load(fullfile(root,sprintf('passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_gainTesting_v%d.mat',run)));ds=S.out.logsout;
 [t,q]=read(ds,'meas_length');D=struct;
 for name={'sampled des_length','filtered_d_length','d_length','theta','des_theta','phi','des_phi','orientationValid','thetaObservable','modelValid','status','limited','pressureCommand_bar','adaptationGamma','parameterEstimate'}
  [tt,x]=read(ds,name{1});if numel(tt)==1,x=repmat(x,numel(t),1);else,x=interp1(tt,x,t,'previous','extrap');end;D.(matlab.lang.makeValidName(name{1}))=x;
 end
 ok=D.orientationValid&D.thetaObservable&D.modelValid;
 thetaRate=[NaN;atan2(sin(diff(D.theta)),cos(diff(D.theta)))./diff(t)];
 pairGood=ok&[false;ok(1:end-1)]&D.status==3;thetaRate(~pairGood)=NaN;
 fprintf('TRIAL%d\n',run);
 for range=[40 290;70 320]
  ix=t>=range(1)&t<range(2)&ok;dv=D.filtered_d_length-D.d_length;
  % Residual from a fitted fundamental quantifies nonsinusoidal content;
  % no fitted signals are substituted into tracking or velocity metrics.
  X=[ones(nnz(ix),1),sin(2*pi*t(ix)/10),cos(2*pi*t(ix)/10)];
  qResidual=q(ix,2:3)-X*(X\q(ix,2:3));
  ir=t>=range(1)&t<range(2)&pairGood;rates=thetaRate(ir);
  fprintf('WIN%d-%d velocityRMSEmm/s%s velocityPeakmm/s%s lengthFundamentalResidualmm%s thetaRateSD%.6f thetaRateMinMax%s thetaSlowPct%.3f phiMinMax%s phiAbove1.05Pct%.3f limitingPct%.3f\n',range(1),range(2),mat2str(1000*sqrt(mean(dv(ix,:).^2)),6),mat2str(1000*max(abs(D.filtered_d_length(ix,:))),6),mat2str(1000*sqrt(mean(qResidual.^2)),6),std(rates),mat2str([min(rates) max(rates)],6),100*mean(rates<.2*(2*pi/10)),mat2str([min(D.phi(ix)) max(D.phi(ix))],6),100*mean(D.phi(ix)>1.05),100*mean(D.limited(t>=range(1)&t<range(2))));
 end
 R{run-3}=struct('t',t,'q',q,'D',D,'thetaRate',thetaRate);
end
fig=figure('Visible','off','Color','w','Position',[60 60 1400 1100]);tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
for panel=1:4
 nexttile;hold on
 for r=1:2
  A=R{r};ix=A.t>=290&A.t<=320;t=A.t(ix);valid=A.D.modelValid(ix)&A.D.orientationValid(ix);
  if panel<=2,y=1000*A.q(ix,panel+1);elseif panel==3,y=1000*A.D.filtered_d_length(ix,2);else,y=A.thetaRate(ix);end
  y(~valid)=NaN;plot(t,y,'LineWidth',1.3);
 end
 A=R{2};ix=A.t>=290&A.t<=320;t=A.t(ix);
 if panel<=2,y=1000*A.D.sampledDes_length(ix,panel);elseif panel==3,y=1000*A.D.d_length(ix,2);else,y=repmat(2*pi/10,numel(t),1);end
 plot(t,y,'k--','LineWidth',1.5);
 labels={'\Delta l_2 (mm)','\Delta l_3 (mm)','l_3 velocity (mm/s)','Apparent \theta rate (rad/s)'};ylabel(labels{panel});
 if panel==1,legend('v4: lambda=[28,25.5]','v5: lambda=[18,16.5]','Desired','Location','northeast','Orientation','horizontal');end
 if panel==4,xlabel('Time (s)');end
 grid on;xlim([290 320]);set(gca,'FontName','Arial','FontSize',15);
end
exportgraphics(fig,fullfile(folder,'v4_v5_last3cycles.png'),'Resolution',200);savefig(fig,fullfile(folder,'v4_v5_last3cycles.fig'));close(fig);diary off
function [t,x]=read(ds,name)
ix=find(strcmp(ds.getElementNames,name),1);v=ds{ix}.Values;t=double(v.Time(:));raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
