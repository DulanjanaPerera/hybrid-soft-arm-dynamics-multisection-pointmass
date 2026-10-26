root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_gainTesting_v7_analysis');
S=load(fullfile(fileparts(folder),'1rad_10s_gainTesting_v7.mat'));ds=S.out.logsout;
[t,th]=read(ds,'parameterEstimate');D=struct;
for name={'adaptationGamma','modelValid','orientationValid','thetaObservable','theta','phi','limited','filtered_d_length','d_length','adaptationActive'}
 [tt,x]=read(ds,name{1});D.(name{1})=interp1(tt,x,t,'previous','extrap');
end
fid=fopen(fullfile(folder,'gamma_stage_audit.txt'),'w');fclose(fid);diary(fullfile(folder,'gamma_stage_audit.txt'));
for stamp=[80 114.6 126.1 170.15 233.15 248.85 254.95 272.6 281.3 301.1 319.95]
 i=find(t>=stamp,1);fprintf('TIME%.3f d%s h%s gamma%s\n',t(i),mat2str(th(i,4:5),6),mat2str(th(i,6:8),6),mat2str(D.adaptationGamma(i,:),5));
end
good=D.modelValid&D.orientationValid&D.thetaObservable;
rate=[NaN;atan2(sin(diff(D.theta)),cos(diff(D.theta)))./diff(t)];pair=good&[false;good(1:end-1)];rate(~pair)=NaN;
for win=[180 220 233.15 254.95 272.6 281.3 301.1;220 233.15 248.85 272.6 281.3 301.1 320]
 ix=t>=win(1)&t<win(2)&good;ir=t>=win(1)&t<win(2)&pair;allix=t>=win(1)&t<win(2);
 fprintf('WIN%.2f-%.2f dstart%s dend%s thetaRateSD%.5f slowPct%.3f phiMax%.5f l3velocityPeakmm/s%.3f limitPct%.3f activePct%.3f\n',win(1),win(2),mat2str(th(find(allix,1),4:5),5),mat2str(th(find(allix,1,'last'),4:5),5),std(rate(ir)),100*mean(rate(ir)<.2*2*pi/10),max(D.phi(ix)),1000*max(abs(D.filtered_d_length(ix,2))),100*mean(D.limited(allix)),100*mean(D.adaptationActive(allix)));
end
fig=figure('Visible','off','Color','w','Position',[50 50 1400 1100]);tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
nexttile;stairs(t,D.adaptationGamma(:,1),'LineWidth',1.5);ylabel('\gamma_k');
nexttile;stairs(t,D.adaptationGamma(:,6),'LineWidth',1.5);ylabel('\gamma_h');
nexttile;stairs(t,D.adaptationGamma(:,4),'LineWidth',1.5);set(gca,'YScale','log');ylabel('\gamma_d (log scale)');
nexttile;plot(t,th(:,4:5),'LineWidth',1.5);ylabel('Estimated d (N s/m)');xlabel('Time (s)');legend('d_2','d_3','Location','northwest');
for ax=findall(fig,'Type','axes')',set(ax,'FontName','Arial','FontSize',16,'XGrid','on','YGrid','on');xlim(ax,[0 330]);xline(ax,254.95,'k--','LineWidth',1.2,'HandleVisibility','off');end
exportgraphics(fig,fullfile(folder,'gamma_and_damping_history.png'),'Resolution',200);exportgraphics(fig,fullfile(folder,'gamma_and_damping_history.pdf'),'ContentType','vector');savefig(fig,fullfile(folder,'gamma_and_damping_history.fig'));close(fig);diary off
function [t,x]=read(ds,name)
ix=find(strcmp(ds.getElementNames,name),1);v=ds{ix}.Values;t=double(v.Time(:));raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
