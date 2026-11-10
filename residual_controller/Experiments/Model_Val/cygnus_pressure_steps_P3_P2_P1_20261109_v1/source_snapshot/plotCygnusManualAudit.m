function plotCygnusManualAudit(folder,description)
%PLOTCYGNUSMANUALAUDIT Saved protocol data only; no model or device access.
if nargin<2,description='Manual motion without actuation';end
T=readtable(fullfile(folder,'tracking.csv'));t=T.elapsed_s;
S=jsondecode(fileread(fullfile(folder,'summary.json')));
fig=figure('Visible','off','Color','w','Position',[20 20 1350 1000]);tiledlayout(3,1,'TileSpacing','compact');
firstAx=nexttile;plot(t,[T.tip_x_mm,T.tip_y_mm,T.tip_z_mm]);grid on;
ylabel('Tip raw XYZ (mm)');legend('X','Y','Z','Location','best');title('Raw sensor positions, no interpolation');
nexttile;plot(t(2:end),diff(t),'k.-');hold on;yline(1/S.transform_update_hz,'b--');grid on;
ylabel('Reply frame interval (s)');legend('Frame-counter-derived interval','Nominal update period','Location','best');title('Intervals between reported frame counters; not GUI-refresh timestamps');
nexttile;stairs(t,double(T.base_status==1),'b');hold on;stairs(t,double(T.tip_status==1),'r--');grid on;ylim([-.05 1.05]);
missing=T.startup==0&T.tip_status==2;scatter(t(missing),zeros(nnz(missing),1),15,'r','filled');
ylabel('Transform valid');legend('0A base','0B tip','Runtime tip missing','Location','best');xlabel('Seconds from first final-session BX reply');
title(sprintf('Runtime MISSING reports: base %d, tip %d',S.runtime_base_missing,S.runtime_tip_missing));
for ax=findall(fig,'Type','axes')',xlim(ax,[0 t(end)]);end
groupsFile=fullfile(folder,'movement_groups.csv');
if isfile(groupsFile)
    G=readtable(groupsFile);bounds=G.group_start_s(2:end);
    for ax=findall(fig,'Type','axes')',for boundary=bounds',xline(ax,boundary,':','Color',[.4 .4 .4]);end,end
    yl=ylim(firstAx);
    for j=1:height(G)
        text(firstAx,(G.response_start_s(j)+G.response_end_s(j))/2,yl(2)-.08*diff(yl),string(G.PMA(j)), ...
            'HorizontalAlignment','center','FontWeight','bold');
    end
end
sgtitle(sprintf('%s: %d runtime tip-MISSING reports',description,S.runtime_tip_missing));
exportgraphics(fig,fullfile(folder,'manual_tracking.png'),'Resolution',180);exportgraphics(fig,fullfile(folder,'manual_tracking.pdf'),'ContentType','vector');savefig(fig,fullfile(folder,'manual_tracking.fig'));close(fig);
end
