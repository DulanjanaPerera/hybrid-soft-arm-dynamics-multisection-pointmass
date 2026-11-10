function plotCygnusManualAudit(folder)
%PLOTCYGNUSMANUALAUDIT Saved protocol data only; no model or device access.
T=readtable(fullfile(folder,'tracking.csv'));t=T.elapsed_s;
fig=figure('Visible','off','Color','w','Position',[20 20 1350 1000]);tiledlayout(3,1,'TileSpacing','compact');
nexttile;plot(t,[T.tip_x_mm,T.tip_y_mm,T.tip_z_mm]);grid on;
ylabel('Tip raw XYZ (mm)');legend('X','Y','Z','Location','best');title('Manual motion: raw sensor positions, no interpolation');
nexttile;plot(t(2:end),diff(t),'k.-');hold on;yline(1/40.033,'b--');grid on;
ylabel('Reply frame interval (s)');legend('Frame-counter-derived interval','Nominal25ms update','Location','best');title('Two175ms gaps occur only near the final stop');
nexttile;stairs(t,double(T.base_status==1),'b');hold on;stairs(t,double(T.tip_status==1),'r--');grid on;ylim([-.05 1.05]);
ylabel('Transform valid');legend('0A base','0B tip','Location','best');xlabel('Seconds from first final-session BX reply');title('Only initial acquisition is invalid; all2444subsequent paired transforms valid');
for ax=findall(fig,'Type','axes')',xlim(ax,[0 t(end)]);end
sgtitle('Manual fast motion, no pressure: no runtime MISSING sensor reports');
exportgraphics(fig,fullfile(folder,'manual_tracking.png'),'Resolution',180);exportgraphics(fig,fullfile(folder,'manual_tracking.pdf'),'ContentType','vector');savefig(fig,fullfile(folder,'manual_tracking.fig'));close(fig);
end
