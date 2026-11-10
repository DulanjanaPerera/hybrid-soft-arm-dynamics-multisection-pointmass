function r=plotCygnusSineAudit(folder)
%PLOTCYGNUSSINEAUDIT Offline raw pose + inferred standard-sine timing.
% No pressure data exist here: fitted template is NOT an observed command.
T=readtable(fullfile(folder,'tracking.csv'));t=T.elapsed_s;
b=[T.base_qw,T.base_qx,T.base_qy,T.base_qz];p=[T.tip_qw,T.tip_qx,T.tip_qy,T.tip_qz];
q=[sum(b.*p,2),b(:,1).*p(:,2:4)-p(:,1).*b(:,2:4)-cross(b(:,2:4),p(:,2:4),2)];
q=q./vecnorm(q,2,2);ok=T.base_status==1&T.tip_status==1;
ref=median(q(ok&t>=1&t<=5,:),1);ref=ref/norm(ref);
% Total relative orientation change from a straight baseline; NOT calibrated
% continuum-arm phi. No gap interpolation, no sensor-offset calibration.
angle=2*acos(min(max(abs(q*ref'),0),1));angle(~ok)=NaN;
mask=ok&T.startup==0&isfinite(angle);
offsets=0:.025:40;loss=inf(size(offsets));coefficients=zeros(numel(offsets),4);
for j=1:numel(offsets)
    P=template(t-offsets(j));A=[ones(nnz(mask),1),P(mask,:)];c=A\angle(mask);
    if all(c(2:4)>=0),loss(j)=mean((A*c-angle(mask)).^2);coefficients(j,:)=c';end
end
[best,j]=min(loss);offset=offsets(j);P=template(t-offset);
missing=T.startup==0&T.tip_status==2;
assigned=sum(P(missing,:),2)>0&all(P(missing,1:2)==0,2)&P(missing,3)>0;
period=2*pi/.8;trial=3*period+10;
r=struct('method','Least-squares alignment of standard 3-channel sine timing to raw relative orientation change', ...
    'referenceSource','RL_data_collection/rlPressureSineBlock.m; user confirms standard block', ...
    'offset_s',offset,'orientationFitRMSE_rad',sqrt(best),'fitCoefficients',coefficients(j,:), ...
    'P3TemplateWindows_elapsed_s',[offset+5+4*trial,offset+5+4*trial+3*period;offset+5+5*trial,offset+5+5*trial+3*period], ...
    'runtimeTipMissing',nnz(missing),'missingWithinInferredP3Excitation',nnz(assigned), ...
    'caveat','Motion-fitted offset includes actuator response lag; no actual pressure timestamps or values were recorded. Channel assignment is sequence-based inference, not synchronized proof.');
f=fopen(fullfile(folder,'sine_alignment.json'),'w');fwrite(f,jsonencode(r,PrettyPrint=true));fclose(f);
writetable(table(t,angle,P,missing,'VariableNames',{'elapsed_s','relativeOrientationChange_rad','inferredPressureTemplate_bar','tipMissing'}),fullfile(folder,'motion_alignment.csv'));
fig=figure('Visible','off','Color','w','Position',[20 20 1400 1050]);tiledlayout(4,1,'TileSpacing','compact');
nexttile;plot(t,angle,'b');ylabel('Orientation change (rad)');title('Raw tip-to-base relative orientation change; not calibrated phi');
nexttile;plot(t,P,'LineWidth',1);ylabel('Template (bar)');legend('P1 template','P2 template','P3 template','Location','best');title('Standard sine sequence, offset fitted to motion; NOT measured pressure');
nexttile;stairs(t,T.base_status==1,'b');hold on;stairs(t,double(T.tip_status==1),'r--');scatter(t(missing),zeros(nnz(missing),1),15,'r','filled');ylabel('Tracker valid');ylim([-.1 1.1]);legend('0A/base','0B/tip','Runtime tip missing','Location','best');
nexttile;plot(t,[T.tip_x_mm,T.tip_y_mm,T.tip_z_mm]);ylabel('Tip raw XYZ (mm)');xlabel('Seconds from first BX reply in final tracking session');legend('X','Y','Z','Location','best');
for ax=findall(fig,'Type','axes')',grid(ax,'on');xlim(ax,[0 t(end)]);end
sgtitle('Cygnus log: direct device-reported tip losses in inferred P3 motion group');
exportgraphics(fig,fullfile(folder,'tracking_and_sine_sequence.png'),'Resolution',180);
exportgraphics(fig,fullfile(folder,'tracking_and_sine_sequence.pdf'),'ContentType','vector');savefig(fig,fullfile(folder,'tracking_and_sine_sequence.fig'));close(fig);disp(r);
end
function P=template(t)
P=zeros(numel(t),3);excitation=3*2*pi/.8;trial=excitation+10;elapsed=t-5;index=floor(elapsed/trial);local=elapsed-index*trial;
valid=elapsed>=0&index<6&local<excitation;
for channel=1:3,m=valid&floor(index/2)+1==channel;P(m,channel)=1.5*(1-cos(.8*local(m)));end
end
