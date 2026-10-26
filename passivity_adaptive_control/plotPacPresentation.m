function figures=plotPacPresentation(filePath,timeWindow,outputDir,prefix)
%PLOTPACPRESENTATION Export saved PAC tracking data; no models/hardware run.
% XYZ is in the existing arm-frame sensor convention. Reference includes
% 55 mm base and 50 mm tip-sensor offsets, not a fitted alignment to data.
% Export PNG (300 dpi), vector PDF/SVG and editable FIG for each figure.
if nargin<4,prefix='tracking';end
S=load(filePath,'out');ds=S.out.logsout;
[t,qd]=read(ds,'sampled des_length');
if nargin<2 || isempty(timeWindow),timeWindow=[t(1) t(end)];end
if nargin<3 || isempty(outputDir),outputDir=fullfile(fileparts(filePath),'presentation_figures');end
if ~isfolder(outputDir),mkdir(outputDir);end
sel=t>=timeWindow(1)&t<=timeWindow(2);assert(nnz(sel)>1,'Empty time window.');
td=align(ds,'des_theta',t);pd=align(ds,'des_phi',t);
tm=align(ds,'theta',t);pm=align(ds,'phi',t);
lengths=align(ds,'meas_length',t);pressure=align(ds,'pressureCommand_bar',t);
xyz=align(ds,'position',t);
ov=align(ds,'orientationValid',t)>.5;cv=align(ds,'configurationValid',t)>.5;
obs=align(ds,'thetaObservable',t)>.5;
lengths(~ov,:)=NaN;pm(~ov)=NaN;tm(~(ov&obs))=NaN;xyz(~cv,:)=NaN;
desiredLength=[-sum(qd,2) qd];
L=[.055 .17411 .05];ratio=ones(size(pd));rho=(L(2)/2+L(3))*pd;
nz=abs(pd)>1e-8;ratio(nz)=sin(pd(nz))./pd(nz);
rho(nz)=L(2)*2*sin(pd(nz)/2).^2./pd(nz)+L(3)*sin(pd(nz));
xyzDesired=[rho.*cos(td),rho.*sin(td),L(1)+L(2)*ratio+L(3)*cos(pd)];
blue=[.05 .37 .70];orange=[.85 .33 .10];purple=[.49 .18 .56];
figures=gobjects(5,1);suffix={'xyz_tracking','pressure_trajectories','length_tracking','theta_phi_tracking','xyz_path'};
% Fixed white background and large text for presentation use.
for f=1:5
 figures(f)=figure('Visible','off','Color','w','Position',[80 80 1400 850], ...
  'Name',suffix{f},'NumberTitle','off');
 if f==1
  tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
  for j=1:3
   nexttile;pair(t(sel),1000*xyzDesired(sel,j),1000*xyz(sel,j),blue);
   axisNames={'X (mm)','Y (mm)','Z (mm, down)'};ylabel(axisNames{j});
   if j==1,legend('Desired','Measured','Location','northeast','Orientation','horizontal');end
   if j==3,xlabel('Time (s)');end
  end
 elseif f==2
  ax=axes;plot(t(sel),pressure(sel,:),'LineWidth',1.8);colororder(ax,[blue;orange;purple]);
  ylabel('Pressure command (bar)');xlabel('Time (s)');legend('P_1','P_2','P_3','Location','northeast','Orientation','horizontal');
 elseif f==3
  tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
  for j=1:3
   nexttile;pair(t(sel),1000*desiredLength(sel,j),1000*lengths(sel,j),blue);
   ylabel(sprintf('\\Delta l_%d (mm)',j));
   if j==1,legend('Desired','Measured','Location','northeast','Orientation','horizontal');end
   if j==3,xlabel('Time (s)');end
  end
 elseif f==4
  tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
  nexttile;pair(t(sel),wrapped(td(sel)),wrapped(tm(sel)),blue);ylabel('\theta (rad)');ylim([-pi pi]);yticks([-pi 0 pi]);yticklabels({'-\pi','0','\pi'});
  legend('Desired','Measured','Location','northeast','Orientation','horizontal');
  nexttile;pair(t(sel),pd(sel),pm(sel),blue);ylabel('\phi (rad)');xlabel('Time (s)');
 elseif f==5
  plot3(1000*xyzDesired(sel,1),1000*xyzDesired(sel,2),1000*xyzDesired(sel,3),'k--','LineWidth',2);hold on
  plot3(1000*xyz(sel,1),1000*xyz(sel,2),1000*xyz(sel,3),'Color',blue,'LineWidth',1.6);
  xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm, down)');axis equal;view(40,24);set(gca,'ZDir','reverse');
  legend('Desired','Measured','Location','northeast');
 end
 axesList=findall(figures(f),'Type','axes');
 for ax=axesList'
  set(ax,'FontName','Arial','FontSize',18,'LineWidth',1,'Box','on','XGrid','on','YGrid','on','GridAlpha',.13,'TickDir','out');
  if f~=5
   xlim(ax,timeWindow);
   if ~strcmp(ax.YLabel.String,'\theta (rad)')
    bounds=ylim(ax);pad=.08*diff(bounds);ylim(ax,bounds+[-pad pad]);
   end
  end
 end
 drawnow;
 base=fullfile(outputDir,[prefix '_' suffix{f}]);
 exportgraphics(figures(f),[base '.png'],'Resolution',300);
 exportgraphics(figures(f),[base '.pdf'],'ContentType','vector');
 print(figures(f),[base '.svg'],'-dsvg');savefig(figures(f),[base '.fig']);
 fprintf('EXPORTED %s\n',base);
end
fid=fopen(fullfile(outputDir,[prefix '_notes.txt']),'w');assert(fid>=0);
fprintf(fid,'Source: %s\nWindow: %.3f to %.3f s\n',filePath,timeWindow);
fprintf(fid,'XYZ: logged position in arm frame; desired includes base .055 m, active length .17411 m, tip sensor .05 m. No fitted shift, phase alignment or smoothing.\n');
fprintf(fid,'Lengths: orientation-inferred length CHANGES; not independent strain measurements. Pressure: commands, not measured chamber pressure. Invalid samples masked; theta wraps have gaps.\n');
for name={'lambda_gains','damping_gains'}
 x=align(ds,name{1},t);ii=[1;find(any(diff(x)~=0,2))+1];fprintf(fid,'%s changes [time values]:\n',name{1});fprintf(fid,'%.3f %.6g %.6g\n',[t(ii) x(ii,:)]');
end
good=sel&ov&obs;eth=atan2(sin(td-tm),cos(td-tm));
fprintf(fid,'Theta MAE %.6f rad; phi MAE %.6f rad\n',mean(abs(eth(good)),'omitnan'),mean(abs(pd(good)-pm(good)),'omitnan'));
fprintf(fid,'XYZ RMSE mm: %s\n',mat2str(1000*sqrt(mean((xyz(sel,:)-xyzDesired(sel,:)).^2,1,'omitnan')),6));
fprintf(fid,'Length RMSE mm: %s\n',mat2str(1000*sqrt(mean((lengths(sel,:)-desiredLength(sel,:)).^2,1,'omitnan')),6));
fclose(fid);
end
function pair(t,d,m,c)
plot(t,d,'k--','LineWidth',1.8);hold on;plot(t,m,'Color',c,'LineWidth',1.6);
end
function y=wrapped(x)
y=atan2(sin(x),cos(x));y([false;abs(diff(y))>pi])=NaN;
end
function x=align(ds,name,t)
[tt,x]=read(ds,name);if ~isequal(t,tt),x=interp1(tt,x,t,'previous',NaN);end
end
function [t,x]=read(ds,name)
ix=find(strcmp(ds.getElementNames,name),1);assert(~isempty(ix),'Missing %s',name);
v=ds{ix}.Values;t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
