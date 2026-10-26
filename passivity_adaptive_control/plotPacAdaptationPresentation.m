function figures=plotPacAdaptationPresentation(filePath,outputDir)
%PLOTPACADAPTATIONPRESENTATION Saved parameters/gains; no model/hardware runs.
if ~isfolder(outputDir),mkdir(outputDir);end
S=load(filePath,'out');ds=S.out.logsout;[t,th]=read(ds,'parameterEstimate');
assert(size(th,2)==10,'Expected ten adaptive estimates.');
[tl,lambda]=read(ds,'lambda_gains');[tk,Ks]=read(ds,'damping_gains');
il=find(any(diff(lambda)~=0,2))+1;ik=find(any(diff(Ks)~=0,2))+1;
events=unique([tl(il);tk(ik)]);colors=[.05 .37 .70;.85 .33 .10;.49 .18 .56];
figures=gobjects(3,1);stems={'fullRun_parameter_adaptations','fullRun_feedback_gains','fullRun_adaptation_switches'};
for f=1:3
 figures(f)=figure('Visible','off','Color','w','Position',[80 80 1400 1050]);
 if f==1
  tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
  groups={1:3,4:5,6:8,9:10};labels={'Stiffness (N/m)','Damping (N s/m)','History weights (N)','Bias (N)'};
  legends={{'k_1','k_2','k_3'},{'d_2','d_3'},{'h_1','h_2','h_3'},{'b_2','b_3'}};
  for j=1:4
   ax=nexttile;colororder(ax,colors);plot(t,th(:,groups{j}),'LineWidth',1.6);ylabel(labels{j});
   legend(legends{j},'Location','northeast','Orientation','horizontal');
   if j==4,xlabel('Time (s)');end
  end
 elseif f==2
  tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
  ax=nexttile;colororder(ax,colors);stairs(tl,lambda,'LineWidth',2);ylabel('\lambda (s^{-1})');legend('\lambda_2','\lambda_3','Location','northeast','Orientation','horizontal');
  ax=nexttile;colororder(ax,colors);stairs(tk,Ks,'LineWidth',2);ylabel('K_s (N s/m)');xlabel('Time (s)');legend('K_{s,2}','K_{s,3}','Location','northeast','Orientation','horizontal');
 else
  names={'enableAdaptation','enableCommonPressureAdaptation','adaptationActive'};
  labels={'Adaptation enabled','Common-P enabled','Update active'};
  if any(strcmp(ds.getElementNames,'freezeHysteresisWeights'))
   names{end+1}='freezeHysteresisWeights';labels{end+1}='History weights frozen';
  end
  tiledlayout(numel(names),1,'TileSpacing','compact','Padding','compact');
  for j=1:numel(names)
   nexttile;[tt,x]=read(ds,names{j});
   % Constant sparse logs persist through the recording, not just last event.
   if tt(end)<t(end),tt(end+1)=t(end);x(end+1,:)=x(end,:);end
   stairs(tt,x,'Color',colors(1,:),'LineWidth',1.5);ylabel(labels{j});ylim([-.15 1.3]);yticks([0 1]);
   if j==numel(names),xlabel('Time (s)');end
  end
 end
 for ax=findall(figures(f),'Type','axes')'
  set(ax,'FontName','Arial','FontSize',17,'LineWidth',1,'Box','on','XGrid','on','YGrid','on','GridAlpha',.13,'TickDir','out');xlim(ax,[t(1) t(end)]);
  if f~=3,b=ylim(ax);pad=max(.1,.18*diff(b));ylim(ax,b+[-pad pad]);end
  for event=events'
   xline(ax,event,'k:','LineWidth',1.2,'HandleVisibility','off');
  end
 end
 drawnow;base=fullfile(outputDir,stems{f});
 exportgraphics(figures(f),[base '.png'],'Resolution',300);
 exportgraphics(figures(f),[base '.pdf'],'ContentType','vector');
 print(figures(f),[base '.svg'],'-dsvg');savefig(figures(f),[base '.fig']);
 fprintf('EXPORTED %s\n',base);
end
fid=fopen(fullfile(outputDir,'adaptation_gain_notes.txt'),'w');assert(fid>=0);
fprintf(fid,'Source: %s\nParameter order: k1 k2 k3 d2 d3 h1 h2 h3 b2 b3\n',filePath);
fprintf(fid,'Plots use retained parameterEstimate, not force-evaluation matrix outputs that can zero on invalid data. Estimates are compensation parameters, not unique physical identification.\n');
fprintf(fid,'Dotted lines mark feedback gain changes at %s s.\n',mat2str(events'));
fprintf(fid,'lambda [time l2 l3]:\n');ii=[1;il];fprintf(fid,'%.3f %.6g %.6g\n',[tl(ii) lambda(ii,:)]');
fprintf(fid,'Ks [time Ks2 Ks3]:\n');ii=[1;ik];fprintf(fid,'%.3f %.6g %.6g\n',[tk(ii) Ks(ii,:)]');
fprintf(fid,'Initial parameters %s\nFinal parameters %s\n',mat2str(th(1,:),7),mat2str(th(end,:),7));
names=ds.getElementNames;
for name={'adaptationGamma','parameterScale','parameterMaxRate','parameterLower','parameterUpper','parameterInitial','adaptationOptions','historyOptions','pressureOptions','supervisorOptions'}
 if ~any(strcmp(names,name{1})),fprintf(fid,'Not logged: %s\n',name{1});end
end
fclose(fid);
end
function [t,x]=read(ds,name)
ix=find(strcmp(ds.getElementNames,name),1);assert(~isempty(ix),'Missing %s',name);
v=ds{ix}.Values;t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
