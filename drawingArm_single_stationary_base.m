function [figArm,figStates] = drawingArm_single_stationary_base(t,X,frameInterval,params,videoPath)
% Animate one distributed section at a fixed world-frame base orientation.
if nargin<5, videoPath=''; end
recordVideo=~isempty(videoPath);
assert(params.N==1 && numel(t)==size(X,1) && size(X,2)>=2);
assert(isscalar(frameInterval) && frameInterval>0);
Rbase=params.R_world_from_arm;
assert(isequal(size(Rbase),[3,3]) && all(isfinite(Rbase(:))) ...
    && norm(Rbase.'*Rbase-eye(3),'fro')<1e-10 && det(Rbase)>0);

if recordVideo
    videoPath=char(videoPath);
    assert(~isfile(videoPath),'Video already exists: %s',videoPath);
    video=VideoWriter(videoPath,'MPEG-4');
    video.FrameRate=1/frameInterval;
    video.Quality=95;
    open(video);
    videoCleanup=onCleanup(@() close(video)); %#ok<NASGU>
end

figArm=figure(1); clf(figArm);
figArm.Name='Single-module stationary-base animation';
if recordVideo, figArm.Position=[100 100 960 720]; end
ax=axes(figArm); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'World X (m)'); ylabel(ax,'World Y (m)');
zlabel(ax,'World Z (m)');
view(ax,[11,13]); rotate3d(figArm,'on');
plotRange=max(0.22,1.4*params.L);
xlim(ax,[-plotRange plotRange]);
ylim(ax,[-plotRange plotRange]);
zlim(ax,[-plotRange plotRange]);
triadColors=[1 0 0;0 .6 0;0 .2 1];
for a=1:3
    v=.05*Rbase(:,a);
    plot3(ax,[0 v(1)],[0 v(2)],[0 v(3)], ...
        'Color',triadColors(a,:),'LineWidth',2, ...
        'DisplayName',sprintf('Arm axis %d',a));
end
plot3(ax,0,0,0,'ko','MarkerFaceColor','k', ...
    'DisplayName','Base origin');
segment=plot3(ax,nan,nan,nan,'LineWidth',2, ...
    'DisplayName','Distributed section');
tipMarker=plot3(ax,nan,nan,nan,'o','MarkerSize',8, ...
    'MarkerFaceColor',[.8 .2 .2],'MarkerEdgeColor','none', ...
    'HandleVisibility','off');
legend(ax,'Location','best');
titleHandle=title(ax,'');

figStates=figure(2); clf(figStates);
figStates.Name='Single-module coordinates';
axQ=axes(figStates); hold(axQ,'on'); grid(axQ,'on');
qLines=gobjects(2,1); qDots=gobjects(2,1);
for i=1:2
    qLines(i)=plot(axQ,t,X(:,i),'LineWidth',1.4);
    qDots(i)=plot(axQ,t(1),X(1,i),'o','MarkerSize',6, ...
        'MarkerFaceColor',qLines(i).Color,'MarkerEdgeColor','none');
end
legend(axQ,qLines,{'l_{12}','l_{13}'},'Location','best');
xlabel(axQ,'Time (s)'); ylabel(axQ,'Length change (m)');
title(axQ,'Single-module coordinates');
cursor=plot(axQ,[t(1) t(1)],ylim(axQ),'k--','LineWidth',1.2);

xiBackbone=linspace(0,1,25);
playClock=tic;
for k=1:numel(t)
    set(cursor,'XData',[t(k) t(k)],'YData',ylim(axQ));
    for i=1:2, set(qDots(i),'XData',t(k),'YData',X(k,i)); end
    l=[0,X(k,1:2)];
    xyz=zeros(3,numel(xiBackbone));
    for j=1:numel(xiBackbone)
        [~,~,pLocal]=HTM_nume(l,xiBackbone(j),params.L,params.r);
        xyz(:,j)=Rbase*pLocal;
    end
    set(segment,'XData',xyz(1,:),'YData',xyz(2,:),'ZData',xyz(3,:));
    tip=xyz(:,end);
    set(tipMarker,'XData',tip(1),'YData',tip(2),'ZData',tip(3));
    set(titleHandle,'String',sprintf('%s | frame %d/%d | t=%.3f s', ...
        params.orientationLabel,k,numel(t),t(k)));
    if recordVideo
        drawnow;
        writeVideo(video,getframe(figArm));
    else
        drawnow limitrate;
    end
    remaining=(k-1)*frameInterval-toc(playClock);
    if remaining>0, pause(remaining); end
end
drawnow;
if recordVideo, clear videoCleanup; end
end
