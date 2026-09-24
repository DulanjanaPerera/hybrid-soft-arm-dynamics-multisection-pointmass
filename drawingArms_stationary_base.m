function [figArm,figStates] = drawingArms_stationary_base(t,X,frameInterval,params,videoPath)
% Animate the standard distributed arm at one fixed base orientation.
% Geometry is shown in world coordinates, without the hybrid display flip.
if nargin<5, videoPath=''; end
recordVideo=~isempty(videoPath);
assert(numel(t)==size(X,1) && size(X,2)>=6 && params.N==3);
assert(isscalar(frameInterval) && frameInterval>0);
Rbase=params.R_world_from_arm;
assert(isequal(size(Rbase),[3,3]) && all(isfinite(Rbase(:))) ...
    && norm(Rbase.'*Rbase-eye(3),'fro')<1e-10 && det(Rbase)>0, ...
    'R_world_from_arm must be a proper rotation.');

if recordVideo
    assert(ischar(videoPath) || (isstring(videoPath) && isscalar(videoPath)), ...
        'videoPath must be a filename.');
    videoPath=char(videoPath);
    assert(~exist(videoPath,'file'),'Video already exists: %s',videoPath);
    video=VideoWriter(videoPath,'MPEG-4');
    video.FrameRate=1/frameInterval;
    video.Quality=95;
    open(video);
    videoCleanup=onCleanup(@() close(video)); %#ok<NASGU>
end

xiBackbone=linspace(0,1,25);
figArm=figure(1); clf(figArm);
figArm.Name='Stationary-base standard arm animation';
if recordVideo, figArm.Position=[100 100 960 720]; end
ax=axes(figArm); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'World X (m)'); ylabel(ax,'World Y (m)');
zlabel(ax,'World Z (m)');
view(ax,[11,13]); rotate3d(figArm,'on');
xlim(ax,[-1 1]); ylim(ax,[-1 1]); zlim(ax,[-1.5 1.1]);
triadLength=.05;
triadColors=[1 0 0;0 .6 0;0 .2 1];
for a=1:3
    v=triadLength*Rbase(:,a);
    plot3(ax,[0 v(1)],[0 v(2)],[0 v(3)], ...
        'Color',triadColors(a,:),'LineWidth',2, ...
        'DisplayName',sprintf('Arm axis %d',a));
end
plot3(ax,0,0,0,'ko','MarkerFaceColor','k','DisplayName','Base origin');
segment=gobjects(3,1); tipMarker=gobjects(3,1);
for n=1:3
    segment(n)=plot3(ax,nan,nan,nan,'LineWidth',2, ...
        'DisplayName',sprintf('Section %d',n));
    col=segment(n).Color;
    tipMarker(n)=plot3(ax,nan,nan,nan,'s','MarkerSize',7, ...
        'MarkerFaceColor',col,'MarkerEdgeColor','none', ...
        'HandleVisibility','off');
end
finalTip=plot3(ax,nan,nan,nan,'o','MarkerSize',8, ...
    'MarkerFaceColor',[.8 .2 .2],'MarkerEdgeColor','none', ...
    'HandleVisibility','off');
legend(ax,'Location','best');
titleHandle=title(ax,'');

figStates=figure(2); clf(figStates);
figStates.Name='Stationary-base standard arm coordinates';
axQ=axes(figStates); hold(axQ,'on'); grid(axQ,'on');
qLines=gobjects(6,1); qDots=gobjects(6,1);
for i=1:6
    qLines(i)=plot(axQ,t,X(:,i),'LineWidth',1.4);
    qDots(i)=plot(axQ,t(1),X(1,i),'o','MarkerSize',6, ...
        'MarkerFaceColor',qLines(i).Color,'MarkerEdgeColor','none');
end
legend(axQ,qLines,{'l_{12}','l_{13}','l_{22}','l_{23}', ...
    'l_{32}','l_{33}'},'Location','best');
xlabel(axQ,'Time (s)'); ylabel(axQ,'Length change (m)');
title(axQ,'Standard arm coordinates');
cursor=plot(axQ,[t(1) t(1)],ylim(axQ),'k--','LineWidth',1.2);

playClock=tic;
for k=1:numel(t)
    set(cursor,'XData',[t(k) t(k)],'YData',ylim(axQ));
    for i=1:6, set(qDots(i),'XData',t(k),'YData',X(k,i)); end
    q=reshape(X(k,1:6),2,3).';
    Rsection=eye(3); origin=zeros(3,1);
    for n=1:3
        l=[0,q(n,:)];
        xyz=zeros(3,numel(xiBackbone));
        for j=1:numel(xiBackbone)
            [~,~,pLocal]=HTM_nume(l,xiBackbone(j),params.L,params.r);
            xyz(:,j)=Rbase*(origin+Rsection*pLocal);
        end
        set(segment(n),'XData',xyz(1,:),'YData',xyz(2,:),'ZData',xyz(3,:));
        [~,Rtip,pTip]=HTM_nume(l,1,params.L,params.r);
        origin=origin+Rsection*pTip;
        Rsection=Rsection*Rtip;
        tip=Rbase*origin;
        set(tipMarker(n),'XData',tip(1),'YData',tip(2),'ZData',tip(3));
    end
    tip=Rbase*origin;
    set(finalTip,'XData',tip(1),'YData',tip(2),'ZData',tip(3));
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
