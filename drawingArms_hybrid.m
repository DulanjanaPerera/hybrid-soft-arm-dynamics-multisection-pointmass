function [figArm,figStates]=drawingArms_hybrid(t,X,frameInterval,params,videoPath)
% Animate the distributed arm; optionally record Figure 1 as an MP4.
% Purple diamonds mark added point masses. Empty videoPath disables recording.
if nargin<5, videoPath=''; end
recordVideo=~isempty(videoPath);
assert(numel(t)==size(X,1) && size(X,2)>=6 && params.N==3);
assert(isscalar(frameInterval) && frameInterval>0);
if recordVideo
    assert(ischar(videoPath) || (isstring(videoPath) && isscalar(videoPath)), ...
        'videoPath must be a filename.');
    videoPath=char(videoPath);
    assert(~exist(videoPath,'file'), 'Video already exists: %s',videoPath);
    video=VideoWriter(videoPath,'MPEG-4');
    video.FrameRate=1/frameInterval;
    video.Quality=95;
    open(video);
    videoCleanup=onCleanup(@() close(video)); % Finalize MP4 even if drawing fails.
end
flip=diag([1,-1,-1]); % Same display flip as the original drawingArms.m.
xiBackbone=linspace(0,1,25); N=params.N;
figArm=figure(1); clf(figArm);
if recordVideo, figArm.Position=[100 100 960 720]; end
ax=axes(figArm); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');
view(ax,[11,13]); rotate3d(figArm,'on');
[low,high]=animationBounds(X,params,flip);
margin=max(.06,.08*max(high-low));
xlim(ax,[low(1)-margin high(1)+margin]);
ylim(ax,[low(2)-margin high(2)+margin]);
zlim(ax,[low(3)-margin high(3)+margin]);
segment=gobjects(N,1); tipMarker=gobjects(N,1);
payloadMarker=gobjects(N,1);
for n=1:N
    segment(n)=plot3(ax,nan,nan,nan,'LineWidth',2);
    col=segment(n).Color;
    tipMarker(n)=plot3(ax,nan,nan,nan,'s','MarkerSize',7, ...
        'MarkerFaceColor',col,'MarkerEdgeColor','none');
    if params.addedMass(n)>0
        payloadMarker(n)=plot3(ax,nan,nan,nan,'d','MarkerSize',11, ...
            'MarkerFaceColor',[.72 .1 .65],'MarkerEdgeColor','k', ...
            'DisplayName',sprintf('Added mass %d: %.1f g at xi=%.2f', ...
            n,1000*params.addedMass(n),params.addedXi(n)));
    end
end
finalTip=plot3(ax,nan,nan,nan,'o','MarkerSize',8, ...
    'MarkerFaceColor',[.8 .2 .2],'MarkerEdgeColor','none');
active=payloadMarker(isgraphics(payloadMarker));
if ~isempty(active), legend(ax,active,'Location','best'); end
titleHandle=title(ax,'');

figStates=figure(2); clf(figStates);
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
title(axQ,'Hybrid arm coordinates');
cursor=plot(axQ,[t(1) t(1)],ylim(axQ),'k--','LineWidth',1.2);
start=tic;
for k=1:numel(t)
    set(cursor,'XData',[t(k) t(k)],'YData',ylim(axQ));
    for i=1:6, set(qDots(i),'XData',t(k),'YData',X(k,i)); end
    q=reshape(X(k,1:6),2,N).';
    R=eye(3); origin=zeros(3,1);
    for n=1:N
        l=[0,q(n,:)]; xyz=zeros(3,numel(xiBackbone));
        for j=1:numel(xiBackbone)
            [~,~,pLocal]=HTM_nume(l,xiBackbone(j),params.L,params.r);
            xyz(:,j)=flip*(origin+R*pLocal);
        end
        set(segment(n),'XData',xyz(1,:),'YData',xyz(2,:),'ZData',xyz(3,:));
        if params.addedMass(n)>0
            [~,~,pLocal]=HTM_nume(l,params.addedXi(n),params.L,params.r);
            p=flip*(origin+R*pLocal);
            set(payloadMarker(n),'XData',p(1),'YData',p(2),'ZData',p(3));
        end
        [~,Rtip,pTip]=HTM_nume(l,1,params.L,params.r);
        nextOrigin=origin+R*pTip; tip=flip*nextOrigin;
        set(tipMarker(n),'XData',tip(1),'YData',tip(2),'ZData',tip(3));
        origin=nextOrigin; R=R*Rtip;
    end
    tip=flip*origin;
    set(finalTip,'XData',tip(1),'YData',tip(2),'ZData',tip(3));
    infT = params.times/length(t);
    Frq = 1/infT;
    if params.exampleOnly, prefix='Illustrative distributed arm + added mass';
    else, prefix='Distributed arm + added point mass'; end
    set(titleHandle,'String',sprintf('%s | frame %d/%d | t=%.2f s \nInference time %.5f s  |  Frequency %.3f Hz', ...
        prefix,k,numel(t),t(k), infT, Frq));
    if recordVideo
        drawnow; % Render every frame before capture; limitrate can skip frames.
        writeVideo(video,getframe(figArm));
    else
        drawnow limitrate;
    end
    remaining=(k-1)*frameInterval-toc(start);
    if remaining>0, pause(remaining); end
end
drawnow;
if recordVideo, clear videoCleanup; end
end

function [low,high]=animationBounds(X,params,flip)
% Stable axes that contain all three backbones over the whole animation.
low=inf(3,1); high=-inf(3,1);
for k=1:size(X,1)
    q=reshape(X(k,1:6),2,3).';
    R=eye(3); origin=zeros(3,1);
    for n=1:3
        l=[0,q(n,:)];
        for xi=linspace(0,1,9)
            [~,~,p]=HTM_nume(l,xi,params.L,params.r);
            x=flip*(origin+R*p);
            low=min(low,x); high=max(high,x);
        end
        [~,Rt,pt]=HTM_nume(l,1,params.L,params.r);
        origin=origin+R*pt; R=R*Rt;
    end
end
end
