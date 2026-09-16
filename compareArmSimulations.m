function compareArmSimulations(dataFolder, dt, playbackSpeed)
%COMPAREARMSIMULATIONS Three synchronized simulations in one window.
% Run from the folder containing the MAT files and HTM_nume:
%   compareArmSimulations(pwd)
% Optional: compareArmSimulations(pwd, 1/60, 0.5) for half-speed playback.
% Uses each recording's own params. Holds the last pose after a run ends.
if nargin < 1 || isempty(dataFolder), dataFolder = pwd; end
if nargin < 2 || isempty(dt), dt = 1/60; end
if nargin < 3 || isempty(playbackSpeed), playbackSpeed = 1; end
validateattributes(dt, {'numeric'}, {'scalar','positive','finite'});
validateattributes(playbackSpeed, {'numeric'}, {'scalar','positive','finite'});
files = {'original_sim_data.mat', 'comparison_after_mass_scaling.mat', ...
         'comparison_after_C_correction.mat'};
names = {'Original', 'Mass scaling corrected', 'C matrix corrected'};
runs = cell(1,3);
tStart = inf; tEnd = -inf; qMin = inf; qMax = -inf;
for a = 1:3
    d = load(fullfile(dataFolder,files{a}), 't','X','params');
    assert(all(isfield(d,{'t','X','params'})), 'Missing t, X or params in %s.',files{a});
    d.t = d.t(:);
    assert(numel(d.t) >= 2 && all(diff(d.t)>0), 'Time must increase in %s.',files{a});
    assert(size(d.X,1)==numel(d.t) && size(d.X,2)>=2*d.params.N, ...
        'Unexpected state dimensions in %s.',files{a});
    d.X = d.X(:,1:2*d.params.N);
    assert(all(isfinite(d.X(:))) && all(isfinite(d.t)), 'Nonfinite data in %s.',files{a});
    runs{a} = d;
    tStart = min(tStart,d.t(1)); tEnd = max(tEnd,d.t(end));
    qMin = min(qMin,min(d.X(:))); qMax = max(qMax,max(d.X(:)));
end
tAnim = (tStart:dt:tEnd).';
if tAnim(end)<tEnd, tAnim(end+1)=tEnd; end
qPad = max(0.05*(qMax-qMin),1e-6);
qLimits = [qMin-qPad,qMax+qPad];
fig = figure('Name','Arm simulation comparison','NumberTitle','off', ...
    'Color','w','Position',[80 80 1400 800]);
layout = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
heading = title(layout,'Simulation comparison');
armAx = gobjects(1,3); curveAx = gobjects(1,3);
armTitle = gobjects(1,3); cursor = gobjects(1,3);
segments = cell(1,3); tips = cell(1,3); dots = cell(1,3);
xi = linspace(0,1,25);
flipZ = diag([1,-1,-1]); % Same downward display orientation as drawingArms.
for a = 1:3
    d = runs{a}; N = d.params.N;
    d.sampleTime = min(max(tAnim,d.t(1)),d.t(end));
    d.qAnim = interp1(d.t,d.X,d.sampleTime,'pchip');
    runs{a} = d;
    armAx(a) = nexttile(layout,a);
    ax = armAx(a); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    view(ax,[11,13]); xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');
    % Identical scales, matching the existing drawingArms view.
    xlim(ax,[-1 1]); ylim(ax,[-1 1]); zlim(ax,[-1.5 1.1]);
    armTitle(a) = title(ax,names{a});
    colors = lines(N);
    segments{a} = gobjects(N,1); tips{a} = gobjects(N,1);
    for s = 1:N
        segments{a}(s) = plot3(ax,nan,nan,nan,'LineWidth',2,'Color',colors(s,:));
        tips{a}(s) = plot3(ax,nan,nan,nan,'s','Color',colors(s,:), ...
            'MarkerFaceColor',colors(s,:),'MarkerSize',6);
    end
    curveAx(a) = nexttile(layout,a+3);
    ax = curveAx(a); hold(ax,'on'); grid(ax,'on');
    colors = lines(2*N); curves = gobjects(2*N,1); dots{a} = gobjects(2*N,1);
    labels = cell(1,2*N);
    for j = 1:2*N
        curves(j) = plot(ax,d.t,d.X(:,j),'Color',colors(j,:),'LineWidth',1);
        labels{j} = sprintf('q_%d',j);
        dots{a}(j) = plot(ax,nan,nan,'o','Color',colors(j,:), ...
            'MarkerFaceColor',colors(j,:),'HandleVisibility','off');
    end
    xlim(ax,[tStart,tEnd]); ylim(ax,qLimits);
    xlabel(ax,'Time (s)'); ylabel(ax,'Length change (m)');
    legend(ax,curves,labels,'Location','best','NumColumns',2);
    cursor(a) = plot(ax,[tStart tStart],qLimits,'k--','HandleVisibility','off');
end
linkaxes(curveAx,'xy');
clock = tic;
for k = 1:numel(tAnim)
    if ~isgraphics(fig), return; end
    tk = tAnim(k);
    for a = 1:3
        d = runs{a}; N = d.params.N;
        q = reshape(d.qAnim(k,:),2,N).';
        R = eye(3); P = zeros(3,1);
        for s = 1:N
            l = [0,q(s,:)];
            points = zeros(3,numel(xi));
            for j = 1:numel(xi)
                [~,~,p] = HTM_nume(l,xi(j),d.params.L,d.params.r);
                points(:,j) = flipZ*(P+R*p);
            end
            set(segments{a}(s),'XData',points(1,:),'YData',points(2,:),'ZData',points(3,:));
            set(tips{a}(s),'XData',points(1,end),'YData',points(2,end),'ZData',points(3,end));
            [~,Rt,pt] = HTM_nume(l,1,d.params.L,d.params.r);
            P = P+R*pt; R = R*Rt;
        end
        state = '';
        if tk > d.t(end)
            state = sprintf(' | ended at %.2f s; held',d.t(end));
        elseif tk < d.t(1)
            state = ' | not started; held';
        end
        set(armTitle(a),'String',[names{a},state]);
        ts = d.sampleTime(k);
        set(cursor(a),'XData',[ts ts]);
        for j = 1:2*N
            set(dots{a}(j),'XData',ts,'YData',d.qAnim(k,j));
        end
    end
    set(heading,'String',sprintf('Simulation comparison | t = %.2f s',tk));
    drawnow;
    if ~isgraphics(fig), return; end
    if k < numel(tAnim)
        waitTime = (tAnim(k+1)-tStart)/playbackSpeed-toc(clock);
        if waitTime > 0, pause(waitTime); end
    end
end
end
