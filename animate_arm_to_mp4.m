function animate_arm_to_mp4(t, X, params, filename, varargin)
% animate_arm_to_mp4
% ------------------------------------------------------------
% Draws and records a continuum arm simulation as an MP4.
%
% Inputs:
%   t        : [Nt x 1] time vector from ODE solver
%   X        : [Nt x 4N] state matrix [l; dl]
%   params   : struct with fields
%              .N, .L, .r, .cog_xi, .mi
%   filename : string, e.g. 'arm_sim.mp4'
%
% Optional (name-value):
%   'dt'     : animation timestep (default = 0.01)
%   'fps'    : video FPS (default = 30)
%   'flipZ'  : true/false (default = true)
%
% ------------------------------------------------------------

%% ---------------- Parse inputs ----------------
p = inputParser;
addParameter(p,'dt',0.01);
addParameter(p,'fps',30);
addParameter(p,'flipZ',true);
parse(p,varargin{:});

dt_anim = p.Results.dt;
fps     = p.Results.fps;
useFlipZ = p.Results.flipZ;

N = params.N;
L = params.L;
r = params.r;

%% ---------------- Resample solution ----------------
t_anim = (t(1):dt_anim:t(end)).';
X_anim = interp1(t, X, t_anim, 'pchip');
Nt = numel(t_anim);

xi = linspace(0,1,25);
Nx = numel(xi);

%% ---------------- Video writer ----------------
vw = VideoWriter(filename,'MPEG-4');
vw.FrameRate = fps;
vw.Quality   = 100;
open(vw);

%% ---------------- Figure setup ----------------
fig = figure('Color','w','Position',[100 100 1200 500]);

% ===== Arm =====
ax1 = subplot(1,2,1);
hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');
xlabel(ax1,'X'); ylabel(ax1,'Y'); zlabel(ax1,'Z');
view(ax1,[11 13]);
rotate3d(ax1,'on');
xlim(ax1,[-1 1]); ylim(ax1,[-1 1]); zlim(ax1,[-1.5 1.1]);

h_seg = gobjects(N,1);
h_tip_seg = gobjects(N,1);

for s = 1:N
    h_seg(s) = plot3(ax1,NaN,NaN,NaN,'LineWidth',2);
    col = h_seg(s).Color;
    h_tip_seg(s) = plot3(ax1,NaN,NaN,NaN,'s', ...
        'MarkerFaceColor',col,'MarkerEdgeColor','none');
end

h_tip = plot3(ax1,NaN,NaN,NaN,'o', ...
    'MarkerFaceColor',[.8 .2 .2],'MarkerEdgeColor','none');

title(ax1,'Continuum arm');

% ===== Length plot =====
ax2 = subplot(1,2,2);
hold(ax2,'on'); grid(ax2,'on');
xlabel(ax2,'time (s)');
ylabel(ax2,'length change (m)');
title(ax2,'Length states');

h_l = gobjects(2*N,1);
for i = 1:2*N
    h_l(i) = plot(ax2,t_anim,X_anim(:,i),'LineWidth',1.2);
end

yl = ylim(ax2);
h_cursor = plot(ax2,[t_anim(1) t_anim(1)],yl,'k--','LineWidth',1.2);

h_dot = gobjects(2*N,1);
for i = 1:2*N
    h_dot(i) = plot(ax2,t_anim(1),X_anim(1,i),'o',...
        'MarkerFaceColor',h_l(i).Color,'MarkerEdgeColor','none');
end

%% ---------------- Animation loop ----------------
Rflip = eye(3);
if useFlipZ
    Rflip = [1 0 0; 0 -1 0; 0 0 -1];
end

for k = 1:Nt

    tk = t_anim(k);
    l_mat = reshape(X_anim(k,1:2*N).',[2,N]).';

    % ----- update length plot -----
    set(h_cursor,'XData',[tk tk],'YData',ylim(ax2));
    for i = 1:2*N
        set(h_dot(i),'XData',tk,'YData',X_anim(k,i));
    end

    % ----- draw arm -----
    Rglob = eye(3);
    Pglob = zeros(3,1);

    for s = 1:N
        length_s = [0, l_mat(s,1), l_mat(s,2)];

        Xs = zeros(Nx,1); Ys = zeros(Nx,1); Zs = zeros(Nx,1);
        for j = 1:Nx
            [~,~,Ploc] = HTM_nume(length_s,xi(j),L,r);
            P = Pglob + Rglob*Ploc;
            P = Rflip*P;
            Xs(j)=P(1); Ys(j)=P(2); Zs(j)=P(3);
        end

        set(h_seg(s),'XData',Xs,'YData',Ys,'ZData',Zs);

        [~,Rtip,Ptip] = HTM_nume(length_s,1,L,r);
        Ptip_g = Rflip*(Pglob + Rglob*Ptip);
        set(h_tip_seg(s),'XData',Ptip_g(1),'YData',Ptip_g(2),'ZData',Ptip_g(3));

        Pglob = Pglob + Rglob*Ptip;
        Rglob = Rglob*Rtip;
    end

    P_end = Rflip*Pglob;
    set(h_tip,'XData',P_end(1),'YData',P_end(2),'ZData',P_end(3));

    drawnow limitrate

    % ----- write frame -----
    frame = getframe(fig);
    writeVideo(vw,frame);
end

close(vw);
fprintf('Saved video: %s\n', filename);

end
