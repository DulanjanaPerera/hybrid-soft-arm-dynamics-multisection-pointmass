function drawingArms(t, X, dt, params)

loc = params.loc;
m = params.m;
stiff = params.stiff;
damp = params.damp;
times = params.times;

%% ===================== RESAMPLE FOR SMOOTH ANIMATION =====================
dt_anim = dt;                        % desired animation timestep
t_anim = (t(1):dt_anim:t(end)).';      % uniform time grid

% Interpolate states onto uniform grid
% - 'pchip' preserves shape better than linear and avoids overshoot vs spline
X_anim = interp1(t, X, t_anim, 'pchip');

Nt = numel(t_anim);

%% ===================== DRAW + PLOT STATES (SIMULTANEOUS ANIMATION) =====================
N = params.N;
L = params.L;
r = params.r;

xi = linspace(0,1,25);
Nx = length(xi);

useFlipZ = true;  % set true if you want 180 deg flip about X

%% --------------------- Figure 1: Arm drawing ---------------------
figure(1); clf
ax1 = axes; hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');
xlabel(ax1,'X'); ylabel(ax1,'Y'); zlabel(ax1,'Z');
rotate3d(ax1,'on');
view([11,13])
str1 = sprintf('%d-section continuum arm with CoG at %.2f', N, loc);
ht1 = title(ax1,str1);

xlim(ax1, [-1 1]); ylim(ax1, [-1 1]); zlim(ax1, [-1.5 1.1]);

h_seg = gobjects(N,1);
h_tip_seg = gobjects(N,1);
segColor = zeros(N,3);

for s = 1:N
    h_seg(s) = plot3(ax1, NaN, NaN, NaN, 'LineWidth', 2);
    segColor(s,:) = h_seg(s).Color;
    h_tip_seg(s) = plot3(ax1, NaN, NaN, NaN, 's', ...
        'MarkerFaceColor', segColor(s,:), 'MarkerEdgeColor','none', 'MarkerSize', 7);
end

h_tip = plot3(ax1, NaN, NaN, NaN, 'o', ...
    'MarkerFaceColor',[.8 .2 .2], 'MarkerEdgeColor','none', 'MarkerSize', 7);

%% --------------------- Figure 2: Length-change plot ---------------------
figure(2); clf
ax2 = axes; hold(ax2,'on'); grid(ax2,'on');
axis(ax2,'tight');

% Plot original (optional) OR interpolated (recommended)
% If you want the plotted curves to match the animation points, plot t_anim/X_anim.
h_l = gobjects(2*N,1);
for i = 1:2*N
    h_l(i) = plot(ax2, t_anim, X_anim(:,i), 'LineWidth', 1.5);
end

str = sprintf('length change of %d Sections', N);
title(ax2, str);
xlabel(ax2, 'time (s)');
ylabel(ax2, 'length change (m)');

leg = strings(2*N,1);
for ksec = 1:N
    leg(2*(ksec-1)+1) = "l_{" + ksec + "1}";
    leg(2*(ksec-1)+2) = "l_{" + ksec + "2}";
end
legend(ax2, leg, 'Location','best');

yl = ylim(ax2);
h_cursor = plot(ax2, [t_anim(1) t_anim(1)], yl, 'k--', 'LineWidth', 1.2);

h_dot = gobjects(2*N,1);
for i = 1:2*N
    col = h_l(i).Color;
    h_dot(i) = plot(ax2, t_anim(1), X_anim(1,i), 'o', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
end

%% --------------------- Animation loop (uniform timestep) ---------------------
% Optional: try to keep "real time" pacing (best-effort)
t_start_wall = tic;

for k = 1:Nt

    tk = t_anim(k);

    % --- best-effort real-time pacing ---
    % comment this out if you just want fastest possible animation.
    target_elapsed = tk - t_anim(1);
    while toc(t_start_wall) < target_elapsed
        % busy-wait very lightly (no pause jitter); you can replace with pause(0)
    end

    % ===== Update Figure 2 =====
    set(h_cursor, 'XData', [tk tk], 'YData', ylim(ax2));
    for i = 1:2*N
        set(h_dot(i), 'XData', tk, 'YData', X_anim(k,i));
    end

    % ===== Update Figure 1 =====
    l_mat = reshape(X_anim(k, 1:2*N).', [2, N]).';  % N x 2

    Rglob = eye(3);
    Pglob = zeros(3,1);
    final_tip = [NaN;NaN;NaN];

    for s = 1:N
        length_s = [0, l_mat(s,1), l_mat(s,2)];

        Xs = zeros(Nx,1); Ys = zeros(Nx,1); Zs = zeros(Nx,1);
        for j = 1:Nx
            [~, ~, Ploc] = HTM_nume(length_s, xi(j), L, r);
            Ppoint = Pglob + Rglob * Ploc;

            if useFlipZ
                Ppoint = [1 0 0; 0 -1 0; 0 0 -1] * Ppoint;
            end

            Xs(j) = Ppoint(1); Ys(j) = Ppoint(2); Zs(j) = Ppoint(3);
        end

        set(h_seg(s), 'XData', Xs, 'YData', Ys, 'ZData', Zs);

        [~, Rtip, Ptip] = HTM_nume(length_s, 1, L, r);
        Pglob_tip = Pglob + Rglob * Ptip;

        if useFlipZ
            Pglob_tip = [1 0 0; 0 -1 0; 0 0 -1] * Pglob_tip;
        end

        set(h_tip_seg(s), 'XData', Pglob_tip(1), 'YData', Pglob_tip(2), 'ZData', Pglob_tip(3));

        % advance (unflipped)
        Pglob = Pglob + Rglob * Ptip;
        Rglob = Rglob * Rtip;

        final_tip = Pglob;
        if useFlipZ
            final_tip = [1 0 0; 0 -1 0; 0 0 -1] * final_tip;
        end
    end

    set(h_tip, 'XData', final_tip(1), 'YData', final_tip(2), 'ZData', final_tip(3));

    set(ht1, 'String', sprintf(['%d-section arm (m=%.2f kg, xi=%.2f, sim time=%.4f s)\n' ...
        'stiff=%.1f | damp=%.1f\nframe %d / %d   t=%.2f s'], ...
        N, m, loc, times, stiff, damp, k, Nt, tk));

    drawnow limitrate
end

end