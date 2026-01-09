clear
N = 3;
loc = 0.1;
m = 1;
stiff = 2.2e1;
damp = 10;



params.N = N;
params.L = 0.278;
params.r = 0.013;
% params.cog_xi = loc*ones(N,1);
params.cog_xi = [0.5; 0.5; 0.5];
% params.mi = m*ones(N,1);
params.mi = [0.01; 1; 0.01];
params.g = [0; 0; -9.81];
params.K = stiff * eye(2*N);
params.tau = zeros(2*N,1);
params.D = damp * eye(2*N);
params.lKbounds = [-0.02; 0.02; 1e6]; % lmin, lmax, Kmax
params.mu = 2000;

% initial length changes
% l0 = [
%      0.010, 0.010;
%      0.015, 0.015;
%      0.001, 0.001;
%      -0.01, -0.01;
%     ];

l0 = -0.005* ones(N,2);
fl0 = reshape(l0',[2*N,1]); % flattern the matrix
dl0 = zeros(2*N,1); % initial length change speed flatterned
X0 = [fl0;dl0];

% Integrate
tspan = [0 10];

% Stiffer ODE solver is used as the matrices are tend to go ill-condition
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);
tic
% [t, X] = ode15s(@(t,X) armS_dynamics_nume_N3(t,X,params), tspan, X0, opts);
[t, X] = ode15s(@(t,X) armS_dynamics_N3_entry_mex_mex(t, X, ...
    params.L, params.r, params.cog_xi, params.mi, params.g, params.K, params.D, params.tau, params.mu, params.lKbounds), tspan, X0, opts);
times = toc

%% ===================== DRAW + PLOT STATES (SIMULTANEOUS ANIMATION) =====================
N = params.N;
L = params.L;
r = params.r;

xi = linspace(0,1,25);
Nt = length(t);
Nx = length(xi);

% Optional: rotate drawing 180deg about X if you want Z flipped
% Rx = [1 0 0; 0 -1 0; 0 0 -1];
useFlipZ = true;  % set true if you want 180 deg flip about X

%% --------------------- Figure 1: Arm drawing ---------------------
figure(1); clf
ax1 = axes; hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');
xlabel(ax1,'X'); ylabel(ax1,'Y'); zlabel(ax1,'Z');
rotate3d(ax1,'on');
str1 = sprintf('%d-section continuum arm with CoG at %.2f', N, loc);
ht1 = title(ax1,str1);

xlim(ax1, [-1 1]); ylim(ax1, [-1 1]); zlim(ax1, [-1.5 1.1]);

% one line per section (auto colors) + section tip marker
h_seg = gobjects(N,1);
h_tip_seg = gobjects(N,1);
segColor = zeros(N,3);

for s = 1:N
    h_seg(s) = plot3(ax1, NaN, NaN, NaN, 'LineWidth', 2);
    segColor(s,:) = h_seg(s).Color;
    h_tip_seg(s) = plot3(ax1, NaN, NaN, NaN, 's', ...
        'MarkerFaceColor', segColor(s,:), 'MarkerEdgeColor','none', 'MarkerSize', 7);
end

% overall tip (end of arm)
h_tip = plot3(ax1, NaN, NaN, NaN, 'o', ...
    'MarkerFaceColor',[.8 .2 .2], 'MarkerEdgeColor','none', 'MarkerSize', 7);

%% --------------------- Figure 2: Length-change plot ---------------------
figure(2); clf
ax2 = axes; hold(ax2,'on'); grid(ax2,'on');
axis(ax2,'tight');

% Plot ALL length-change signals (first 2N columns of X)
h_l = gobjects(2*N,1);
for i = 1:2*N
    h_l(i) = plot(ax2, t, X(:,i), 'LineWidth', 1.5);
end

str = sprintf('length change of %d Sections', N);
title(ax2, str);
xlabel(ax2, 'time (s)');
ylabel(ax2, 'length change (m)');

% Legend labels
leg = strings(2*N,1);
for ksec = 1:N
    leg(2*(ksec-1)+1) = "l_{" + ksec + "1}";
    leg(2*(ksec-1)+2) = "l_{" + ksec + "2}";
end
legend(ax2, leg, 'Location','best');

% Moving time cursor
yl = ylim(ax2);
h_cursor = plot(ax2, [t(1) t(1)], yl, 'k--', 'LineWidth', 1.2);

% Moving markers (dots) on each curve at current time
h_dot = gobjects(2*N,1);
for i = 1:2*N
    % use same color as the corresponding curve
    col = h_l(i).Color;
    h_dot(i) = plot(ax2, t(1), X(1,i), 'o', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
end

%% --------------------- Animation loop (updates both figures) ---------------------
for k = 1:Nt

    % ===== Update Figure 2 (time cursor + dots) =====
    tk = t(k);
    set(h_cursor, 'XData', [tk tk], 'YData', ylim(ax2));  % keep cursor spanning y-limits

    for i = 1:2*N
        set(h_dot(i), 'XData', tk, 'YData', X(k,i));
    end

    % ===== Update Figure 1 (arm drawing) =====
    l_mat = reshape(X(k, 1:2*N).', [2, N]).';  % N x 2 : [l2 l3]

    Rglob = eye(3);
    Pglob = zeros(3,1);
    final_tip = [NaN;NaN;NaN];

    for s = 1:N
        length_s = [0, l_mat(s,1), l_mat(s,2)];

        % sample this section backbone in GLOBAL
        Xs = zeros(Nx,1); Ys = zeros(Nx,1); Zs = zeros(Nx,1);
        for j = 1:Nx
            [~, ~, Ploc] = HTM_nume(length_s, xi(j), L, r);
            Ppoint = Pglob + Rglob * Ploc;

            if useFlipZ
                Ppoint = [1 0 0; 0 -1 0; 0 0 -1] * Ppoint; %#ok<MINV>
            end

            Xs(j) = Ppoint(1); Ys(j) = Ppoint(2); Zs(j) = Ppoint(3);
        end

        set(h_seg(s), 'XData', Xs, 'YData', Ys, 'ZData', Zs);

        % section tip
        [~, Rtip, Ptip] = HTM_nume(length_s, 1, L, r);
        Pglob_tip = Pglob + Rglob * Ptip;

        if useFlipZ
            Pglob_tip = [1 0 0; 0 -1 0; 0 0 -1] * Pglob_tip; %#ok<MINV>
        end

        set(h_tip_seg(s), 'XData', Pglob_tip(1), 'YData', Pglob_tip(2), 'ZData', Pglob_tip(3));

        % advance to next section (NOTE: advance uses unflipped transforms)
        % so do NOT use flipped Pglob_tip here; use the original:
        Pglob = Pglob + Rglob * Ptip;
        Rglob = Rglob * Rtip;

        final_tip = Pglob;
        if useFlipZ
            final_tip = [1 0 0; 0 -1 0; 0 0 -1] * final_tip; %#ok<MINV>
        end
    end

    set(h_tip, 'XData', final_tip(1), 'YData', final_tip(2), 'ZData', final_tip(3));

    set(ht1, 'String', sprintf('%d-section arm (m=%.2f kg, xi = %.2f, sim time = %.4f s) \n stiff = %.1f | damp = %.1f \n frame %d / %d   t = %.3f s', ...
        N, m, loc, times, stiff, damp, k, Nt, t(k)));

    drawnow
    % pause(0.001);
end
%% ==============================================================================