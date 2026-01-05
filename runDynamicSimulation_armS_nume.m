clear
N = 4;
params.N = N;
params.L = 0.278;
params.r = 0.013;
params.cog_xi = 0.1*ones(1,N);
params.mi = 0.2*ones(N,1);
params.g = [0;0;-9.81];
params.K = 1.2e3 * eye(2*N);
params.tau = zeros(2*N,1);
params.D = 1000* eye(2*N);

% initial length changes
% l0 = [
%      0.010, 0.010;
%      0.015, 0.015;
%      0.001, 0.001;
%      -0.01, -0.01;
%     ];

l0 = 0.005* ones(N,2);
fl0 = reshape(l0',[2*N,1]); % flattern the matrix
dl0 = zeros(2*N,1); % initial length change speed flatterned
X0 = [fl0;dl0];

% Integrate
tspan = [0 10];

% Stiffer ODE solver is used as the matrices are tend to go ill-condition
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);
tic
[t, X] = ode15s(@(t,X) armS_dynamics_nume(t,X,params), tspan, X0, opts);
toc

%% ===================== DRAW N-SECTION ARM (colored segments) =====================
% Assumes you already have: t, X from ode15s, and params with N,L,r
N = params.N;
L = params.L;
r = params.r;

xi = linspace(0,1,25);
Nt = length(t);
Nx = length(xi);

% --- Figure setup ---
figure(1); clf
ax = axes; hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
rotate3d(ax,'on');
ht = title(ax,'N-section continuum arm backbone');

% axis limits (adjust as needed)
xlim(ax, [-1 1]); ylim(ax, [-1 1]); zlim(ax, [-0.1 1.7]);

% --- One line per section (different colors automatically) ---
h_seg = gobjects(N,1);
h_tip_seg = gobjects(N,1);

for s = 1:N
    h_seg(s) = plot3(ax, NaN, NaN, NaN, 'LineWidth', 2);          % line color auto
    col = h_seg(s).Color;                                        % capture auto color
    h_tip_seg(s) = plot3(ax, NaN, NaN, NaN, 's', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
end

% overall tip (end of arm)
h_tip = plot3(ax, NaN, NaN, NaN, 'o', ...
    'MarkerFaceColor',[.8 .2 .2], 'MarkerEdgeColor','none', 'MarkerSize', 7);

% --- Animation loop ---
for k = 1:Nt

    % generalized coordinates at this time step (N x 2): [l2 l3]
    l_mat = reshape(X(k, 1:2*N).', [2, N]).';

    % global transform at base of section 1
    Rglob = eye(3);
    Pglob = zeros(3,1);

    % keep track of final tip
    final_tip = [NaN;NaN;NaN];

    for s = 1:N
        length_s = [0, l_mat(s,1), l_mat(s,2)];

        % sample this section backbone in GLOBAL
        Xs = zeros(Nx,1); Ys = zeros(Nx,1); Zs = zeros(Nx,1);
        for j = 1:Nx
            [~, ~, Ploc] = HTM_nume(length_s, xi(j), L, r) ;
            Ppoint =  (Pglob + Rglob * Ploc);
            Xs(j) = Ppoint(1); Ys(j) = Ppoint(2); Zs(j) = Ppoint(3);
        end

        % update this section line
        set(h_seg(s), 'XData', Xs, 'YData', Ys, 'ZData', Zs);

        % compute section tip transform (xi = 1) and update section-tip marker
        [~, Rtip, Ptip] = HTM_nume(length_s, 1, L, r);
        Pglob_tip = (Pglob + Rglob * Ptip);              % tip in global
        set(h_tip_seg(s), 'XData', Pglob_tip(1), 'YData', Pglob_tip(2), 'ZData', Pglob_tip(3));

        % advance base frame to next section
        Pglob = Pglob_tip;
        Rglob = Rglob * Rtip;

        final_tip = Pglob_tip;
    end

    % update overall tip marker (end effector)
    set(h_tip, 'XData', final_tip(1), 'YData', final_tip(2), 'ZData', final_tip(3));

    set(ht, 'String', sprintf('N-section arm — frame %d / %d   t = %.3f s', ...
        k, Nt, t(k)));
    % view([45, 45])
    drawnow
    pause(0.001);
end
% ==============================================================================

figure(2)
hold on
for i=1:N
    plot(t, X(:,i));
end
hold off
grid on
axis tight
str = sprintf('length change of %d Sections', N);
title(str)
xlabel 'time (s)'
ylabel 'length change (m)'

leg = strings(2*N,1);
for k=1:N
    leg(2*(k-1)+1) = "l_{" + k + "1}";
    leg(2*(k-1)+2) = "l_{" + k + "2}";
end

legend(leg)