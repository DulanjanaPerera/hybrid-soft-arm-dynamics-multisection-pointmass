% Simulate the standard distributed-mass arm at fixed base orientations.
% Each row of orientationRPYDeg is a separate stationary-base solve.
% Angles are [roll pitch yaw] in degrees; R = Rz(yaw)*Ry(pitch)*Rx(roll).
% This script leaves results, X0, params, and orientationRPYDeg in the workspace.

params.N = 3;
params.L = 0.278;                 % Section length (m)
params.r = 0.013;                 % Arm geometry parameter (m)
params.mi = [0.1;0.1;0.1];        % Distributed section masses (kg)
params.gWorld = [0;0;9.81];     % Gravity in world coordinates (m/s^2)
params.K = 2200*eye(6);
params.D = 600*eye(6);
params.mu = 2000;
params.lKbounds = [-0.02;0.02;1e6];

pressureBar = 0;
area = pi*(params.r/2)^2;
params.tau = -(area*pressureBar*1e5)*ones(6,1);

% The base remains stationary during each solve. Add or edit rows to compare
% orientations; measured 3x3 rotations can be used directly in a custom RHS.
orientationRPYDeg = [0 0 0;
                     90 0 0;
                     0 90 0];

q0BySection = [-0.008 0.003;
               -0.004 -0.006;
                0.002 -0.005];
dq0BySection = zeros(3,2);
q0 = reshape(q0BySection.',6,1);
dq0 = reshape(dq0BySection.',6,1);
X0 = [q0;dq0];

durationSeconds = 0.1;
outputIntervalSeconds = 0.01;
tspan = 0:outputIntervalSeconds:durationSeconds;
odeOptions = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
useMex = true;                    % false uses MATLAB source
showPlots = true;                 % coordinate histories and final shapes

assert(isequal(size(params.mi),[3,1]) && all(isfinite(params.mi)) ...
    && all(params.mi>0),'mi must contain three positive masses.');
assert(isequal(size(params.gWorld),[3,1]) && all(isfinite(params.gWorld)));
assert(isequal(size(params.K),[6,6]) && isequal(size(params.D),[6,6]) ...
    && isequal(size(params.tau),[6,1]));
assert(isequal(size(X0),[12,1]) && all(isfinite(X0)));
assert(size(orientationRPYDeg,2)==3 && ...
    all(isfinite(orientationRPYDeg(:))), ...
    'Orientations must have [roll pitch yaw] rows in degrees.');
assert(numel(tspan)>=2 && all(diff(tspan)>0));
if useMex
    assert(exist('armS_stationary_base_mex','file')==3, ...
        ['Stationary-base MEX missing. Run build_armS_stationary_base_mex ' ...
         'or set useMex=false.']);
end

numCases = size(orientationRPYDeg,1);
results = repmat(struct('label','','R_world_from_arm',eye(3), ...
    'g_arm',zeros(3,1),'t',[],'X',[],'elapsedSeconds',0),numCases,1);
for k=1:numCases
    angles = deg2rad(orientationRPYDeg(k,:));
    cr=cos(angles(1)); sr=sin(angles(1));
    cp=cos(angles(2)); sp=sin(angles(2));
    cy=cos(angles(3)); sy=sin(angles(3));
    Rx=[1 0 0;0 cr -sr;0 sr cr];
    Ry=[cp 0 sp;0 1 0;-sp 0 cp];
    Rz=[cy -sy 0;sy cy 0;0 0 1];
    R_world_from_arm=Rz*Ry*Rx;
    results(k).label=sprintf('RPY [%g %g %g] deg',orientationRPYDeg(k,:));
    results(k).R_world_from_arm=R_world_from_arm;
    results(k).g_arm=R_world_from_arm.'*params.gWorld;

    if useMex
        rhs=@(tt,x) armS_stationary_base_mex(tt,x,params.L,params.r, ...
            params.mi,R_world_from_arm,params.gWorld,params.K, ...
            params.D,params.tau,params.mu,params.lKbounds);
    else
        rhs=@(tt,x) armS_stationary_base_entry(tt,x,params.L,params.r, ...
            params.mi,R_world_from_arm,params.gWorld,params.K, ...
            params.D,params.tau,params.mu,params.lKbounds);
    end
    solveClock=tic;
    [results(k).t,results(k).X]=ode15s(rhs,tspan,X0,odeOptions);
    results(k).elapsedSeconds=toc(solveClock);
    assert(all(isfinite(results(k).X(:))) && ...
        abs(results(k).t(end)-tspan(end))<1e-12, ...
        'Simulation did not reach the requested endpoint.');
    fprintf('%s: %d output points, %.3f s solve, final |X-X0| %.6e\n', ...
        results(k).label,numel(results(k).t), ...
        results(k).elapsedSeconds,norm(results(k).X(end,:)'-X0));
end

if showPlots
    colors=lines(numCases);
    coordinateLabels={'l_{12}','l_{13}','l_{22}','l_{23}', ...
        'l_{32}','l_{33}'};
    figure('Name','Stationary-base standard arm coordinates');
    tiledlayout(3,2);
    for i=1:6
        ax=nexttile; hold(ax,'on'); grid(ax,'on');
        for k=1:numCases
            plot(ax,results(k).t,results(k).X(:,i), ...
                'Color',colors(k,:),'LineWidth',1.3, ...
                'DisplayName',results(k).label);
        end
        xlabel(ax,'Time (s)'); ylabel(ax,'Length change (m)');
        title(ax,coordinateLabels{i});
        if i==1, legend(ax,'Location','best'); end
    end

    figure('Name','Stationary-base final arm shapes');
    ax=axes; hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    xlabel(ax,'World X (m)'); ylabel(ax,'World Y (m)');
    zlabel(ax,'World Z (m)'); view(ax,3);
    for k=1:numCases
        q=results(k).X(end,1:6).';
        Rroot=results(k).R_world_from_arm;
        origin=zeros(3,1); Rsection=eye(3);
        for n=1:3
            l=[0,q(2*n-1:2*n).'];
            xyz=zeros(3,25);
            for j=1:25
                xi=(j-1)/24;
                [~,~,pLocal]=HTM_nume(l,xi,params.L,params.r);
                xyz(:,j)=Rroot*(origin+Rsection*pLocal);
            end
            if n==1
                plot3(ax,xyz(1,:),xyz(2,:),xyz(3,:), ...
                    'Color',colors(k,:),'LineWidth',1.8, ...
                    'DisplayName',results(k).label);
            else
                plot3(ax,xyz(1,:),xyz(2,:),xyz(3,:), ...
                    'Color',colors(k,:),'LineWidth',1.8, ...
                    'HandleVisibility','off');
            end
            [~,Rtip,pTip]=HTM_nume(l,1,params.L,params.r);
            origin=origin+Rsection*pTip;
            Rsection=Rsection*Rtip;
        end
    end
    plot3(ax,0,0,0,'ko','MarkerFaceColor','k', ...
        'DisplayName','Base origin');
    legend(ax,'Location','best');
end
