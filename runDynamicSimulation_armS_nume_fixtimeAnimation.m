clear
N = 3;
loc = 0.5;
m = 0.1;
A = pi*(0.013/2)^2;
stiff = 2.2e3;
damp = 400;


params.loc = loc;
params.stiff = stiff;
params.damp = damp;
params.m = m;
params.N = N;
params.L = 0.278;
params.r = 0.013;
params.cog_xi = loc*ones(N,1);
% params.cog_xi = [0.5; 0.5; 0.5];
params.mi = m*ones(N,1);
% params.mi = [0.1; 0.1; 0.1];
params.g = [0; 0; 9.81];

params.K = stiff * eye(2*N);
% params.K(5,5) = 1000;
% params.K(6,6) = 1000;

% params.tau = zeros(2*N,1);
params.tau = -(A * 2.0 * 1e5)*ones(2*N,1);
params.D = damp * eye(2*N);
% params.D(5,5) = 100;
% params.D(6,6) = 100;

params.lKbounds = [-0.02; 0.02; 1e6]; % lmin, lmax, Kmax
params.mu = 2000;

l0 = -0.000* ones(N,2);
l0(1,:) = -0.01;
l0(2,:) = -0.01;
fl0 = reshape(l0',[2*N,1]); % flattern the matrix
dl0 = 0.000001 * ones(2*N,1); % initial length change speed flatterned
X0 = [fl0; dl0];

% Integrate
tspan = [0 10];

% Stiffer ODE solver is used as the matrices are tend to go ill-condition
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-2);
tic
% [t, X] = ode15s(@(t,X) armS_dynamics_nume(t,X,params), tspan, X0, opts);
% [t, X] = ode15s(@(t,X) armS_dynamics_N3_entry_mex(t, X, ...
%     params.L, params.r, params.cog_xi, params.mi, params.g, params.K, params.D, params.tau, params.mu, params.lKbounds), tspan, X0, opts);
[t, X] = ode15s(@(t,X) armS_dynamics_N3_entry_mex_mex(t, X, ...
    params.L, params.r, params.cog_xi, params.mi, params.g, params.K, params.D, params.tau, params.mu, params.lKbounds), tspan, X0, opts);
times = toc

params.times = times;
%% ==============================================================================
drawingArms(t, X, 0.001, params)