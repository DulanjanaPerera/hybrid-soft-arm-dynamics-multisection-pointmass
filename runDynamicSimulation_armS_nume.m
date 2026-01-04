clear
N = 4;
params.N = N;
params.L = 0.278;
params.r = 0.013;
params.cog_xi = 0.5*ones(1,N);
params.mi = 0.1*ones(N,1);
params.g = [0;0;9.81];
params.K = 1.2e3 * eye(2*N);
params.tau = zeros(2*N,1);
params.D = 30* eye(2*N);

% initial length changes
l0 = [
     0.010, 0.010;
     0.015, 0.015;
     0.001, 0.001;
     -0.01, -0.01;
    ];
fl0 = reshape(l0',[2*N,1]); % flattern the matrix
dl0 = zeros(2*N,1); % initial length change speed flatterned
X0 = [fl0;dl0];

% Integrate
tspan = [0 10];

% Stiffer ODE solver is used as the matrices are tend to go ill-condition
[t, X] = ode15s(@(t,X) armS_dynamics_nume(t,X,params), tspan, X0);