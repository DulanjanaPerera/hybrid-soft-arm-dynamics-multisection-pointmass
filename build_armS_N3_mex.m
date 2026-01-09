function build_armS_N3_mex()
N  = 3;
nd = 2*N;

% Example inputs with correct sizes/types
t = 0;
X = zeros(2*nd,1);

L = 0.278;
r = 0.013;
cog_xi = 0.5*ones(N,1);
mi = 0.1*ones(N,1);
g = [0;0;-9.81];
Kmin = 2.2e1*eye(nd);
D = 10*eye(nd);
tau = zeros(nd,1);
mu = 2000;
lKbounds = [-0.02; 0.02; 1e6];

cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.TargetLang = 'C++';
cfg.EnableOpenMP = false; % keep false unless you explicitly parallelize

codegen -config cfg armS_dynamics_N3_entry_mex ...
    -args {t, X, L, r, cog_xi, mi, g, Kmin, D, tau, mu, lKbounds};
end
