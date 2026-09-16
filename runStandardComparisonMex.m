function [t,X,params] = runStandardComparisonMex(referenceFile,matlabFile)
% Same initial state, parameters and solver settings as interpreted run.
if nargin < 1, referenceFile = 'comparison_after_C_correction.mat'; end
if nargin < 2, matlabFile = 'comparison_standard_distributed.mat'; end
assert(exist('armS_standard_mex','file')==3,'Run build_armS_standard_mex first.');
saved = load(referenceFile,'t','X','params');
params = saved.params;
mi = params.mi(:); g = params.g(:); tau = params.tau(:);
bounds = params.lKbounds(:);
X0 = saved.X(1,:).';
rhs = @(tt,xx) armS_standard_mex(tt,xx,params.L,params.r, ...
    mi,g,params.K,params.D,tau,params.mu,bounds);
% Initial-state agreement check also uses this recording's actual parameters.
y = armS_standard_dynamics(saved.t(1),X0,params);
z = rhs(saved.t(1),X0);
assert(all(isfinite(z)) && max(abs(y-z)./(1+abs(y)))<1e-8, ...
    'MEX differs from MATLAB at the initial state. Rebuild before continuing.');
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
rhs(saved.t(1),X0); % Warm up before timing.
clock = tic;
[t,X] = ode15s(rhs,saved.t(:),X0,opts);
params.times = toc(clock);
assert(numel(t)==numel(saved.t) && all(isfinite(X(:))), ...
    'Simulation did not return all requested finite samples.');
save('comparison_standard_distributed_mex.mat','t','X','params');
fprintf('MEX integration time: %.4f s\n',params.times);
if isfile(matlabFile)
    baseline = load(matlabFile,'t','X','params');
    same = isequal(t,baseline.t(:)) && isequal(X0,baseline.X(1,:).');
    fields = {'L','r','mi','g','K','D','tau','mu','lKbounds'};
    for j=1:numel(fields)
        f=fields{j}; same=same && isequal(params.(f),baseline.params.(f));
    end
    if same
        delta = X-baseline.X;
        fprintf('MATLAB/MEX trajectory max |dq|: %.3e m; max |dvelocity|: %.3e m/s\n', ...
            max(max(abs(delta(:,1:6)))),max(max(abs(delta(:,7:12)))));
        fprintf('Recorded MATLAB time: %.4f s; timing ratio: %.2fx\n', ...
            baseline.params.times,baseline.params.times/params.times);
    else
        warning('Baseline setup differs; trajectory comparison skipped.');
    end
end
drawingArms(t,X,1/60,params,0.001);
end
