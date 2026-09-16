function [t,X,params] = runStandardComparison(referenceFile)
% Use precisely the saved point-mass run's parameters, start state and times.
% The interpreted model is for validation; compile only after checks pass.
if nargin<1, referenceFile='comparison_after_C_correction.mat'; end
saved=load(referenceFile,'t','X','params');
params=saved.params;
validate_armS_standard(params);
opts=odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
clock=tic;
[t,X]=ode15s(@(t,x) armS_standard_dynamics(t,x,params), ...
    saved.t(:),saved.X(1,:).',opts);
params.times=toc(clock);
save('comparison_standard_distributed.mat','t','X','params');
fprintf('Saved comparison_standard_distributed.mat (%.2f s computation).\n',params.times);
drawingArms(t,X,1/60,params,0.001);
end
