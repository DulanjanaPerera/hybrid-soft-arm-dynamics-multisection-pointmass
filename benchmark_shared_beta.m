function timing=benchmark_shared_beta()
% Same one-second IVP and solver settings, no plotting or code generation.
s=load('comparison_after_C_correction.mat','t','X','params');
f=load('shared_beta_fit.mat','result'); p=s.params; p.beta=f.result.betaMatrix;
t=s.t(s.t<=1); t=t(:); x0=s.X(1,:).';
opts=odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
rhs={@(tt,x) armS_standard_dynamics(tt,x,p), ...
    @(tt,x) armS_standard_mex(tt,x,p.L,p.r,p.mi(:),p.g(:), ...
        p.K,p.D,p.tau(:),p.mu,p.lKbounds(:)), ...
    @(tt,x) armS_dynamics_nume(tt,x,p), ...
    @(tt,x) armS_dynamics_N3_entry_mex_mex(tt,x,p.L,p.r, ...
        p.cog_xi(:),p.mi(:),p.g(:),p.K,p.D,p.tau(:),p.mu, ...
        p.lKbounds(:),p.beta)};
names={'standard MATLAB','standard MEX','fitted point MATLAB','fitted point MEX'};
times=zeros(4,3); paths=cell(1,4);
for i=1:4
    rhs{i}(t(1),x0); % warm both execution mode and binary
    for rep=1:3
        timer=tic;
        [~,paths{i}]=ode15s(rhs{i},t,x0,opts);
        times(i,rep)=toc(timer);
    end
    fprintf('%s: median %.4f s; runs %.4f %.4f %.4f s\n', ...
        names{i},median(times(i,:)),times(i,:));
end
timing.names=names;
timing.times=times;
timing.medians=median(times,2);
timing.maxMatlabMexStateDifference=[max(abs(paths{1}-paths{2}),[],'all'), ...
    max(abs(paths{3}-paths{4}),[],'all')];
fprintf('Maximum state MATLAB/MEX difference: standard %.3e, fitted point %.3e\n', ...
    timing.maxMatlabMexStateDifference);
save('shared_beta_timing.mat','timing');
end
