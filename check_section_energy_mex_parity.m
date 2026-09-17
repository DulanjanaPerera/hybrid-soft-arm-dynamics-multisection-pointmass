function parity=check_section_energy_mex_parity()
% Compare interpreted MATLAB and MEX trajectories at representative fits.
d=load(fullfile('cog_beta_energy_fit','fit.mat'),'result'); f=d.result;
s=load('comparison_after_C_correction.mat','X'); x0=s.X(1,:).';
p=f.parameters; parity.locations=[.1,.5,1];
parity.time=linspace(0,1,301).';
opts=odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
for j=1:3
    g=round(10*parity.locations(j)); p.cog_xi=parity.locations(j)*ones(3,1);
    beta=f.beta(:,:,4,g); p.beta=beta;
    args={p.L,p.r,p.cog_xi,p.mi(:),p.g(:),p.K,p.D,p.tau(:),p.mu,p.lKbounds(:),beta};
    mat=@(t,x) armS_dynamics_N3_entry_mex(t,x,args{:});
    mex=@(t,x) armS_dynamics_N3_entry_mex_mex(t,x,args{:});
    tic; [tm,Xm]=ode15s(mat,parity.time,x0,opts); parity.matlabSeconds(j)=toc;
    tic; [tx,Xx]=ode15s(mex,parity.time,x0,opts); parity.mexSeconds(j)=toc;
    assert(isequal(tm,tx) && all(isfinite([Xm(:);Xx(:)])));
    parity.maxCoordinateError_m(j)=max(abs(Xm(:,1:6)-Xx(:,1:6)),[],'all');
    parity.maxVelocityError_mps(j)=max(abs(Xm(:,7:12)-Xx(:,7:12)),[],'all');
    parity.maxScaledStateError(j)=max(abs(Xm-Xx)./(1+abs(Xm)),[],'all');
    assert(parity.maxScaledStateError(j)<1e-7,'MATLAB/MEX trajectory mismatch.');
    fprintf('xi %.1f: MATLAB %.3fs MEX %.3fs, max scaled state gap %.3e\n', ...
        parity.locations(j),parity.matlabSeconds(j),parity.mexSeconds(j), ...
        parity.maxScaledStateError(j));
end
save(fullfile('cog_beta_energy_fit','mex_parity.mat'),'parity');
end
