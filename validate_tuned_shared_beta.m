function report=validate_tuned_shared_beta(fitFile)
% Held-out force and trajectory checks for the shared-beta fit.
if nargin<1, fitFile='shared_beta_fit.mat'; end
s=load(fitFile,'result'); fit=s.result;
p=fit.parameters; betaFit=fit.betaMatrix; betaOne=ones(3,3);
qHeld=fit.heldPoses; n=size(qHeld,2);
rngState=rng; cleanup=onCleanup(@() rng(rngState)); %#ok<NASGU>
rng(1887); velocities=.04*rand(6,n)-.02;
forceError=zeros(n,2); forceReference=zeros(n,1);
massRcond=zeros(n,3); cholFlags=zeros(n,3);
for k=1:n
    q=qHeld(:,k); dq=velocities(:,k);
    [Ms,Cs]=armS_standard_core(q,dq,p);
    fref=Cs*dq; forceReference(k)=norm(fref);
    [~,cholFlags(k,1)]=chol((Ms+Ms.')/2);
    massRcond(k,1)=rcond(Ms);
    for b=1:2
        if b==1, beta=betaOne; else, beta=betaFit; end
        [M,C]=pointCore(q,dq,p,beta);
        forceError(k,b)=norm(C*dq-fref);
        massRcond(k,b+1)=rcond(M);
        [~,cholFlags(k,b+1)]=chol((M+M.')/2);
    end
end
report.forceRms=sqrt(mean(forceError.^2,1));
report.forceRelative=sqrt(sum(forceError.^2)/sum(forceReference.^2));
report.forceReferenceRms=sqrt(mean(forceReference.^2));
report.forceErrors=forceError;
report.heldVelocities=velocities;
report.heldMassRcond=massRcond;
report.heldCholFlags=cholFlags;
assert(all(cholFlags(:)==0),'A held-out mass matrix is not positive definite.');
fprintf('Held-out C*dq RMS force: %.3e -> %.3e; aggregate relative %.3f -> %.3f\n', ...
    report.forceRms,report.forceRelative);
fprintf('Minimum held-out rcond [standard, beta=1, fitted]: %.3e %.3e %.3e\n', ...
    min(massRcond,[],1));

ref=load('comparison_after_C_correction.mat','t','X','params');
t=ref.t(:); q0=ref.X(1,1:6).'; dq0=ref.X(1,7:12).';
starts=[q0,[-.008;.003;-.004;-.006;.002;-.005]];
names={'reference bend','asymmetric bend'};
inputs=[p.tau(:),p.tau(:)+[.5;-.5;0;0;0;0]];
inputNames={'free','differential'};
opts=odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
cases=struct([]); idx=0;
for u=1:2
    for k=1:2
        idx=idx+1; x0=[starts(:,k);dq0]; tau=inputs(:,u);
        std=@(tt,x) armS_standard_mex(tt,x,p.L,p.r,p.mi(:),p.g(:), ...
            p.K,p.D,tau,p.mu,p.lKbounds(:));
        point=@(beta) @(tt,x) armS_dynamics_N3_entry_mex_mex(tt,x, ...
            p.L,p.r,p.cog_xi(:),p.mi(:),p.g(:),p.K,p.D,tau,p.mu, ...
            p.lKbounds(:),beta);
        [ts,xs]=ode15s(std,t,x0,opts);
        [to,xo]=ode15s(point(betaOne),t,x0,opts);
        [tf,xf]=ode15s(point(betaFit),t,x0,opts);
        assert(isequal(ts,to,tf) && all(isfinite([xs(:);xo(:);xf(:)])));
        c.name=sprintf('%s / %s',names{k},inputNames{u});
        c.standard=xs; c.betaOne=xo; c.fitted=xf;
        c.coordinateMax=[max(abs(xo(:,1:6)-xs(:,1:6)),[],'all'), ...
            max(abs(xf(:,1:6)-xs(:,1:6)),[],'all')];
        c.coordinateRms=[sqrt(mean((xo(:,1:6)-xs(:,1:6)).^2,'all')), ...
            sqrt(mean((xf(:,1:6)-xs(:,1:6)).^2,'all'))];
        tipError=zeros(numel(t),2);
        for j=1:numel(t)
            tipStd=tipPosition(xs(j,1:6).',p);
            tipError(j,1)=norm(tipPosition(xo(j,1:6).',p)-tipStd);
            tipError(j,2)=norm(tipPosition(xf(j,1:6).',p)-tipStd);
        end
        c.tipMax=max(tipError,[],1);
        c.tipRms=sqrt(mean(tipError.^2,1));
        c.tipError=tipError;
        if idx==1, cases=c; else, cases(idx)=c; end
        fprintf('%s: max q %.3f -> %.3f mm; max tip %.1f -> %.1f mm\n', ...
            c.name,1000*c.coordinateMax,1000*c.tipMax);
    end
end
report.time=t; report.cases=cases; report.betaFit=betaFit;
if strcmp(fitFile,'shared_beta_fit.mat')
    outputFile='shared_beta_validation.mat';
else
    [~,stem]=fileparts(fitFile);
    outputFile=[strrep(stem,'fit','validation'),'.mat'];
end
save(outputFile,'report');
end

function [M,C]=pointCore(q,dq,p,beta)
[M,C]=armS_core_N3_mex(0,reshape(q,2,3).',reshape(dq,2,3).', ...
    p.L,p.r,p.cog_xi,p.mi,p.g,p.K,beta);
end

function tip=tipPosition(q,p)
R=eye(3); tip=zeros(3,1);
for n=1:3
    [~,Rn,pn]=HTM_nume([0,q(2*n-1:2*n).'],1,p.L,p.r);
    tip=tip+R*pn; R=R*Rn;
end
end
