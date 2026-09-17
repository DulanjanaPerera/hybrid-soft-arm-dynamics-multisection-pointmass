function report=validate_section_energy_beta()
% Validate energy-fitted beta against fixed standard dynamics and prior sweep.
d=load(fullfile('cog_beta_energy_fit','fit.mat'),'result'); fit=d.result;
old=load(fullfile('cog_sweep_fixed_beta','cog_sweep.mat'),'sweep'); s=old.sweep;
p=fit.parameters; outdir=fullfile(pwd,'cog_beta_energy_fit');
assert(isequal(fit.locations,s.locations));
report.created=char(datetime('now')); report.locations=fit.locations;
report.labels={'all ones','fixed midpoint shared','new CoG shared','new CoG section'};
report.beta=fit.beta; report.fit=fit;
sample=1:5:numel(s.time); report.time=s.time(sample);
report.forceRelative=zeros(10,4); report.forceRms=zeros(10,4);
report.metrics=struct([]); row=0;
for g=1:10
    pp=p; pp.cog_xi=fit.locations(g)*ones(3,1);
    [report.forceRelative(g,:),report.forceRms(g,:)]=forceCheck(pp,fit,g);
    check=coreCheck(pp,fit.beta(:,:,4,g),fit.heldPoses(:,2));
    if g==1, report.coreChecks=check; else, report.coreChecks(g)=check; end
    for k=1:4
        c=s.cases(k); stdX=c.standardX(sample,:);
        for b=1:4
            beta=fit.beta(:,:,b,g); pp.beta=beta;
            if fit.condition(g,b).minRcond<1e-8
                r.X=[]; r.beta=beta; r.xi=fit.locations(g);
                r.skipped='Fit-pose rcond(M) below 1e-8; trajectory numerically unreliable';
                r.qDifference=nan(numel(sample),6);
                r.tipDifference=nan(numel(sample),3);
                r.kineticDifference=nan(numel(sample),3);
                r.gravityDifference=nan(numel(sample),3);
                r.totalDifference=nan(numel(sample),3);
                r.sameStateKineticDifference=nan(numel(sample),3);
                r.totalChangeDifference=nan(numel(sample),3);
                r.minRcond=nan;
                report.cases(k,g,b)=r;
                for n=1:3
                    row=row+1; item=emptyMetric(k,r.xi,b,n);
                    if row==1, report.metrics=item; else, report.metrics(row)=item; end
                end
                continue
            end
            if b==2
                X=s.cases(k).runs(g).X;
            elseif b>2 && isequal(beta,fit.beta(:,:,1,g))
                X=report.cases(k,g,1).X;
            elseif b==4 && isequal(beta,fit.beta(:,:,3,g))
                X=report.cases(k,g,3).X;
            else
                rhs=@(tt,x) armS_dynamics_N3_entry_mex_mex(tt,x, ...
                    pp.L,pp.r,pp.cog_xi,pp.mi(:),pp.g(:),pp.K,pp.D,c.tau, ...
                    pp.mu,pp.lKbounds(:),pp.beta);
                [times,X]=ode15s(rhs,s.time,c.x0, ...
                    odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3));
                assert(isequal(times,s.time) && all(isfinite(X(:))));
            end
            % RHS uses runtime beta and must match interpreted MATLAB entry.
            x0=c.x0; rhsValue=armS_dynamics_N3_entry_mex_mex(s.time(1),x0, ...
                pp.L,pp.r,pp.cog_xi,pp.mi(:),pp.g(:),pp.K,pp.D,c.tau, ...
                pp.mu,pp.lKbounds(:),pp.beta);
            matValue=armS_dynamics_N3_entry_mex(s.time(1),x0, ...
                pp.L,pp.r,pp.cog_xi,pp.mi(:),pp.g(:),pp.K,pp.D,c.tau, ...
                pp.mu,pp.lKbounds(:),pp.beta);
            assert(max(abs(rhsValue-matValue)./(1+abs(matValue)))<1e-8);
            r.X=X; r.beta=beta; r.xi=fit.locations(g); r.skipped='';
            r.qDifference=X(sample,1:6)-stdX(:,1:6);
            r.tipDifference=zeros(numel(sample),3);
            r.kineticDifference=zeros(numel(sample),3);
            r.gravityDifference=zeros(numel(sample),3);
            r.totalDifference=zeros(numel(sample),3);
            r.sameStateKineticDifference=zeros(numel(sample),3);
            if b==2
                prev=c.runs(g); r.tipDifference=prev.tipDifference(sample,:);
                r.kineticDifference=prev.energyDifference.kinetic(sample,:);
                r.gravityDifference=prev.energyDifference.gravity(sample,:);
                r.totalDifference=prev.energyDifference.total(sample,:);
                r.sameStateKineticDifference=prev.sameStateEnergyDifference.kinetic(sample,:);
            else
                for j=1:numel(sample)
                    point=arm_module_observables(X(sample(j),:).',pp,'pointmass');
                    same=arm_module_observables(stdX(j,:).',pp,'pointmass');
                    stdtip=squeeze(c.standard.tip(sample(j),:,:));
                    r.tipDifference(j,:)=sqrt(sum((point.tip-stdtip).^2,1));
                    r.kineticDifference(j,:)=point.kinetic-c.standard.kinetic(sample(j),:);
                    r.gravityDifference(j,:)=point.gravity-c.standard.gravity(sample(j),:);
                    r.totalDifference(j,:)=point.total-c.standard.total(sample(j),:);
                    r.sameStateKineticDifference(j,:)=same.kinetic-c.standard.kinetic(sample(j),:);
                end
            end
            r.totalChangeDifference=r.totalDifference-r.totalDifference(1,:);
            r.minRcond=inf;
            for j=round(linspace(1,size(X,1),11))
                q=X(j,1:6).'; v=X(j,7:12).';
                M=pointCore(q,v,pp,beta);
                [~,flag]=chol((M+M.')/2);
                assert(flag==0,'Trajectory M not SPD at xi %.1f, case %d, beta %d',r.xi,k,b);
                r.minRcond=min(r.minRcond,rcond(M));
            end
            report.cases(k,g,b)=r;
            for n=1:3
                row=row+1; ix=2*n-1:2*n;
                item.caseIndex=k; item.xi=r.xi; item.candidateIndex=b; item.module=n;
                item.coordinateMax_m=max(abs(r.qDifference(:,ix)),[],'all');
                item.tipMax_m=max(r.tipDifference(:,n));
                for field={'kinetic','gravity','total','totalChange','sameStateKinetic'}
                    f=field{1}; delta=r.([f,'Difference'])(:,n);
                    item.([f,'Rms_J'])=sqrt(trapz(report.time,delta.^2)/(report.time(end)-report.time(1)));
                    item.([f,'Max_J'])=max(abs(delta));
                end
                item.minRcond=r.minRcond;
                if row==1, report.metrics=item; else, report.metrics(row)=item; end
            end
        end
    end
    fprintf('xi %.1f: held C*dq relative [%.3f %.3f %.3f %.3f], min trajectory rcond %.3e\n', ...
        fit.locations(g),report.forceRelative(g,:),min([report.metrics([report.metrics.xi]==fit.locations(g)).minRcond],[],'omitnan'));
    save(fullfile(outdir,'validation.mat'),'report','-v7.3');
end

function item=emptyMetric(k,xi,b,n)
item.caseIndex=k; item.xi=xi; item.candidateIndex=b; item.module=n;
item.coordinateMax_m=nan; item.tipMax_m=nan;
for field={'kinetic','gravity','total','totalChange','sameStateKinetic'}
    f=field{1}; item.([f,'Rms_J'])=nan; item.([f,'Max_J'])=nan;
end
item.minRcond=nan;
end
end

function [relative,rms]=forceCheck(p,fit,g)
qset=fit.heldPoses; V=fit.heldVelocities;
numer=zeros(1,4); denom=0; count=0;
for j=1:size(qset,2)
    q=qset(:,j);
    for vj=1:size(V,2)
        v=V(:,vj,j); [~,Cs]=armS_standard_core(q,v,p);
        target=Cs*v; denom=denom+sum(target.^2); count=count+1;
        for b=1:4
            [~,C]=pointCore(q,v,p,fit.beta(:,:,b,g));
            numer(b)=numer(b)+sum((C*v-target).^2);
        end
    end
end
relative=sqrt(numer/denom); rms=sqrt(numer/count);
end

function c=coreCheck(p,beta,q)
v=[.031;-.017;.024;-.016;.013;-.021];
[M,C,~,dM]=pointCore(q,v,p,beta);
c.symmetry=norm(M-M.','fro'); c.skew=norm(sum(dM.*reshape(v,1,1,6),3)-C-C.','fro');
c.derivative=zeros(1,6); h=1e-6;
for j=1:6
    qp=q; qm=q; qp(j)=qp(j)+h; qm(j)=qm(j)-h;
    Mp=pointCore(qp,v,p,beta); Mm=pointCore(qm,v,p,beta);
    c.derivative(j)=norm((Mp-Mm)/(2*h)-dM(:,:,j),'fro')/max(1,norm(dM(:,:,j),'fro'));
end
c.minRcond=rcond(M); [~,c.cholFlag]=chol((M+M.')/2);
assert(c.symmetry<1e-9 && c.skew<1e-8 && max(c.derivative)<1e-6 && c.cholFlag==0);
end

function [M,C,G,dM]=pointCore(q,v,p,beta)
[M,C,G,dM]=armS_core_N3_mex(0,reshape(q,2,3).',reshape(v,2,3).', ...
    p.L,p.r,p.cog_xi,p.mi,p.g,p.K,beta);
end
