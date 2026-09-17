function sweep=run_cog_sweep()
% Fixed shared beta; vary all three selected mass positions together.
fit=load('shared_beta_fit_bound5.mat','result');
ref=load('comparison_after_C_correction.mat','t','X','params');
p=ref.params; p.beta=fit.result.betaMatrix;
assert(isequal(p.mi,fit.result.parameters.mi) && p.N==3);
outdir=fullfile(pwd,'cog_sweep_fixed_beta');
if ~isfolder(outdir), mkdir(outdir); end
verifyObservables(p);
sweep.beta=p.beta; sweep.parameters=p;
sweep.locations=(1:10)/10;
sweep.time=linspace(ref.t(1),ref.t(end),1501).'; % 300 Hz, includes saved times
sweep.definition='pointmass minus standard; all section positions swept together';
sweep.created=char(datetime('now'));
starts=[ref.X(1,1:6).',[-.008;.003;-.004;-.006;.002;-.005]];
names={'Reference bend','Asymmetric bend'};
inputs=[p.tau(:),p.tau(:)+[.5;-.5;0;0;0;0]];
inputNames={'free','differential'};
t=sweep.time; nt=numel(t); opts=odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
summary=struct([]); row=0; caseIndex=0;
for u=1:2
    for start=1:2
        caseIndex=caseIndex+1;
        x0=[starts(:,start);ref.X(1,7:12).']; tau=inputs(:,u);
        c=struct(); c.name=[names{start},' / ',inputNames{u}]; c.x0=x0; c.tau=tau;
        std=@(tt,x) armS_standard_mex(tt,x,p.L,p.r,p.mi(:),p.g(:), ...
            p.K,p.D,tau,p.mu,p.lKbounds(:));
        [ts,c.standardX]=ode15s(std,t,x0,opts);
        assert(isequal(ts,t));
        c.standard=observablesSeries(c.standardX,p,'standard');
        for g=1:10
            pp=p; pp.cog_xi=sweep.locations(g)*ones(3,1);
            rhs=@(tt,x) armS_dynamics_N3_entry_mex_mex(tt,x,pp.L,pp.r, ...
                pp.cog_xi,pp.mi(:),pp.g(:),pp.K,pp.D,tau,pp.mu,pp.lKbounds(:),pp.beta);
            % Verify runtime CoG and beta against interpreted entry at each location.
            y=armS_dynamics_N3_entry_mex(t(1),x0,pp.L,pp.r,pp.cog_xi, ...
                pp.mi(:),pp.g(:),pp.K,pp.D,tau,pp.mu,pp.lKbounds(:),pp.beta);
            z=rhs(t(1),x0);
            assert(max(abs(y-z)./(1+abs(y)))<1e-8);
            r=struct(); r.xi=sweep.locations(g); timer=tic;
            [tp,r.X]=ode15s(rhs,t,x0,opts); r.integrationSeconds=toc(timer);
            assert(isequal(tp,t) && all(isfinite(r.X(:))));
            r.point=observablesSeries(r.X,pp,'pointmass');
            r.sameStatePoint=observablesSeries(c.standardX,pp,'pointmass');
            r.coordinateDifference=r.X(:,1:6)-c.standardX(:,1:6);
            r.tipVectorDifference=r.point.tip-c.standard.tip;
            r.tipDifference=squeeze(sqrt(sum(r.tipVectorDifference.^2,2)));
            for field={'kinetic','gravity','elastic','total'}
                f=field{1};
                r.energyDifference.(f)=r.point.(f)-c.standard.(f);
                r.sameStateEnergyDifference.(f)=r.sameStatePoint.(f)-c.standard.(f);
            end
            r.energyDifference.totalChange=r.energyDifference.total-r.energyDifference.total(1,:);
            % Sample mass conditioning along each simulated path (31 checkpoints).
            sample=unique(round(linspace(1,nt,31))); r.massRcond=zeros(size(sample));
            r.cholFlag=zeros(size(sample));
            for z=1:numel(sample)
                x=r.X(sample(z),:).';
                M=armS_core_N3_mex(t(sample(z)),reshape(x(1:6),2,3).', ...
                    reshape(x(7:12),2,3).',pp.L,pp.r,pp.cog_xi,pp.mi,pp.g,pp.K,pp.beta);
                r.massRcond(z)=rcond(M); [~,r.cholFlag(z)]=chol((M+M.')/2);
            end
            assert(all(r.cholFlag==0),'Indefinite M for %s xi %.1f',c.name,r.xi);
            r.conditionTimes=t(sample);
            for n=1:3
                row=row+1; ix=2*n-1:2*n;
                item.caseIndex=caseIndex; item.caseName=c.name; item.xi=r.xi; item.module=n;
                qdiff=r.coordinateDifference(:,ix);
                item.coordinateMaxEach_m=max(abs(qdiff),[],1);
                item.coordinateRmsEach_m=sqrt(trapz(t,qdiff.^2,1)/(t(end)-t(1)));
                item.coordinateMax_m=max(abs(qdiff),[],'all');
                item.coordinateRms_m=sqrt(trapz(t,sum(qdiff.^2,2))/(2*(t(end)-t(1))));
                [item.tipMax_m,peak]=max(r.tipDifference(:,n)); item.tipMaxTime_s=t(peak);
                item.tipRms_m=sqrt(trapz(t,r.tipDifference(:,n).^2)/(t(end)-t(1)));
                for field={'kinetic','gravity','elastic','total','totalChange'}
                    f=field{1}; delta=r.energyDifference.(f)(:,n);
                    item.([f,'MaxAbs_J'])=max(abs(delta));
                    item.([f,'Rms_J'])=sqrt(trapz(t,delta.^2)/(t(end)-t(1)));
                end
                delta=r.sameStateEnergyDifference.kinetic(:,n);
                item.sameStateKineticRms_J=sqrt(trapz(t,delta.^2)/(t(end)-t(1)));
                item.initialTotalDifference_J=r.energyDifference.total(1,n);
                item.minRcond=min(r.massRcond);
                item.maxCoordinateMagnitude_m=max(abs(r.X(:,ix)),[],'all');
                if row==1, summary=item; else, summary(row)=item; end
            end
            if g==1, c.runs=r; else, c.runs(g)=r; end
            fprintf('%s, xi=%.1f: max q %.2f mm, tips [%.1f %.1f %.1f] mm, min rcond %.2e\n', ...
                c.name,r.xi,1000*max(abs(r.coordinateDifference),[],'all'), ...
                1000*max(r.tipDifference,[],1),min(r.massRcond));
        end
        if caseIndex==1, sweep.cases=c; else, sweep.cases(caseIndex)=c; end
        sweep.summary=summary;
        save(fullfile(outdir,'cog_sweep.mat'),'sweep','-v7.3');
    end
end
fid=fopen(fullfile(outdir,'summary.json'),'w');
fwrite(fid,jsonencode(summary,'PrettyPrint',true),'char'); fclose(fid);
fprintf('All 40 trajectories completed. Data: %s\n',outdir);
end

function data=observablesSeries(X,p,model)
nt=size(X,1); data.tip=zeros(nt,3,3);
for f={'kinetic','gravity','elastic','total'}, data.(f{1})=zeros(nt,3); end
for k=1:nt
    obs=arm_module_observables(X(k,:).',p,model);
    data.tip(k,:,:)=reshape(obs.tip,1,3,3);
    for f={'kinetic','gravity','elastic','total'}, data.(f{1})(k,:)=obs.(f{1}); end
end
end

function verifyObservables(p)
q=[-.01;-.01;-.008;.004;.006;-.003]; dq=[.01;-.02;.015;.003;-.005;.008];
for model={'standard','pointmass'}
    for xi=[.1,.5,1]
        p.cog_xi=xi*ones(3,1);
        obs=arm_module_observables([q;dq],p,model{1});
        for n=1:3
            pn=p; pn.mi=zeros(3,1); pn.mi(n)=p.mi(n);
            if strcmp(model{1},'standard')
                M=armS_standard_core(q,dq,pn);
            else
                M=armS_core_N3_mex(0,reshape(q,2,3).',reshape(dq,2,3).', ...
                    pn.L,pn.r,pn.cog_xi,pn.mi,pn.g,pn.K,pn.beta);
            end
            expected=.5*dq.'*M*dq;
            assert(abs(obs.kinetic(n)-expected)<1e-10*max(1,abs(expected)), ...
                'Module kinetic energy does not match the production mass contribution.');
        end
    end
end
fprintf('Module energy formulas agree with isolated production mass contributions.\n');
end
