function results = comparePointmassStandardExperiment()
% Matched dynamics experiment; point mass stays at each section midpoint.
ref = load('comparison_after_C_correction.mat','t','X','params');
p = ref.params;
p.cog_xi = 0.5*ones(3,1);
assert(p.N==3 && exist('armS_standard_mex','file')==3 && ...
    exist('armS_dynamics_N3_entry_mex_mex','file')==3);
t = ref.t(:);
q0 = ref.X(1,1:6).';
dq0 = ref.X(1,7:12).';
starts = [q0, zeros(6,1), [-.008;.003;-.004;-.006;.002;-.005]];
names = {'reference bend','straight','asymmetric bend'};
tau0 = p.tau(:);
% A modest, constant differential input on the first section.
inputs = [tau0, tau0 + [0.5;-0.5;0;0;0;0]];
inputNames = {'free response','differential input'};
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
results.parameters = p;
results.time = t;
results.inputNames = inputNames;
results.startNames = names;
results.cases = struct([]);
caseNo = 0;
for u=1:2
    for k=1:3
        caseNo=caseNo+1;
        x0=[starts(:,k);dq0];
        tau=inputs(:,u);
        std=@(tt,x) armS_standard_mex(tt,x,p.L,p.r,p.mi(:),p.g(:), ...
            p.K,p.D,tau,p.mu,p.lKbounds(:));
        pm=@(tt,x) armS_dynamics_N3_entry_mex_mex(tt,x,p.L,p.r, ...
            p.cog_xi(:),p.mi(:),p.g(:),p.K,p.D,tau,p.mu,p.lKbounds(:));
        [ts,xs]=ode15s(std,t,x0,opts);
        [tp,xp]=ode15s(pm,t,x0,opts);
        assert(isequal(ts,tp) && all(isfinite(xs(:))) && all(isfinite(xp(:))));
        dq=xs(:,1:6)-xp(:,1:6);
        dv=xs(:,7:12)-xp(:,7:12);
        tips=zeros(numel(t),3,2);
        for j=1:numel(t)
            tips(j,:,1)=tipPosition(xs(j,1:6).',p).';
            tips(j,:,2)=tipPosition(xp(j,1:6).',p).';
        end
        tipError=vecnorm(tips(:,:,1)-tips(:,:,2),2,2);
        c.name=sprintf('%s / %s',names{k},inputNames{u});
        c.x0=x0; c.tau=tau; c.standard=xs; c.pointmass=xp;
        c.tipStandard=tips(:,:,1); c.tipPointmass=tips(:,:,2);
        c.maxCoordinate=max(abs(dq),[],'all');
        c.rmsCoordinate=sqrt(mean(dq(:).^2));
        c.maxVelocity=max(abs(dv),[],'all');
        c.maxTip=max(tipError);
        c.rmsTip=sqrt(mean(tipError.^2));
        c.finalTip=tipError(end);
        if caseNo==1
            results.cases=c;
        else
            results.cases(caseNo)=c;
        end
        fprintf('%s: max |q| %.3e m, RMS |q| %.3e m, max tip %.3e m, RMS tip %.3e m\n', ...
            c.name,c.maxCoordinate,c.rmsCoordinate,c.maxTip,c.rmsTip);
    end
end
% Instantaneous comparisons at identical states, independent of drift.
states=[starts, results.cases(1).standard(round(end/2),1:6).', ...
    results.cases(1).pointmass(round(end/2),1:6).'];
static=zeros(size(states,2),4);
conditioning=zeros(size(states,2),2);
accelerationNorms=zeros(size(states,2),3);
for k=1:size(states,2)
    q=states(:,k); dq=dq0;
    [Ms,Cs,Gs]=armS_standard_core(q,dq,p);
    [Mp,Cp,Gp]=armS_core_N3_mex(0,reshape(q,2,3).', ...
        reshape(dq,2,3).',p.L,p.r,p.cog_xi,p.mi,p.g,p.K);
    static(k,1)=norm(Ms-Mp,'fro')/norm(Ms,'fro');
    conditioning(k,:)=[rcond(Ms),rcond(Mp)];
    static(k,2)=norm(Gs-Gp)/max(norm(Gs),eps);
    static(k,3)=norm(Cs-Cp,'fro')/max(norm(Cs,'fro'),eps);
    K=p.K+diag(0.5*p.lKbounds(3)*(2+ ...
        tanh(p.mu*(q-p.lKbounds(2)))-tanh(p.mu*(q-p.lKbounds(1)))));
    as=Ms\(tau0-(Cs+p.D)*dq-Gs-K*q);
    ap=Mp\(tau0-(Cp+p.D)*dq-Gp-K*q);
    static(k,4)=norm(as-ap)/max(norm(as),eps);
    accelerationNorms(k,:)=[norm(as),norm(ap),norm(as-ap)];
end
results.staticRelativeErrors=static;
results.massMatrixRcond=conditioning;
results.accelerationNorms=accelerationNorms;
disp('Static columns: relative M, G, C, and acceleration differences');
disp(static);
save('comparison_pointmass_standard_experiment.mat','results');
figure('Visible','off','Color','w','Position',[100 100 1100 650]);
tiledlayout(2,1);
nexttile; hold on; grid on;
for k=1:numel(results.cases)
    plot(t,max(abs(results.cases(k).standard(:,1:6)- ...
        results.cases(k).pointmass(:,1:6)),[],2));
end
ylabel('Maximum coordinate error (m)');
nexttile; hold on; grid on;
for k=1:numel(results.cases)
    plot(t,vecnorm(results.cases(k).tipStandard- ...
        results.cases(k).tipPointmass,2,2));
end
xlabel('Time (s)'); ylabel('Tip position error (m)');
legend({results.cases.name},'Location','eastoutside');
exportgraphics(gcf,'comparison_pointmass_standard_errors.png','Resolution',150);
close(gcf);
end

function tip=tipPosition(q,p)
R=eye(3); tip=zeros(3,1);
for n=1:3
    [~,Rn,pn]=HTM_nume([0,q(2*n-1:2*n).'],1,p.L,p.r);
    tip=tip+R*pn;
    R=R*Rn;
end
end
