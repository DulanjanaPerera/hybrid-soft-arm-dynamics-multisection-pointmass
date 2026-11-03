function acInverseCheck(root,path,folder)
%ACINVERSECHECK Offline self-consistency only, NOT a control/physical force test.
base=fullfile(root,'residual_dynamics_learning','force_model');
addpath(base,fullfile(base,'streaming'),fullfile(base,'fine_tuning'));
hash=sgHash(path);a=load(path);m=a.model;cfg=m.cfg;
D=load(fullfile(base,'fine_tuning','campaign_20261103_v1','data','data.mat'));H=sgStates(m.params,D,D.index,cfg);
assert(~isfolder(folder));mkdir(folder);
protocol=struct('modelPath',path,'modelHash',hash,'selection','seed113 representative of each predeclared architecture; frozen before external evaluation', ...
    'meaning','self-consistent model inverse; no measured force or tracking claim','commonPressure','mean logged target pressure', ...
    'bounds_bar',[0 3],'boundsScope','offline RL command domain; not a PAC runtime setting', ...
    'slew_barps',6,'period_s',.05,'state','pure trial evaluation; final selected state committed once', ...
    'optimizer','fmincon SQP, 2 differential coordinates, common pressure fixed');
save(fullfile(folder,'protocol.mat'),'protocol');
A=pi*.0065^2*1e5*[-1 1 0;-1 0 1];E=[-1 -1;2 -1;-1 2]/3;
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',60, ...
    'OptimalityTolerance',1e-9,'ConstraintTolerance',1e-10,'StepTolerance',1e-11);
T=table;cases=struct([]);
% One untimed warm-up avoids conflating first-call/JIT cost with steady timing.
sgCell(m.params,zeros(8,1),sgZero(cfg,1),cfg);
for i=1:height(D.index)
    S=D.sequences{D.index.run(i)};k=D.index.anchor(i);pTarget=S.p(k,:)';pPrevious=S.p(k-1,:)';
    lb=max(0,pPrevious-.3);ub=min(3,pPrevious+.3);p0=mean(pTarget);
    feasible=all(pTarget>=lb-1e-10&pTarget<=ub+1e-10);
    if ~feasible
        T=[T;table(i,string(S.runId),false,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            'VariableNames',{'anchorId','runId','knownTargetFeasible','forceError_N','boundViolation_bar','solve_ms','exitflag','evaluations','minimumSingular_NperBar','pressurePerturbationGain'})];continue; %#ok<AGROW>
    end
    h=cellfun(@(z)z(:,i),H,'UniformOutput',false);before=h;
    raw=[S.q(k,:) S.v(k,:)]';dt=S.dt(k);
    target=force(pTarget,h,raw,dt,m,A);
    % Deliberately displaced initial guess avoids trivial solves during holds.
    u0=[pPrevious(2)-pPrevious(1);pPrevious(3)-pPrevious(1)]+[.15;-.10];
    cost=@(u)sum((force(p0+E*u,h,raw,dt,m,A)-target).^2);
    clock=tic;[u,~,flag,output]=fmincon(cost,u0,[E;-E],[ub-p0;p0-lb],[],[],[],[],[],opts);
    p=p0+E*u;[hCommit,r]=sgCell(m.params,([raw;p;dt]-m.normalizer.mean)./m.normalizer.std,h,cfg); %#ok<ASGLU>
    elapsed=1e3*toc(clock);err=norm(A*p+r-target);violation=max([0;lb-p;p-ub]);
    assert(isequal(h,before),'Trial solve mutated committed state.');
    [expected,~]=sgCell(m.params,([raw;p;dt]-m.normalizer.mean)./m.normalizer.std,before,cfg);
    assert(isequal(hCommit,expected),'Selected state must advance exactly once.');
    J=zeros(2);eps=1e-5;
    for d=1:2,step=zeros(2,1);step(d)=eps;J(:,d)=(force(p+E*step,h,raw,dt,m,A)-force(p-E*step,h,raw,dt,m,A))/(2*eps);end
    singular=min(svd(J));delta=[1;-1]*.001;
    pert=@(u)sum((force(p0+E*u,h,raw,dt,m,A)-target-delta).^2);
    u2=fmincon(pert,u,[E;-E],[ub-p0;p0-lb],[],[],[],[],[],opts);
    gain=norm(E*(u2-u))/norm(delta);
    T=[T;table(i,string(S.runId),true,err,violation,elapsed,flag,output.funcCount,singular,gain, ...
        'VariableNames',{'anchorId','runId','knownTargetFeasible','forceError_N','boundViolation_bar','solve_ms','exitflag','evaluations','minimumSingular_NperBar','pressurePerturbationGain'})]; %#ok<AGROW>
    cases=[cases;struct('anchorId',i,'pTarget',pTarget,'pPrevious',pPrevious,'pSolved',p,'target_N',target, ...
        'jacobian_NperBar',J,'perturbedPressure',p0+E*u2,'preState',{h},'committedState',{hCommit})]; %#ok<AGROW>
end
T.pressureRecoveryError_bar=NaN(height(T),1);T.jacobianCondition=NaN(height(T),1);
for j=1:numel(cases)
    k=cases(j).anchorId;T.pressureRecoveryError_bar(k)=norm(cases(j).pSolved-cases(j).pTarget);
    T.jacobianCondition(k)=cond(cases(j).jacobian_NperBar);
end
writetable(T,fullfile(folder,'inverse_cases.csv'));save(fullfile(folder,'inverse_cases.mat'),'cases','protocol','-v7.3');
good=T.knownTargetFeasible;summary=struct('cases',height(T),'knownFeasible',nnz(good), ...
    'forceErrorMedian_N',median(T.forceError_N(good)),'forceErrorMax_N',max(T.forceError_N(good)), ...
    'maxBoundViolation_bar',max(T.boundViolation_bar(good)), ...
    'latencyMedian_ms',median(T.solve_ms(good)),'latencyP95_ms',prctile(T.solve_ms(good),95), ...
    'fractionOver50ms',mean(T.solve_ms(good)>50),'statePurityPass',true);
save(fullfile(folder,'summary.mat'),'summary');writetable(struct2table(summary),fullfile(folder,'summary.csv'));assert(strcmp(sgHash(path),hash));disp(summary);
copyfile(mfilename('fullpath')+".m",fullfile(folder,'acInverseCheck.m'));
end
function y=force(p,h,raw,dt,m,A)
[~,r]=sgCell(m.params,([raw;p;dt]-m.normalizer.mean)./m.normalizer.std,h,m.cfg);y=A*p+r;
end
