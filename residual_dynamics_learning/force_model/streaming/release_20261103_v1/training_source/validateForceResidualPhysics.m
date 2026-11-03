function validateForceResidualPhysics(folder)
%VALIDATEFORCERESIDUALPHYSICS Existing-RHS, force-placement, gradients, solver.
assert(~isfolder(folder));mkdir(folder);cfg=forceResidualConfig();n=cfg.nominal;
rng(11);Q=[zeros(2,1),.02*(2*rand(2,40)-1),[.032 -.017;.001 -.016]];
V=.02*(2*rand(size(Q))-1);P=forceResidualPhysics(Q,V,n);err=zeros(1,size(Q,2));
for k=1:size(Q,2)
    [M,C,G]=armS_single_core(Q(:,k),V(:,k),n);
    reference=[M(1,1);M(1,2);M(2,2);C(1,1);C(1,2);C(2,1);C(2,2);G];
    err(k)=max(abs(P(:,k)-reference)./max(1,abs(reference)));
end
assert(max(err)<1e-9,'Differentiable physics differs from existing nominal.');
alpha=pi*.0065^2*1e5;q=[.004;-.003];v=[.001;-.002];p=[1;1;1];r=[.8;-.3];
[a,b]=forceResidualStep(q,v,p,r,.05,cfg);
changed=p+[0;r/alpha];[c,d]=forceResidualStep(q,v,changed,zeros(2,1),.05,cfg);
assert(max(abs([a-c;b-d]))<1e-12,'Residual not in generalized-force RHS.');
[~,g]=dlfeval(@forceGradient,dlarray(r),q,v,p,cfg);epsilon=1e-5;fd=zeros(2,1);
for k=1:2
    delta=zeros(2,1);delta(k)=epsilon;
    plus=forceResidualStep(q,v,p,r+delta,.05,cfg);minus=forceResidualStep(q,v,p,r-delta,.05,cfg);
    fd(k)=(sum(plus.^2)-sum(minus.^2))/(2*epsilon);
end
gradientError=max(abs(extractdata(g)-fd)./max(1e-9,abs(fd)));
assert(gradientError<1e-4,'Force-through-integration gradient failed.');
% Nonlinear RHS reference, including a case near the corrected implied-l1 bound.
cases={[zeros(2,1);zeros(2,1)],[.031;-.004;0;0],[-.016;-.016;0;0]};
pressures={[0;.9;.3],[0;2.5;.4],[2.8;0;0]};convergence=table;
for index=1:numel(cases)
    x=cases{index};p=pressures{index};params=n;params.tau=alpha*[-1 1 0;-1 0 1]*p;
    [~,ref]=ode15s(@(t,x)armS_single_dynamics(t,x,params),[0 .05],x, ...
        odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.0005));
    for sub=[1 2 4 8]
        test=cfg;test.substeps=sub;[q1,~]=forceResidualStep(x(1:2),x(3:4),p,[0;0],.05,test);
        error_mm=1000*norm(q1-ref(end,1:2)');
        convergence=[convergence;table(index,sub,error_mm,'VariableNames', ...
            {'caseId','substeps','positionErrorVsOde15s_mm'})]; %#ok<AGROW>
    end
end
writetable(convergence,fullfile(folder,'integration_convergence.csv'));
disp(convergence);
assert(all(convergence.positionErrorVsOde15s_mm(convergence.substeps==8)<.05));
writetable(table(max(err),gradientError,'VariableNames',{'maxNominalRelativeError','forceGradientRelativeError'}),fullfile(folder,'verification.csv'));
save(fullfile(folder,'verification.mat'),'cfg','err','gradientError','convergence');
fprintf('PASS: exact nominal coefficients (%g), force RHS sign/units, integration gradients (%g), nonlinear ODE convergence.\n',max(err),gradientError);
end
function [y,g]=forceGradient(r,q,v,p,cfg)
next=forceResidualStep(q,v,p,r,.05,cfg);y=sum(next.^2,'all');g=dlgradient(y,r);
end
