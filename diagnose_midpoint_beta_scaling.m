function diagnosis=diagnose_midpoint_beta_scaling()
% Test whether scaling each midpoint section's entire mass contribution
% can match the distributed model at two bent initial configurations.
d=load('comparison_pointmass_standard_experiment.mat','results');
r=d.results; p=r.parameters;
if ~isfield(p,'beta'), p.beta=ones(3,3); end
ids=[1 3 4 6];
fprintf('Four bent trajectories (midpoint beta = 1):\n');
for j=ids
    c=r.cases(j);
    fprintf('%s: max q %.3f mm, RMS q %.3f mm, max tip %.1f mm, RMS tip %.1f mm\n', ...
        c.name,1000*c.maxCoordinate,1000*c.rmsCoordinate, ...
        1000*c.maxTip,1000*c.rmsTip);
end
qs=[r.cases(1).x0(1:6),r.cases(3).x0(1:6)];
A=[]; target=[]; blocks=cell(1,2); targets=cell(1,2);
for k=1:2
    q=qs(:,k); dq=zeros(6,1);
    Ms=armS_standard_core(q,dq,p);
    B=zeros(36,3); Bg=zeros(6,3);
    for n=1:3
        m=zeros(3,1); m(n)=p.mi(n);
        [Mi,~,Gi]=armS_core_N3_mex(0,reshape(q,2,3).', ...
            zeros(3,2),p.L,p.r,p.cog_xi,m,p.g,p.K,p.beta);
        B(:,n)=Mi(:); Bg(:,n)=Gi;
    end
    blocks{k}=struct('M',B,'G',Bg,'Ms',Ms);
    [~,~,Gs]=armS_standard_core(q,dq,p);
    blocks{k}.Gs=Gs;
    A=[A;B]; target=[target;Ms(:)]; %#ok<AGROW>
end
beta=lsqnonneg(A,target);
diagnosis.beta=beta;
fprintf('Best nonnegative per-section whole-mass scales fitted to both initial M matrices: [%.4f %.4f %.4f]\n',beta);
for k=1:2
    b=blocks{k};
    originalM=reshape(b.M*ones(3,1),6,6);
    fittedM=reshape(b.M*beta,6,6);
    diagnosis.relativeMassError(k,:)=[norm(originalM-b.Ms,'fro'), ...
        norm(fittedM-b.Ms,'fro')]/norm(b.Ms,'fro');
    diagnosis.relativeGravityError(k,:)=[norm(b.G*ones(3,1)-b.Gs), ...
        norm(b.G*beta-b.Gs)]/max(norm(b.Gs),1);
    diagnosis.rcond(k,:)=[rcond(b.Ms),rcond(originalM),rcond(fittedM)];
    fprintf('pose %d: M error %.3f -> %.3f, G error %.3f -> %.3f, rcond [std old fitted] %.2e %.2e %.2e\n', ...
        k,diagnosis.relativeMassError(k,:), ...
        diagnosis.relativeGravityError(k,:),diagnosis.rcond(k,:));
end
diagnosis.caseIds=ids;
diagnosis.cases=r.cases(ids);
save('midpoint_beta_scaling_diagnosis.mat','diagnosis');
figure('Visible','off','Color','w','Position',[100 100 1050 600]);
tiledlayout(2,1);
nexttile; hold on; grid on;
for j=ids
    c=r.cases(j);
    plot(r.time,1000*max(abs(c.standard(:,1:6)-c.pointmass(:,1:6)),[],2));
end
ylabel('Maximum coordinate difference (mm)');
nexttile; hold on; grid on;
for j=ids
    c=r.cases(j);
    plot(r.time,1000*vecnorm(c.tipStandard-c.tipPointmass,2,2));
end
xlabel('Time (s)'); ylabel('Tip position difference (mm)');
legend({r.cases(ids).name},'Location','eastoutside');
exportgraphics(gcf,'comparison_four_bent_cases.png','Resolution',150);
close(gcf);
end
