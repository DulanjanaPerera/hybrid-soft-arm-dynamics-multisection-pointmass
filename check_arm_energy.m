function results=check_arm_energy(pointX,beta)
% Zero-input mechanical-energy check for both translational-mass models.
s=load('comparison_standard_distributed.mat','t','X','params');
p=load('comparison_after_C_correction.mat','t','X','params');
if nargin==0
    p.params.beta=ones(3,3);
    outputFile='arm_energy_check.mat';
else
    assert(nargin==2 && isequal(size(pointX),size(p.X)) && ...
        isequal(size(beta),[3,3]));
    p.X=pointX;
    p.params.beta=beta;
    outputFile='arm_energy_check_candidate.mat';
end
assert(isequal(s.t,p.t) && all(p.params.cog_xi==0.5));
assert(all(p.params.tau==0) && all(s.params.tau==0));
t=s.t(:); N=numel(t);
results.time=t;
results.standard=energySeries(s.X,s.params,false);
results.pointmass=energySeries(p.X,p.params,true);
q=s.X(1,1:6).'; dq=s.X(1,7:12).'; step=1e-7;
results.gradientErrors=zeros(2,2);
for model=1:2
    pointmass=model==2;
    if pointmass
        prm=p.params;
        [~,~,G]=armS_core_N3_mex(0,reshape(q,2,3).', ...
            reshape(dq,2,3).',prm.L,prm.r,prm.cog_xi,prm.mi,prm.g,prm.K,prm.beta);
    else
        prm=s.params;
        [~,~,G]=armS_standard_core(q,dq,prm);
    end
    gf=zeros(6,1); kf=zeros(6,1);
    for h=1:6
        e=zeros(6,1); e(h)=step;
        gf(h)=(gravityPotential(q+e,prm,pointmass)- ...
            gravityPotential(q-e,prm,pointmass))/(2*step);
        kf(h)=(elasticPotential(q+e,prm)-elasticPotential(q-e,prm))/(2*step);
    end
    Kq=prm.K*q+0.5*prm.lKbounds(3)*(2+ ...
        tanh(prm.mu*(q-prm.lKbounds(2)))- ...
        tanh(prm.mu*(q-prm.lKbounds(1)))).*q;
    results.gradientErrors(model,:)=[norm(gf-G)/max(norm(G),1), ...
        norm(kf-Kq)/max(norm(Kq),1)];
end
fprintf('Potential-gradient relative errors [gravity elastic]: standard %.3e %.3e; pointmass %.3e %.3e\n', ...
    results.gradientErrors(1,:),results.gradientErrors(2,:));
for name={'standard','pointmass'}
    a=results.(name{1});
    balance=a.total-a.total(1)+cumtrapz(t,a.dissipation);
    a.balanceResidual=balance;
    a.maxIncrease=max(diff(a.total));
    a.maxBalanceResidual=max(abs(balance));
    a.energyDrop=a.total(1)-a.total(end);
    results.(name{1})=a;
    fprintf('%s: energy drop %.6g J, max one-sample increase %.3e J, max balance residual %.3e J (%.3g of drop)\n', ...
        name{1},a.energyDrop,a.maxIncrease,a.maxBalanceResidual, ...
        a.maxBalanceResidual/max(a.energyDrop,eps));
end
save(outputFile,'results');
end

function out=energySeries(X,p,pointmass)
N=size(X,1); out.kinetic=zeros(N,1); out.gravity=zeros(N,1);
out.elastic=zeros(N,1); out.dissipation=zeros(N,1);
for k=1:N
    q=X(k,1:6).'; dq=X(k,7:12).';
    if pointmass
        [M,~,~]=armS_core_N3_mex(0,reshape(q,2,3).', ...
            reshape(dq,2,3).',p.L,p.r,p.cog_xi,p.mi,p.g,p.K,p.beta);
    else
        M=armS_standard_core(q,dq,p);
    end
    out.kinetic(k)=0.5*dq.'*M*dq;
    out.gravity(k)=gravityPotential(q,p,pointmass);
    out.elastic(k)=elasticPotential(q,p);
    out.dissipation(k)=dq.'*p.D*dq;
end
out.total=out.kinetic+out.gravity+out.elastic;
end

function U=gravityPotential(q,p,pointmass)
R=eye(3); P=zeros(3,1); U=0;
for n=1:3
    l=[0,q(2*n-1:2*n).'];
    if pointmass
        [~,~,center]=HTM_nume(l,p.cog_xi(n),p.L,p.r);
    else
        center=integratedPosition_nume(l,p.L,p.r);
    end
    U=U+p.mi(n)*p.g.'*(P+R*center);
    [~,Rt,pt]=HTM_nume(l,1,p.L,p.r);
    P=P+R*pt; R=R*Rt;
end
end

function U=elasticPotential(q,p)
% Integrate the actual nonlinear force K_i(q_i)*q_i from 0 to q_i.
[x,w]=gauss32(); U=0;
for i=1:6
    z=q(i)*x;
    Ki=p.K(i,i)+0.5*p.lKbounds(3)*(2+ ...
        tanh(p.mu*(z-p.lKbounds(2)))- ...
        tanh(p.mu*(z-p.lKbounds(1))));
    U=U+q(i)*sum(w.*Ki.*z);
end
end

function [x,w]=gauss32()
k=(1:31).'; b=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1));
x=(diag(D)+1)/2; w=(V(1,:).^2).';
end
