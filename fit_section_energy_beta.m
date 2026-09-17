function result=fit_section_energy_beta()
% Fit section translational kinetic energy at identical q,dq, for each xi.
% The first section has no upstream B: only its beta_v3 is identifiable.
f=load('shared_beta_fit_bound5.mat','result'); old=f.result;
base=load('comparison_after_C_correction.mat','params'); p=base.params;
outdir=fullfile(pwd,'cog_beta_energy_fit'); if ~isfolder(outdir), mkdir(outdir); end
result.created=char(datetime('now')); result.locations=(1:10)/10;
result.trainPoses=old.trainPoses; result.heldPoses=old.heldPoses;
result.fixedBeta=old.betaMatrix; result.parameters=p;
result.velocitySeed=19017; result.velocityScale=0.06;
result.trainVelocities=makeVelocities(size(old.trainPoses,2),19017);
result.heldVelocities=makeVelocities(size(old.heldPoses,2),29017);
result.bounds=[1 100;0 100;1 100];
for gi=1:10
    xi=result.locations(gi); pp=p; pp.cog_xi=xi*ones(3,1);
    train=makeData(result.trainPoses,result.trainVelocities,pp);
    held=makeData(result.heldPoses,result.heldVelocities,pp);
    scales=max(sqrt(mean(train.target.^2,1)), ...
        .05*max(sqrt(mean(train.target.^2,1))));
    result.scales_J(gi,:)=scales;
    ranks=zeros(1,3); singular=zeros(3,3);
    for n=1:3
        A=train.features{n}./scales(n);
        v=svd(A,0); singular(1:numel(v),n)=v;
        ranks(n)=sum(v>max(size(A))*eps(max(v)));
    end
    assert(isequal(ranks,[1 3 3]),'Energy design is rank deficient at xi %.1f',xi);
    result.designRank(gi,:)=ranks;
    result.designSingular(:,:,gi)=singular;
    candidate=zeros(3,3,4); candidate(:,:,1)=ones(3,3);
    candidate(:,:,2)=old.betaMatrix;
    candidate(:,:,3)=fitShared(train,scales,old.betaShared);
    candidate(:,:,4)=fitSections(train,scales,old.betaMatrix);
    result.beta(:,:,:,gi)=candidate;
    for ci=1:4
        b=candidate(:,:,ci);
        result.train(gi,ci)=score(train,b,scales);
        result.held(gi,ci)=score(held,b,scales);
        result.condition(gi,ci)=checkMass([result.trainPoses,result.heldPoses],pp,b);
    end
    fprintf('xi %.1f: held normalized kinetic RMS [ones,fixed,shared,section] %.4g %.4g %.4g %.4g; ranks [%d %d %d]\n', ...
        xi,[result.held(gi,:).aggregateRms],ranks);
    save(fullfile(outdir,'fit.mat'),'result','-v7.3');
end
end

function V=makeVelocities(n,seed)
state=rng; cleanup=onCleanup(@() rng(state)); %#ok<NASGU>
rng(seed); V=.06*(2*rand(6,10,n)-1);
% Include pure local and mixed upstream motion for identifiable directions.
for j=1:6, V(:,j,:)=0; V(j,j,:)=.06; end
end

function data=makeData(poses,V,p)
np=size(poses,2); nv=size(V,2); ns=np*nv;
data.target=zeros(ns,3); data.base=zeros(ns,3);
data.features={zeros(ns,3),zeros(ns,3),zeros(ns,3)};
row=0;
for k=1:np
    q=poses(:,k); mats=zeros(6,6,3,5); ref=zeros(6,6,3);
    for n=1:3
        pn=p; pn.mi=zeros(3,1); pn.mi(n)=p.mi(n);
        ref(:,:,n)=armS_standard_core(q,zeros(6,1),pn);
        mats(:,:,n,1)=pointM(q,pn,ones(3,3));
        for j=1:3
            beta=ones(3,3); beta(n,j)=2;
            mats(:,:,n,j+1)=pointM(q,pn,beta)-mats(:,:,n,1);
        end
    end
    for j=1:nv
        row=row+1; v=V(:,j,k);
        for n=1:3
            data.target(row,n)=.5*v.'*ref(:,:,n)*v;
            data.base(row,n)=.5*v.'*mats(:,:,n,1)*v;
            for h=1:3
                data.features{n}(row,h)=.5*v.'*mats(:,:,n,h+1)*v;
            end
        end
    end
end
end

function M=pointM(q,p,beta)
M=armS_core_N3_mex(0,reshape(q,2,3).',zeros(3,2), ...
    p.L,p.r,p.cog_xi,p.mi,p.g,p.K,beta);
end

function B=fitShared(data,scales,initial)
A=[]; y=[];
for n=1:3
    A=[A;data.features{n}/scales(n)]; %#ok<AGROW>
    y=[y;(data.target(:,n)-data.base(:,n))/scales(n)]; %#ok<AGROW>
end
b=solveCone(A,y,initial(:)); B=repmat(b(:).',3,1);
end

function B=fitSections(data,scales,initial)
B=ones(3,3);
for n=1:3
    A=data.features{n}/scales(n);
    y=(data.target(:,n)-data.base(:,n))/scales(n);
    if n==1
        d=max(0,min(99,(A(:,3).'*y)/(A(:,3).'*A(:,3))));
        B(n,3)=1+d;
    else
        B(n,:)=solveCone(A,y,initial(n,:).').';
    end
end
end

function b=solveCone(A,y,initial)
lb=[1;0;1]; ub=[100;100;100];
obj=@(x) sum((A*(x-1)-y).^2)/size(A,1);
opts=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxIterations',400,'OptimalityTolerance',1e-12,'StepTolerance',1e-12);
starts=[initial(:),[1.4;1.7;2.5],[2;2;3],[1.1;1.1;1.2]];
best=inf; b=initial(:);
for k=1:size(starts,2)
    [x,val,flag]=fmincon(obj,starts(:,k),[],[],[],[],lb,ub,@cone,opts);
    if flag>0 && val<best, b=x; best=val; end
end
assert(isfinite(best),'Energy fit failed.');
end

function [c,ceq]=cone(b)
c=(b(2)-1)^2-(b(1)-1)*(b(3)-1); ceq=[];
end

function s=score(data,B,scales)
e=zeros(size(data.target));
for n=1:3
    e(:,n)=data.base(:,n)+data.features{n}*(B(n,:).'-1)-data.target(:,n);
end
s.rms_J=sqrt(mean(e.^2,1));
s.normalizedRms=s.rms_J./scales;
s.aggregateRms=sqrt(mean((e./scales).^2,'all'));
s.targetRms_J=sqrt(mean(data.target.^2,1));
end

function c=checkMass(poses,p,beta)
c.minRcond=inf; c.maxSymmetry=0; c.cholFlag=0;
for k=1:size(poses,2)
    M=pointM(poses(:,k),p,beta);
    c.minRcond=min(c.minRcond,rcond(M));
    c.maxSymmetry=max(c.maxSymmetry,norm(M-M.','fro'));
    [~,flag]=chol((M+M.')/2); c.cholFlag=max(c.cholFlag,flag);
end
assert(c.cholFlag==0,'Mass matrix not SPD at a sampled pose.');
end
