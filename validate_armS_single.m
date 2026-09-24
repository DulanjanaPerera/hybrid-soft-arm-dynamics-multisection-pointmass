function validate_armS_single()
% Compare one-section dynamics with quadrature and the isolated first section.
poses=[zeros(2,1),[-.001;-.001],[-.008;.003]];
dq=[.01;-.02];
Rset=stationaryBaseTestRotations();
L=.278; r=.013; mi=.1; gWorld=[0;0;9.81];
K=2200*eye(2); D=600*eye(2); tau=[.1;-.1];
mu=2000; bounds=[-.02;.02;1e6];
params=struct('N',1,'L',L,'r',r,'mi',mi,'g',gWorld, ...
    'K',K,'D',D,'tau',tau,'mu',mu,'lKbounds',bounds);
worstQuadratureM=0; worstQuadratureG=0; worstThreeSection=0;
worstIdentity=0; worstDerivative=0; worstSkew=0; worstStatic=0;
for p=1:size(poses,2)
    q=poses(:,p); X=[q;dq];
    params.g=gWorld;
    original=armS_single_entry(0,X,L,r,mi,gWorld,K,D,tau,mu,bounds);
    [M0,C0,~,dM0]=armS_single_core(q,dq,params);
    fd=zeros(2,2,2); step=1e-7;
    for h=1:2
        delta=zeros(2,1); delta(h)=step;
        Mp=armS_single_core(q+delta,dq,params);
        Mm=armS_single_core(q-delta,dq,params);
        fd(:,:,h)=(Mp-Mm)/(2*step);
    end
    worstDerivative=max(worstDerivative, ...
        norm(dM0(:)-fd(:))/max(norm(fd(:)),eps));
    Mdot=dM0(:,:,1)*dq(1)+dM0(:,:,2)*dq(2);
    Z=Mdot-2*C0;
    worstSkew=max(worstSkew, ...
        norm(Z+Z.','fro')/max(norm(Mdot,'fro')+2*norm(C0,'fro'),eps));
    for k=1:numel(Rset)
        R=Rset{k}; params.g=R.'*gWorld;
        [M,C,G,dM]=armS_single_core(q,dq,params);
        [Mq,Gq]=oneSectionWorldQuadrature(q,L,r,mi,R,gWorld);
        [~,pd]=chol((M+M.')/2);
        assert(pd==0,'One-section mass matrix is not positive definite.');
        worstQuadratureM=max(worstQuadratureM, ...
            norm(M-Mq,'fro')/max(norm(Mq,'fro'),eps));
        worstQuadratureG=max(worstQuadratureG, ...
            norm(G-Gq)/max(norm(Gq),1));

        params3=params;
        params3.N=3; params3.mi=[mi;0;0];
        [M3,C3,G3,dM3]=armS_standard_core( ...
            [q;zeros(4,1)],[dq;zeros(4,1)],params3);
        match=max([norm(M-M3(1:2,1:2),'fro'), ...
            norm(C-C3(1:2,1:2),'fro'),norm(G-G3(1:2)), ...
            norm(dM(:)-reshape(dM3(1:2,1:2,1:2),[],1))]);
        worstThreeSection=max(worstThreeSection,match);
        if k==1
            rhs=armS_single_stationary_base_entry( ...
                0,X,L,r,mi,R,gWorld,K,D,tau,mu,bounds);
            worstIdentity=max(worstIdentity,max(abs(rhs-original)));
        end
        Keff=K;
        for i=1:2
            Keff(i,i)=K(i,i)+.5*bounds(3)*(2 ...
                +tanh(mu*(q(i)-bounds(2))) ...
                -tanh(mu*(q(i)-bounds(1))));
        end
        rhsStatic=armS_single_stationary_base_entry( ...
            0,[q;zeros(2,1)],L,r,mi,R,gWorld,K,D, ...
            G+Keff*q,mu,bounds);
        worstStatic=max(worstStatic,norm(rhsStatic));
    end
end
assert(worstQuadratureM<1e-5 && worstQuadratureG<1e-8 ...
    && worstThreeSection<1e-10 && worstIdentity<1e-12 ...
    && worstDerivative<1e-5 && worstSkew<1e-10 ...
    && worstStatic<1e-8,'One-section validation failed.');
fprintf(['PASS single module: quadrature M %.3e, G %.3e; ' ...
    'three-section match %.3e; identity RHS %.3e; ' ...
    'dM %.3e; skew %.3e; static RHS %.3e\n'], ...
    worstQuadratureM,worstQuadratureG,worstThreeSection, ...
    worstIdentity,worstDerivative,worstSkew,worstStatic);
end

function [M,G]=oneSectionWorldQuadrature(q,L,r,mi,R,gWorld)
% Independent 16-node Gauss integration of world point Jacobians.
k=(1:15).'; b=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1));
xi=(diag(D)+1)/2; w=(V(1,:).^2).';
M=zeros(2); G=zeros(2,1);
l=[0,q.'];
for s=1:16
    J=LocalJacob_nume(l,xi(s),L,r);
    Jworld=R*J;
    M=M+mi*w(s)*(Jworld.'*Jworld);
    G=G+mi*w(s)*Jworld.'*gWorld;
end
end
