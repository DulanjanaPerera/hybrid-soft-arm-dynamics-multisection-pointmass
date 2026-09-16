function validate_armS_standard(params)
% Check derivative indexing, Christoffel contraction, and spatial integration.
% Run validate_armS_standard(params) before building a MEX.
poses=[zeros(6,1),[-.001;-.001;-.001;-.001;-1e-6;-1e-6], ...
       [-.008;.003;-.004;-.006;.002;-.005]];
dq=[.01;-.02;.015;.003;-.005;.008];
step=1e-7;
for k=1:size(poses,2)
    q=poses(:,k);
    [M,C,G,dM]=armS_standard_core(q,dq,params);
    fd=zeros(6,6,6);
    for h=1:6
        delta=zeros(6,1); delta(h)=step;
        Mp=armS_standard_core(q+delta,dq,params);
        Mm=armS_standard_core(q-delta,dq,params);
        fd(:,:,h)=(Mp-Mm)/(2*step);
    end
    Mdot=zeros(6); for h=1:6, Mdot=Mdot+dM(:,:,h)*dq(h); end
    Z=Mdot-2*C;
    [Mq,Gq]=quadratureReference(q,params);
    symmetry=norm(M-M.','fro')/max(norm(M,'fro'),eps);
    derivative=norm(dM(:)-fd(:))/max(norm(fd(:)),eps);
    skew=norm(Z+Z.','fro')/max(norm(Mdot,'fro')+2*norm(C,'fro'),eps);
    massIntegral=norm(M-Mq,'fro')/max(norm(Mq,'fro'),eps);
    gravityIntegral=norm(G-Gq)/max(norm(Gq),1);
    [~,pd]=chol((M+M.')/2);
    fprintf('Pose %d: symmetry %.3e, dM %.3e, skew %.3e, integral M %.3e, G %.3e\n', ...
        k,symmetry,derivative,skew,massIntegral,gravityIntegral);
    assert(pd==0,'Mass matrix is not positive definite.');
    assert(symmetry<1e-10 && derivative<1e-5 && skew<1e-10, ...
        'Mass derivative / Christoffel checks failed.');
    assert(massIntegral<1e-5 && gravityIntegral<1e-8, ...
        'Integral check failed. Check integral expressions and HTM orthogonality.');
end
fprintf('All standard-model checks passed at the three test poses.\n');
end

function [M,G]=quadratureReference(q,params)
% Independent assembly directly from global point Jacobians, 16 Gauss nodes.
k=(1:15).'; b=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1));
xi=(diag(D)+1)/2; w=(V(1,:).^2).';
R=eye(3); Pq=zeros(3,6); Rq=zeros(3,3,6);
M=zeros(6); G=zeros(6,1);
for n=1:3
    old=1:2*(n-1); cur=2*n-1:2*n; l=[0,q(cur).'];
    for s=1:16
        [~,~,p]=HTM_nume(l,xi(s),params.L,params.r);
        pj=LocalJacob_nume(l,xi(s),params.L,params.r);
        J=zeros(3,6);
        for a=old, J(:,a)=Pq(:,a)+Rq(:,:,a)*p; end
        J(:,cur)=R*pj;
        M=M+params.mi(n)*w(s)*(J.'*J);
        G=G+params.mi(n)*w(s)*J.'*params.g;
    end
    [~,Rt,p]=HTM_nume(l,1,params.L,params.r);
    [pj,rj]=LocalJacob_nume(l,1,params.L,params.r);
    for a=old
        Pq(:,a)=Pq(:,a)+Rq(:,:,a)*p;
        Rq(:,:,a)=Rq(:,:,a)*Rt;
    end
    Pq(:,cur)=R*pj;
    for a=1:2, Rq(:,:,cur(a))=R*rj(:,3*(a-1)+(1:3)); end
    R=R*Rt;
end
end
