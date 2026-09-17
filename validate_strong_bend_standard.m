function results=validate_strong_bend_standard()
% Probe the existing Taylor kinematics beyond the original validation poses.
s=load('comparison_after_C_correction.mat','params'); p=s.params;
poses=[zeros(6,1), repmat([.01;-.01],3,1), ...
    repmat([.02;-.02],3,1), repmat([.03;-.03],3,1), ...
    repmat([.04;-.04],3,1), [.03;0;-.025;.015;.02;-.03]];
dq=[.01;-.02;.015;.003;-.005;.008];
results=zeros(size(poses,2),7);
for k=1:size(poses,2)
    q=poses(:,k);
    [M,C,G,dM]=armS_standard_core(q,dq,p);
    [Mq,Gq,orth]=directReference(q,p);
    fd=zeros(6,6,6);
    for h=1:6
        e=zeros(6,1); e(h)=1e-7;
        fd(:,:,h)=(armS_standard_core(q+e,dq,p)- ...
            armS_standard_core(q-e,dq,p))/(2e-7);
    end
    Mdot=zeros(6);
    for h=1:6, Mdot=Mdot+dM(:,:,h)*dq(h); end
    z=Mdot-2*C;
    results(k,:)=[max(abs(q)),orth, ...
        norm(M-Mq,'fro')/norm(Mq,'fro'), ...
        norm(G-Gq)/max(norm(Gq),1), ...
        norm(dM(:)-fd(:))/max(norm(fd(:)),eps), ...
        norm(z+z.','fro')/max(norm(Mdot,'fro')+2*norm(C,'fro'),eps), ...
        rcond(M)];
    fprintf('pose %d max q %.0f mm: orth %.3e, M quad %.3e, G quad %.3e, dM %.3e, skew %.3e, rcond %.3e\n', ...
        k,1000*results(k,1),results(k,2:end));
end
save('strong_bend_validation.mat','poses','results');
end

function [M,G,orth]=directReference(q,p)
k=(1:15).'; b=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1));
xi=(diag(D)+1)/2; w=(V(1,:).^2).';
R=eye(3); Pq=zeros(3,6); Rq=zeros(3,3,6);
M=zeros(6); G=zeros(6,1); orth=0;
for n=1:3
    old=1:2*(n-1); cur=2*n-1:2*n; l=[0,q(cur).'];
    for j=1:16
        [~,~,pos]=HTM_nume(l,xi(j),p.L,p.r);
        pj=LocalJacob_nume(l,xi(j),p.L,p.r);
        J=zeros(3,6);
        for a=old, J(:,a)=Pq(:,a)+Rq(:,:,a)*pos; end
        J(:,cur)=R*pj;
        M=M+p.mi(n)*w(j)*(J.'*J);
        G=G+p.mi(n)*w(j)*J.'*p.g;
    end
    [~,Rt,pos]=HTM_nume(l,1,p.L,p.r);
    orth=max(orth,norm(Rt.'*Rt-eye(3),'fro'));
    [pj,rj]=LocalJacob_nume(l,1,p.L,p.r);
    for a=old
        Pq(:,a)=Pq(:,a)+Rq(:,:,a)*pos;
        Rq(:,:,a)=Rq(:,:,a)*Rt;
    end
    Pq(:,cur)=R*pj;
    for a=1:2, Rq(:,:,cur(a))=R*rj(:,3*(a-1)+(1:3)); end
    R=R*Rt;
end
end
