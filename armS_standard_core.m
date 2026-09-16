function [M,C,G,dM] = armS_standard_core(q,dq,params)
% Distributed uniform section mass; translational kinetic energy only.
% q,dq: [l12;l13;l22;l23;l32;l33] and their velocities (six entries).
% Integrals use normalized material coordinate xi in [0,1], dm=mi*dXi.
% G uses the SAME potential-gradient/sign convention as armS_core_N3_mex.
% No CoG location or empirical beta coefficients enter this model.
% R is treated as a rotation for the local-local block (R'*R=I).
% With polynomial HTMs this carries their rotation-orthogonality error.
nq = 6;
assert(params.N==3 && numel(q)==nq && numel(dq)==nq);
q=q(:); dq=dq(:);
M=zeros(nq); dM=zeros(nq,nq,nq); G=zeros(nq,1);
R=eye(3);
% Derivatives of the section-base position and orientation in the base frame.
Pq=zeros(3,nq); Rq=zeros(3,3,nq);
Pqq=zeros(3,nq,nq); Rqq=zeros(3,3,nq,nq);
for n=1:3
    old=1:2*(n-1); cur=2*n-1:2*n;
    l=[0,q(cur).'];
    [mu,muq,muqq]=integratedPosition_nume(l,params.L,params.r);
    [S,Sq]=integratedPositionProduct_compact(l,params.L,params.r);
    [F,Fq]=integratedPositionDerivativeProduct_compact(l,params.L,params.r);
    [E,Eq]=integratedJacobianProduct_compact(l,params.L,params.r);
    m=params.mi(n);
    % J_i(xi)=Pq_i+Rq_i*p for upstream coordinates;
    % J_a(xi)=R*p_q_a for this section's coordinates.
    for i=1:2*(n-1)
        A=Pq(:,i); B=Rq(:,:,i);
        G(i)=G(i)+m*(A+B*mu).'*params.g;
        for j=1:2*(n-1)
            D=Pq(:,j); H=Rq(:,:,j);
            M(i,j)=M(i,j)+m*(A.'*D+A.'*H*mu+D.'*B*mu+trace(B.'*H*S));
            for h=1:2*n
                Ah=Pqq(:,i,h); Bh=Rqq(:,:,i,h);
                Dh=Pqq(:,j,h); Hh=Rqq(:,:,j,h);
                muh=zeros(3,1); Sh=zeros(3);
                if h>2*(n-1)
                    a=h-2*(n-1); muh=muq(:,a); Sh=Sq(:,:,a);
                end
                dM(i,j,h)=dM(i,j,h)+m*(Ah.'*D+A.'*Dh ...
                    +Ah.'*H*mu+A.'*Hh*mu+A.'*H*muh ...
                    +Dh.'*B*mu+D.'*Bh*mu+D.'*B*muh ...
                    +trace((Bh.'*H+B.'*Hh)*S+B.'*H*Sh));
            end
        end
        for a=1:2
            j=cur(a);
            v=m*(A.'*R*muq(:,a)+trace(B.'*R*F(:,:,a)));
            M(i,j)=M(i,j)+v; M(j,i)=M(j,i)+v;
            for h=1:2*n
                uq_h=zeros(3,1); Fh=zeros(3);
                if h>2*(n-1)
                    b=h-2*(n-1); uq_h=muqq(:,a,b); Fh=Fq(:,:,a,b);
                end
                v=m*(Pqq(:,i,h).'*R*muq(:,a)+A.'*Rq(:,:,h)*muq(:,a) ...
                    +A.'*R*uq_h+trace((Rqq(:,:,i,h).'*R+B.'*Rq(:,:,h))*F(:,:,a) ...
                    +B.'*R*Fh));
                dM(i,j,h)=dM(i,j,h)+v; dM(j,i,h)=dM(j,i,h)+v;
            end
        end
    end
    M(cur,cur)=M(cur,cur)+m*E;
    for a=1:2
        dM(cur,cur,cur(a))=dM(cur,cur,cur(a))+m*Eq(:,:,a);
        G(cur(a))=G(cur(a))+m*(R*muq(:,a)).'*params.g;
    end
    % Advance base derivatives by differentiating Pnew=P+R*p, Rnew=R*Rt.
    [~,Rt,p]=HTM_nume(l,1,params.L,params.r);
    [pj,rj,pjj,rjj]=LocalJacob_nume(l,1,params.L,params.r);
    dp=zeros(3,nq); dr=zeros(3,3,nq);
    ddp=zeros(3,nq,nq); ddr=zeros(3,3,nq,nq);
    for a=1:2
        rows=3*(a-1)+(1:3);
        dp(:,cur(a))=pj(:,a); dr(:,:,cur(a))=rj(:,rows);
        for b=1:2
            cols=3*(b-1)+(1:3);
            ddp(:,cur(a),cur(b))=pjj(rows,b);
            ddr(:,:,cur(a),cur(b))=rjj(rows,cols);
        end
    end
    newPq=zeros(3,nq); newRq=zeros(3,3,nq);
    newPqq=zeros(3,nq,nq); newRqq=zeros(3,3,nq,nq);
    for i=1:2*n
        newPq(:,i)=Pq(:,i)+Rq(:,:,i)*p+R*dp(:,i);
        newRq(:,:,i)=Rq(:,:,i)*Rt+R*dr(:,:,i);
        for h=1:2*n
            newPqq(:,i,h)=Pqq(:,i,h)+Rqq(:,:,i,h)*p ...
                +Rq(:,:,i)*dp(:,h)+Rq(:,:,h)*dp(:,i)+R*ddp(:,i,h);
            newRqq(:,:,i,h)=Rqq(:,:,i,h)*Rt+Rq(:,:,i)*dr(:,:,h) ...
                +Rq(:,:,h)*dr(:,:,i)+R*ddr(:,:,i,h);
        end
    end
    R=R*Rt; Pq=newPq; Rq=newRq; Pqq=newPqq; Rqq=newRqq;
end
C=christoffelSymbol(3,dM,dq);
end
