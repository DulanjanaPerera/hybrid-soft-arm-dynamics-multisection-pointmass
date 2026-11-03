function [q,v]=forceResidualStep(q,v,p,residual,dt,cfg)
%FORCERESIDUALSTEP Second-order L-stable SDIRK; nonlinear implicit stages.
% M(q)ddq+(C(q,v)+D)v+Kq+bound(q)=alpha*A*p+residual (N).
% Two-state batch solve/Newton is differentiable; no learned position jumps.
h=reshape(dt,1,[])/cfg.substeps;n=cfg.nominal;
force=pi*(.0065)^2*1e5*[-1 1 0;-1 0 1]*p+residual;
gamma=1-1/sqrt(2);
for sub=1:cfg.substeps
    [q1,v1]=stage(q,v,force,gamma*h,n,cfg);
    baseQ=q+(1-gamma)/gamma*(q1-q);baseV=v+(1-gamma)/gamma*(v1-v);
    [q,v]=stage(baseQ,baseV,force,gamma*h,n,cfg);
end
end
function [next,newV]=stage(q,v,force,h,n,cfg)
next=q;newV=v;
for outer=1:cfg.stageIterations
    P=forceResidualPhysics(next,newV,n);m11=P(1,:);m12=P(2,:);m22=P(3,:);
    d11=P(4,:)+n.D(1,1);d12=P(5,:)+n.D(1,2);
    d21=P(6,:)+n.D(2,1);d22=P(7,:)+n.D(2,2);
    a11=m11+h.*d11+h.^2*n.K(1,1);a12=m12+h.*d12+h.^2*n.K(1,2);
    a21=m12+h.*d21+h.^2*n.K(2,1);a22=m22+h.*d22+h.^2*n.K(2,2);
    rhs1=m11.*(q(1,:)+h.*v(1,:))+m12.*(q(2,:)+h.*v(2,:))+ ...
        h.*(d11.*q(1,:)+d12.*q(2,:))+h.^2.*(force(1,:)-P(8,:));
    rhs2=m12.*(q(1,:)+h.*v(1,:))+m22.*(q(2,:)+h.*v(2,:))+ ...
        h.*(d21.*q(1,:)+d22.*q(2,:))+h.^2.*(force(2,:)-P(9,:));
    for iteration=1:cfg.newtonIterations
        [b11,b12,b22,f1,f2]=bound(next,n);
        e1=a11.*next(1,:)+a12.*next(2,:)+h.^2.*f1-rhs1;
        e2=a21.*next(1,:)+a22.*next(2,:)+h.^2.*f2-rhs2;
        j11=a11+h.^2.*b11;j12=a12+h.^2.*b12;
        j21=a21+h.^2.*b12;j22=a22+h.^2.*b22;
        determinant=j11.*j22-j12.*j21;
        delta=[(j22.*e1-j12.*e2)./determinant;(-j21.*e1+j11.*e2)./determinant];
        scale=min(1,cfg.newtonMoveLimit_m./(max(abs(delta),[],1)+1e-30));
        next=next-delta.*scale;
    end
    newV=(next-q)./h;
end
end
function [j11,j12,j22,f1,f2]=bound(q,n)
l=[-sum(q,1);q];lo=n.lKbounds(1);hi=n.lKbounds(2);strength=n.lKbounds(3);
t1=tanh(n.mu*(l-hi));t2=tanh(n.mu*(l-lo));
c=.5*strength*(2+t1-t2);derivative=c+.5*strength*n.mu*l.*(t2.^2-t1.^2);
f=c.*l;f1=f(2,:)-f(1,:);f2=f(3,:)-f(1,:);
j11=derivative(1,:)+derivative(2,:);j12=derivative(1,:);j22=derivative(1,:)+derivative(3,:);
end
