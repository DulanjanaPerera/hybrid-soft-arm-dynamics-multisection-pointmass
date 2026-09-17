function out=arm_module_observables(X,p,model)
% Physical section energies, including motion induced by upstream sections.
% model is 'standard' or 'pointmass'. Potential sign matches the model G.
q=X(1:6); dq=X(7:12); q=q(:); dq=dq(:);
out.tip=zeros(3,3); out.kinetic=zeros(1,3); out.gravity=zeros(1,3);
R=eye(3); origin=zeros(3,1); V=zeros(3,1); Rd=zeros(3);
vbody=zeros(3,1); omega=zeros(3);
isPoint=strcmp(model,'pointmass');
for n=1:3
    ix=2*n-1:2*n; l=[0,q(ix).']; u=dq(ix);
    [~,Rt,pt]=HTM_nume(l,1,p.L,p.r);
    [Pt,Rtq]=LocalJacob_nume(l,1,p.L,p.r);
    Rtdot=Rtq(:,1:3)*u(1)+Rtq(:,4:6)*u(2);
    if isPoint
        [~,~,center]=HTM_nume_mex(l,p.cog_xi(n),p.L,p.r);
        P=LocalJacob_nume_mex(l,p.cog_xi(n),p.L,p.r);
        b=p.beta(n,:); B=omega*center; local=P*u;
        out.kinetic(n)=0.5*p.mi(n)*(vbody.'*vbody+ ...
            2*vbody.'*B+2*vbody.'*local+b(1)*(B.'*B)+ ...
            2*b(2)*(B.'*local)+b(3)*(local.'*local));
    else
        [center,P]=integratedPosition_nume(l,p.L,p.r);
        S=integratedPositionProduct_compact(l,p.L,p.r);
        F=integratedPositionDerivativeProduct_compact(l,p.L,p.r);
        E=integratedJacobianProduct_compact(l,p.L,p.r);
        Fdot=F(:,:,1)*u(1)+F(:,:,2)*u(2);
        % Exactly the implemented standard M, including its local R'*R=I
        % simplification. Do not silently replace it with different energy.
        out.kinetic(n)=0.5*p.mi(n)*(V.'*V+2*V.'*Rd*center+ ...
            2*V.'*R*P*u+trace(Rd.'*Rd*S)+ ...
            2*trace(Rd.'*R*Fdot)+u.'*E*u);
    end
    out.gravity(n)=p.mi(n)*p.g(:).'*(origin+R*center);
    V=V+Rd*pt+R*Pt*u;
    Rd=Rd*Rt+R*Rtdot;
    vbody=Rt.'*(vbody+omega*pt+Pt*u);
    omega=Rt.'*(omega*Rt+Rtdot);
    origin=origin+R*pt; R=R*Rt;
    out.tip(:,n)=origin;
end
out.elastic=sectionElastic(q,p);
out.total=out.kinetic+out.gravity+out.elastic;
end

function U=sectionElastic(q,p)
persistent nodes weights
if isempty(nodes)
    k=(1:63).'; b=k./sqrt(4*k.^2-1);
    [V,D]=eig(diag(b,1)+diag(b,-1));
    nodes=(diag(D)+1)/2; weights=(V(1,:).^2).';
end
z=nodes*q.';
K=diag(p.K).'+0.5*p.lKbounds(3)*(2+ ...
    tanh(p.mu*(z-p.lKbounds(2)))-tanh(p.mu*(z-p.lKbounds(1))));
perCoordinate=q.'.*sum(weights.*K.*z,1);
% Resolve the narrow length-limit transition adaptively when crossed.
for i=1:6
    if q(i)<p.lKbounds(1)+.003 || q(i)>p.lKbounds(2)-.003
        force=@(x) (p.K(i,i)+0.5*p.lKbounds(3)*(2+ ...
            tanh(p.mu*(x-p.lKbounds(2)))- ...
            tanh(p.mu*(x-p.lKbounds(1))))).*x;
        perCoordinate(i)=integral(force,0,q(i),'AbsTol',1e-11,'RelTol',1e-10);
    end
end
U=sum(reshape(perCoordinate,2,3),1);
end
