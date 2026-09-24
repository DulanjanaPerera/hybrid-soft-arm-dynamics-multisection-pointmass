function validate_armS_stationary_base()
% Validate stationary orientation against the original RHS and world quadrature.
poses=[zeros(6,1),[-.001;-.001;-.001;-.001;-1e-6;-1e-6], ...
    [-.008;.003;-.004;-.006;.002;-.005]];
dq=[.01;-.02;.015;.003;-.005;.008];
Rset=stationaryBaseTestRotations();
L=.278; r=.013; mi=.1*ones(3,1); gWorld=[0;0;-9.81];
K=2200*eye(6); D=100*eye(6); tau=[.1;-.1;0;0;0;0];
mu=2000; bounds=[-.02;.02;1e6];
params=struct('N',3,'L',L,'r',r,'mi',mi,'g',gWorld, ...
    'K',K,'D',D,'tau',tau,'mu',mu,'lKbounds',bounds);
worstIdentity=0; worstG=0; worstM=0; worstMC=0;
worstDerivative=0; worstSkew=0; worstLinearity=0; worstStatic=0;
worstMex=0; mexChecks=0; hasMex=(exist('armS_stationary_base_mex','file')==3);
for p=1:size(poses,2)
    q=poses(:,p); X=[q;dq];
    original=armS_standard_entry(0,X,L,r,mi,gWorld,K,D,tau,mu,bounds);
    [M0,C0,~,dM0]=armS_standard_core(q,dq,params);
    fd=zeros(6,6,6); step=1e-7;
    for h=1:6
        delta=zeros(6,1); delta(h)=step;
        Mp=armS_standard_core(q+delta,dq,params);
        Mm=armS_standard_core(q-delta,dq,params);
        fd(:,:,h)=(Mp-Mm)/(2*step);
    end
    derivative=norm(dM0(:)-fd(:))/max(norm(fd(:)),eps);
    Mdot=zeros(6);
    for h=1:6, Mdot=Mdot+dM0(:,:,h)*dq(h); end
    Z=Mdot-2*C0;
    skew=norm(Z+Z.','fro')/max(norm(Mdot,'fro')+2*norm(C0,'fro'),eps);
    worstDerivative=max(worstDerivative,derivative);
    worstSkew=max(worstSkew,skew);
    for k=1:numel(Rset)
        R=Rset{k}; gArm=R.'*gWorld; params.g=gArm;
        rhs=armS_stationary_base_entry(0,X,L,r,mi,R,gWorld,K,D,tau,mu,bounds);
        if hasMex
            mexRhs=armS_stationary_base_mex(0,X,L,r,mi,R,gWorld,K,D,tau,mu,bounds);
            worstMex=max(worstMex,max(abs(rhs-mexRhs)./(1+abs(rhs))));
            mexChecks=mexChecks+1;
        end
        [M,C,G,dM]=armS_standard_core(q,dq,params);
        [Mq,Gq]=worldQuadrature(q,L,r,mi,R,gWorld);
        mErr=norm(M-Mq,'fro')/max(norm(Mq,'fro'),eps);
        gErr=norm(G-Gq)/max(norm(Gq),1);
        mcErr=max([norm(M-M0,'fro'),norm(C-C0,'fro'),norm(dM(:)-dM0(:))]);
        worstM=max(worstM,mErr); worstG=max(worstG,gErr);
        worstMC=max(worstMC,mcErr);
        if k==1, worstIdentity=max(worstIdentity,max(abs(rhs-original))); end
        % G is linear in world gravity, independent of the stiffness law.
        gA=[1.3;-2.1;.7]; gB=[-.5;.4;-3.2];
        params.g=R.'*gA; [~,~,GA]=armS_standard_core(q,dq,params);
        params.g=R.'*gB; [~,~,GB]=armS_standard_core(q,dq,params);
        params.g=R.'*(gA+gB); [~,~,Gsum]=armS_standard_core(q,dq,params);
        worstLinearity=max(worstLinearity,norm(Gsum-GA-GB));
        % Static compensation must use +G with the existing RHS sign.
        Keff=K;
        for i=1:6
            Keff(i,i)=K(i,i)+.5*bounds(3)*(2 ...
                +tanh(mu*(q(i)-bounds(2)))-tanh(mu*(q(i)-bounds(1))));
        end
        staticRhs=armS_stationary_base_entry(0,[q;zeros(6,1)], ...
            L,r,mi,R,gWorld,K,D,G+Keff*q,mu,bounds);
        worstStatic=max(worstStatic,norm(staticRhs));
        fprintf('Pose %d orientation %d: M %.3e, G %.3e, M/C/dM %.3e\n', ...
            p,k,mErr,gErr,mcErr);
        assert(mErr<1e-5 && gErr<1e-8 && mcErr<1e-12, ...
            'Orientation or quadrature check failed.');
    end
end
if hasMex
    R=Rset{end}; changedGravity=[1.1;-2.2;-8.7];
    X=[poses(:,end);dq];
    matlabRhs=armS_stationary_base_entry(0,X,L,r,mi,R,changedGravity,K,D,tau,mu,bounds);
    mexRhs=armS_stationary_base_mex(0,X,L,r,mi,R,changedGravity,K,D,tau,mu,bounds);
    worstMex=max(worstMex,max(abs(matlabRhs-mexRhs)./(1+abs(matlabRhs))));
    mexChecks=mexChecks+1;
    assert(worstMex<1e-8,'MATLAB and MEX disagree.');
end
assert(worstIdentity<1e-12 && worstDerivative<1e-5 && worstSkew<1e-10 ...
    && worstLinearity<1e-10 && worstStatic<1e-8, ...
    'Identity, derivative, gravity linearity, or static-compensation check failed.');
assert(norm((Rset{2}.'*gWorld)-gWorld)<1e-12,'Yaw changed gravity.');
badR=eye(3); badR(1,1)=-1;
try
    armS_stationary_base_entry(0,zeros(12,1),L,r,mi,badR,gWorld,K,D,tau,mu,bounds);
    error('Reflection was accepted.');
catch err
    assert(~strcmp(err.message,'Reflection was accepted.'),'Reflection was accepted.');
end
badR=eye(3); badR(1,1)=1.01;
try
    armS_stationary_base_entry(0,zeros(12,1),L,r,mi,badR,gWorld,K,D,tau,mu,bounds);
    error('Nonorthonormal matrix was accepted.');
catch err
    assert(~strcmp(err.message,'Nonorthonormal matrix was accepted.'), ...
        'Nonorthonormal matrix was accepted.');
end
fprintf(['PASS: identity %.3e; world quadrature M %.3e, G %.3e; ' ...
    'orientation M/C/dM %.3e; dM %.3e; skew %.3e; ' ...
    'gravity linearity %.3e; static RHS %.3e; MEX %d checks %.3e\n'], ...
    worstIdentity,worstM,worstG,worstMC,worstDerivative,worstSkew, ...
    worstLinearity,worstStatic,mexChecks,worstMex);
end

function [M,G]=worldQuadrature(q,L,r,mi,Rbase,gWorld)
% Independent 16-node Gauss assembly from world point Jacobians.
k=(1:15).'; b=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1));
xi=(diag(D)+1)/2; w=(V(1,:).^2).';
R=eye(3); Pq=zeros(3,6); Rq=zeros(3,3,6);
M=zeros(6); G=zeros(6,1);
for n=1:3
    old=1:2*(n-1); cur=2*n-1:2*n; l=[0,q(cur).'];
    for s=1:16
        [~,~,p]=HTM_nume(l,xi(s),L,r);
        pj=LocalJacob_nume(l,xi(s),L,r);
        J=zeros(3,6);
        for a=old, J(:,a)=Pq(:,a)+Rq(:,:,a)*p; end
        J(:,cur)=R*pj;
        Jworld=Rbase*J;
        M=M+mi(n)*w(s)*(Jworld.'*Jworld);
        G=G+mi(n)*w(s)*Jworld.'*gWorld;
    end
    [~,Rt,p]=HTM_nume(l,1,L,r);
    [pj,rj]=LocalJacob_nume(l,1,L,r);
    for a=old
        Pq(:,a)=Pq(:,a)+Rq(:,:,a)*p;
        Rq(:,:,a)=Rq(:,:,a)*Rt;
    end
    Pq(:,cur)=R*pj;
    for a=1:2, Rq(:,:,cur(a))=R*rj(:,3*(a-1)+(1:3)); end
    R=R*Rt;
end
end
