function [tau,s,Y,fi,fc,fd,fs,ff,fh,fb,Kh,Dh,valid]= ...
    pacAdaptiveForce(measured,desired,gravity,lambda,Ks,state,inputValid)
%PACADAPTIVEFORCE Nominal M/C/G, structured adaptive K/D/history/bias.
%#codegen
tau=zeros(2,1);s=zeros(2,1);Y=zeros(2,10);
fi=tau;fc=tau;fd=tau;fs=tau;ff=tau;fh=tau;fb=tau;
Kh=zeros(2);Dh=zeros(2);valid=false;
if ~inputValid || ~all(isfinite(measured(:))) || ...
 ~all(isfinite(desired(:))) || ~all(isfinite(gravity(:))) || ...
 ~all(isfinite(lambda(:))) || ~all(isfinite(Ks(:))) || ...
 any(lambda(:)<=0) || any(Ks(:)<=0) || ...
 ~all(isfinite(state.theta)) || ~all(isfinite(state.z))
    return;
end
[p,~,~]=pacNominalParameters(gravity);
q=measured(1:2);dq=measured(3:4);qd=desired(1:2);
v=desired(3:4)-lambda(:).*(q-qd);
a=desired(5:6)-lambda(:).*(dq-desired(3:4));s=dq-v;
[M,C,G]=armS_single_core(q,dq,p);
if ~all(isfinite(M(:))) || ~all(isfinite(C(:))) || ...
 ~all(isfinite(G)) || M(1,1)<=0 || det(M)<=0,return;end
th=state.theta;
Kh=[th(1)+th(3),th(3);th(3),th(2)+th(3)];Dh=diag(th(4:5));
% Project each physical muscle history force into reduced length coordinates.
T=[-1 -1;1 0;0 1];
Y(:,1)=[q(1);0];Y(:,2)=[0;q(2)];Y(:,3)=[q(1)+q(2);q(1)+q(2)];
Y(:,4)=[v(1);0];Y(:,5)=[0;v(2)];
Y(:,6:8)=T.'*diag(state.z);Y(:,9:10)=eye(2);
bound=zeros(2,1);
for i=1:2
 bound(i)=.5*p.lKbounds(3)*(2+tanh(p.mu*(q(i)-p.lKbounds(2)))- ...
     tanh(p.mu*(q(i)-p.lKbounds(1))))*q(i);
end
fi=M*a;fc=C*v;fd=Dh*v;fs=G+Kh*q+bound;ff=-Ks(:).*s;
fh=Y(:,6:8)*th(6:8);fb=th(9:10);
tau=fi+fc+fd+fs+ff+fh+fb;valid=all(isfinite(tau));
if ~valid,tau(:)=0;end
end
