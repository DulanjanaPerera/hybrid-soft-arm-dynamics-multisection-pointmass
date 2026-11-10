function [next,force,Jp]=residualGruCell(network,x,previous)
%RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate evaluation.
% x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history fixed.
%#codegen
z=(x-network.mean)./network.std;
Jp=zeros(2,3);
u=network.input*z+network.bias;w=network.recurrent*previous;
r=.5*(1+tanh(.5*(u(1:128)+w(1:128))));
g=.5*(1+tanh(.5*(u(129:256)+w(129:256))));
c=tanh(u(257:384)+r.*w(257:384));
next=g.*previous+(1-g).*c;
head=network.output*next+network.skip*z+network.outputBias;
a=tanh(network.forceScale*head/network.forceBound);force=network.forceBound*a;
if nargout<3,return;end
W=network.input(:,5:7)./network.std(5:7)';
dr=(r.*(1-r)).*W(1:128,:);dg=(g.*(1-g)).*W(129:256,:);
dc=(1-c.^2).*(W(257:384,:)+w(257:384).*dr);
dh=(1-g).*dc+(previous-c).*dg;
Jp=network.forceScale*(1-a.^2).*(network.output*dh+network.skip(:,5:7)./network.std(5:7)');
end
