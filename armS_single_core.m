function [M,C,G,dM] = armS_single_core(q,dq,params)
% One distributed section, translational soft-material kinetic energy only.
% q,dq = [l12;l13] and their rates; dm = mi*dXi.
assert(params.N==1 && numel(q)==2 && numel(dq)==2 ...
    && numel(params.mi)==1);
q=q(:); dq=dq(:);
l=[0,q.'];
[~,muq]=integratedPosition_nume(l,params.L,params.r);
[E,Eq]=integratedJacobianProduct_compact(l,params.L,params.r);
M=params.mi*E;
dM=params.mi*Eq;
G=params.mi*muq.'*params.g;
C=christoffelSymbol(1,dM,dq);
end
