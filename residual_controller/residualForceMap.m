function [total,correction,raw,next,J]=residualForceMap(p,measured,dt,h,network,cfg,gravityShift)
%RESIDUALFORCEMAP Activation re-expression keeps F_cfg+r_eff=B0*p+r0+deltaG.
%#codegen
J=zeros(2,3);
if nargout>=5,[next,raw,Jr]=residualGruCell(network,[measured(1:4);p;dt],h);
else,[next,raw]=residualGruCell(network,[measured(1:4);p;dt],h);Jr=zeros(2,3);end
B=network.forcePerBar*[-1 1 0;-1 0 1];ep=max(p-cfg.activation_bar,0);
nominal=cfg.forcePerBar_N*[ep(2)-ep(1);ep(3)-ep(1)];
total=B*p+raw+gravityShift;correction=total-nominal;J=B+Jr;
end
