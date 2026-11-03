function [state,force]=sgCell(params,x,state,cfg)
%SGCELL Causal stacked GRU. Inputs are already normalized; state is explicit.
original=x;
for layer=1:numel(cfg.widths)
    H=cfg.widths(layer);tag=sprintf('l%d',layer);old=state{layer};
    u=params.([tag 'Input'])*x+params.([tag 'Bias']);
    w=params.([tag 'Recurrent'])*old;
    reset=.5*(1+tanh(.5*(u(1:H,:)+w(1:H,:))));
    update=.5*(1+tanh(.5*(u(H+1:2*H,:)+w(H+1:2*H,:))));
    candidate=tanh(u(2*H+1:3*H,:)+reset.*w(2*H+1:3*H,:));
    state{layer}=update.*old+(1-update).*candidate;x=state{layer};
end
head=params.output*x+params.skip*original+params.outputBias;
force=cfg.forceBound_N*tanh(cfg.forceScale_N*head/cfg.forceBound_N);
end
