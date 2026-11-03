function [hidden,force]=forceResidualCell(params,x,hidden,cfg)
%FORCERESIDUALCELL GRU reset-after cell, two bounded force outputs in newtons.
H=cfg.hiddenUnits;u=params.input*x+params.bias;w=params.recurrent*hidden;
reset=1./(1+exp(-(u(1:H,:)+w(1:H,:))));
update=1./(1+exp(-(u(H+1:2*H,:)+w(H+1:2*H,:))));
candidate=tanh(u(2*H+1:3*H,:)+reset.*w(2*H+1:3*H,:));
hidden=update.*hidden+(1-update).*candidate;
head=params.output*hidden+params.skip*x+params.outputBias;
force=cfg.forceBound_N*tanh(cfg.forceScale_N*head/cfg.forceBound_N);
end
