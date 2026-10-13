function state=pacAdaptiveInitialState()
%#codegen
cfg=pacAdaptiveConfig();state.theta=cfg.initial;
state.z=zeros(3,1);state.lastSample=-1;
end
