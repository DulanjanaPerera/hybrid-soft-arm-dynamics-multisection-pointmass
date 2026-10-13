function [state,active,atBound]=pacAdaptiveUpdate(state,Y,s,dt,allow)
%PACADAPTIVEUPDATE Box-projected, rate-bounded gradient. No hidden I/O.
% Outputs are estimates for the NEXT controller sample.
%#codegen
cfg=pacAdaptiveConfig();active=false;
atBound=any(state.theta<=cfg.lower | state.theta>=cfg.upper);
if ~allow || dt<=0 || dt>cfg.maxUpdateInterval_s || ...
 ~all(isfinite(Y(:))) || ~all(isfinite(s)) || norm(s)<=cfg.deadband_mPerS,return;end
% zeta=theta/scale; Yscaled=Y*diag(scale), zeta_dot=-Gamma*Yscaled'*s.
rate=-(cfg.scale.^2).*cfg.gamma.*(Y.'*s);
rate=min(max(rate,-cfg.maxRate),cfg.maxRate);
rate((state.theta<=cfg.lower & rate<0) | ...
     (state.theta>=cfg.upper & rate>0))=0;
state.theta=min(max(state.theta+dt*rate,cfg.lower),cfg.upper);
active=any(rate~=0);atBound=any(state.theta<=cfg.lower | state.theta>=cfg.upper);
end
