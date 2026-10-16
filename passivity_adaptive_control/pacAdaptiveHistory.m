function z=pacAdaptiveHistory(z,dq,dt,cfg)
%PACADAPTIVEHISTORY Bounded normalized Bouc-Wen n=1 feature states.
% Internal alpha/beta/gamma fixed; force weights adapt, not hidden-state law.
%#codegen
if nargin<4,cfg=pacAdaptiveConfig();end
if dt<=0 || dt>cfg.maxUpdateInterval_s || ~all(isfinite(dq)),return;end
w=[-dq(1)-dq(2);dq(1);dq(2)]/cfg.historyLength_m;
% Stable substeps for large measured velocity; avoid explicit-Euler explosions.
n=min(200,max(1,ceil(dt*max(abs(w))/.05)));h=dt/n;
for j=1:n
 rate=w-cfg.historyBeta*abs(w).*z-cfg.historyGamma*w.*abs(z);
 z=min(max(z+h*rate,-cfg.historyBound),cfg.historyBound);
end
end
