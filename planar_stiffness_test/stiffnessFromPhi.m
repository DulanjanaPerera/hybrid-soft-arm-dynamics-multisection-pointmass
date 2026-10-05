function [k,inRange] = stiffnessFromPhi(phi,fit)
% Empirical P1 stiffness in N/m, phi in radians. Outside fitted range: NaN.
% Pass fit explicitly to avoid reloading it when repeatedly evaluating.
if nargin<2
 here=fileparts(mfilename('fullpath'));
 s=load(fullfile(here,'results','stiffness_phi_fit.mat'),'fit');fit=s.fit;
end
inRange=isfinite(phi)&phi>=fit.PhiRange_rad(1)&phi<=fit.PhiRange_rad(2);
k=nan(size(phi));k(inRange)=ppval(fit.PP,phi(inRange));
end
