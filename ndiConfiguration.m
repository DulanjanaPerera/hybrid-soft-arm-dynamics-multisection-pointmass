function [theta,phi,residual,valid] = ndiConfiguration(p,phiPrior,thetaPrior)
% Geometry boundary: never pass a missing reading into finite-only IK.
% Priors must be the last valid configuration, not a missing sample.
theta=NaN; phi=NaN; residual=NaN;
valid=all(isfinite(p));
if ~valid, return; end
if ~isfinite(phiPrior), phiPrior=0.0001; end
if ~isfinite(thetaPrior), thetaPrior=0; end
[theta,phi,residual] = f20260219_2_task2config_withLenExt( ...
    p,ndiSensorGeometry(),1,phiPrior,thetaPrior);
end
