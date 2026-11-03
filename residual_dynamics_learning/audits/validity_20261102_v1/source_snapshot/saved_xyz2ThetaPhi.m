function [theta,phi,residual,valid] = xyz2ThetaPhi(rel_p,phi_prev)
persistent lastPhi lastTheta lastOutputTheta lastOutputPhi lastResidual
if isempty(lastPhi)
    lastPhi=0.0001; lastTheta=0;
    lastOutputTheta=NaN; lastOutputPhi=NaN; lastResidual=NaN;
end
% Keep the external memory for compatibility; never use its NaN.
if isfinite(phi_prev), lastPhi=phi_prev; end
[theta,phi,residual,valid]=ndiSupportFnc.configuration(rel_p,lastPhi,lastTheta);
if valid
    lastPhi=phi; lastTheta=theta;
    lastOutputTheta=theta; lastOutputPhi=phi; lastResidual=residual;
else
    % Scope outputs hold, but configurationValid remains false.
    theta=lastOutputTheta; phi=lastOutputPhi; residual=lastResidual;
end
end
