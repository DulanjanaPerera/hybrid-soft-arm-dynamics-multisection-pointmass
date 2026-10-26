function [moduleTipXYZ_m,tipXYZValid] = moduleTipXYZBlock( ...
    sensorRelativeXYZ_m,R_NDI_from_arm,theta,phi, ...
    orientationValid,baseOrientationValid)
%#codegen
% Input XYZ: existing Sensor/position output, metres, AFTER rotateX90.
% Output XYZ: physical module endpoint relative to module base, expressed
% in the CURRENT calibrated module-base frame. Not world-frame position.
% Assumes axial sensor mounting offsets: base .055 m, distal .050 m.
moduleTipXYZ_m=nan(3,1);
tipXYZValid=false;
if ~orientationValid || ~baseOrientationValid || ...
        ~all(isfinite(sensorRelativeXYZ_m(:))) || ...
        ~all(isfinite(R_NDI_from_arm(:))) || ...
        ~isfinite(theta) || ~isfinite(phi)
    return;
end
if norm(R_NDI_from_arm.'*R_NDI_from_arm-eye(3),'fro')>1e-3 || ...
        abs(det(R_NDI_from_arm)-1)>1e-3
    return;
end
% Undo the legacy fixed plotting-frame rotation, then use measured base pose.
R_fixed=[1 0 0;0 0 1;0 -1 0];
sensorXYZ_arm=R_NDI_from_arm.'*R_fixed.'*sensorRelativeXYZ_m(:);
% Measured tip tangent reconstructed from orientation-derived theta/phi.
% This does not reconstruct XYZ from a constant-curvature arc.
tipTangent=[sin(phi)*cos(theta);sin(phi)*sin(theta);cos(phi)];
moduleTipXYZ_m=sensorXYZ_arm-[0;0;0.055]-0.050*tipTangent;
tipXYZValid=true;
end
