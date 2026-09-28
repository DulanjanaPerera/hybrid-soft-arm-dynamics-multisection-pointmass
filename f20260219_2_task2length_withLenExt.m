function [l, fitResidual, info] = f20260219_2_task2length_withLenExt(p_target, L, xi, r, phi0, thetaPrior, straightRadiusTol)
% Estimate three PMA length changes from tip-sensor XYZ (m).
% Use l(2:3) as the two single-module dynamics coordinates.
%
% p_target is relative to the base sensor in the module's base frame.
% L = [base-sensor offset; fixed module length; tip-sensor offset] (m).
% phi0 and thetaPrior may be taken from the previous sample. Set
% straightRadiusTol to the measured lateral NDI noise scale (m).
% Optional outputs report XYZ fit error (m) and solver diagnostics.

    arguments (Input)
        p_target (3,1) double {mustBeFinite}
        L (3,1) double {mustBeFinite,mustBeNonnegative} = [0.05;0.1778;0.05]
        xi (1,1) double {mustBeFinite,mustBePositive} = 1
        r (1,1) double {mustBeFinite,mustBePositive} = 0.013
        phi0 (1,1) double {mustBeFinite} = 0.01
        thetaPrior (1,1) double {mustBeFinite} = 0
        straightRadiusTol (1,1) double {mustBeFinite,mustBeNonnegative} = 1e-4
    end

    [theta,phi,fitResidual,info] = f20260219_2_task2config_withLenExt( ...
        p_target,L,xi,phi0,thetaPrior,straightRadiusTol);

    l = r*phi*[-cos(theta);
        cos(theta)/2 - sqrt(3)*sin(theta)/2;
        cos(theta)/2 + sqrt(3)*sin(theta)/2];
end
