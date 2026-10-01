function [theta_est, phi_est, fitResidual, info] = f20260219_2_task2config_withLenExt(p_target, L, xi, phi0, thetaPrior, straightRadiusTol)
% Estimate constant-curvature direction and bend angle from sensor XYZ.
% p_target is tip-sensor position relative to the base sensor, expressed in
% the continuum module's base frame. The module length L(2) is fixed.
%
% Inputs:
%   p_target         3x1 sensor position [x;y;z] (m)
%   L                [base-sensor offset; module length; tip-sensor offset] (m)
%   xi               backbone fraction at the tip sensor, 0 < xi <= 1
%   phi0             previous/initial bend estimate (rad); candidate in fit
%   thetaPrior       direction to retain when bending is unobservable (rad)
%   straightRadiusTol  lateral noise threshold (m); tune to the NDI setup
%
% Outputs:
%   theta_est        bending direction (rad)
%   phi_est          bend angle in [0,pi) (rad)
%   fitResidual      Euclidean sensor-position fit error (m)
%   info             thetaObservable, exitflag, iterations, method

    arguments (Input)
        p_target (3,1) double {mustBeFinite}
        L (3,1) double {mustBeFinite,mustBeNonnegative} = ndiSupportFnc.sensorGeometry()
        xi (1,1) double {mustBeFinite,mustBePositive} = 1
        phi0 (1,1) double {mustBeFinite} = 0.01
        thetaPrior (1,1) double {mustBeFinite} = 0
        straightRadiusTol (1,1) double {mustBeFinite,mustBeNonnegative} = 1e-4
    end

    if L(2) <= 0
        error('task2config:InvalidModuleLength', 'L(2) must be positive.');
    end
    if xi > 1
        error('task2config:InvalidXi', 'xi must be at most 1.');
    end

    rho = hypot(p_target(1), p_target(2));
    zMeasured = p_target(3);
    if rho <= straightRadiusTol
        % At a straight pose XYZ does not determine the bending direction.
        theta_est = atan2(sin(thetaPrior), cos(thetaPrior));
        phi_est = 0;
        zStraight = L(1) + xi*L(2) + L(3);
        fitResidual = hypot(rho, zMeasured-zStraight);
        info = struct('thetaObservable',false,'exitflag',1, ...
            'iterations',0,'method','straight');
        return;
    end

    theta_est = atan2(p_target(2), p_target(1));
    objective = @(phi) phiErrorSquared(phi, rho, zMeasured, L, xi);
    upper = pi - 1e-6;
    options = optimset('Display','off','TolX',1e-9);
    [phiOpt,~,exitflag,output] = fminbnd(objective, 0, upper, options);

    % Explicitly check the boundaries and the caller's previous estimate.
    candidates = [0, phiOpt, upper, min(max(phi0,0),upper)];
    errors = arrayfun(objective, candidates);
    [bestError, bestIndex] = min(errors);
    phi_est = candidates(bestIndex);
    fitResidual = sqrt(max(bestError,0));
    info = struct('thetaObservable',true,'exitflag',exitflag, ...
        'iterations',output.iterations,'method','bounded_scalar');
end

function err2 = phiErrorSquared(phi, rhoMeasured, zMeasured, L, xi)
    [rhoModel,zModel] = sensorRadialAxial(phi,L,xi);
    err2 = (rhoModel-rhoMeasured)^2 + (zModel-zMeasured)^2;
end

function [rho,z] = sensorRadialAxial(phi,L,xi)
    if phi == 0
        rho = 0;
        z = L(1) + xi*L(2) + L(3);
        return;
    end

    angle = xi*phi;
    % 1-cos(angle) = 2*sin(angle/2)^2 avoids cancellation near zero.
    rho = L(3)*sin(angle) + 2*L(2)*sin(angle/2)^2/phi;
    z = L(1) + L(3)*cos(angle) + L(2)*sin(angle)/phi;
end
