function [theta_est, phi_est] = f20260219_2_task2config_withLenExt(p_target, L, xi, phi0)
% Estimate theta and phi from a 3D vector [x; y] using hybrid method:
% - theta is computed directly from
% - phi is estimated by minimizing position error
%
% Inputs:
%   p_target : 3x1 vector [x; y; z] - measured position
%   L        : length of the section and extensions [3x1] (m)
%			   [from base to base-sensor linear offset; length of continuum arm;
%			    from tip of the arm to tip-sensor linear offset ]
%   xi       : selection factor [constant] {0, 1}
%   phi0     : initial guess for phi (optional)
%
% Outputs:
%   theta_est : estimated theta (rad)
%   phi_est   : estimated phi   (rad)

    arguments (Input)
        p_target (3,1) double
        L (3,1) double {mustBePositive} = [0.05; 0.1778; 0.05]
        xi (1,1) double {mustBePositive} = 1.0
        phi0 (1,1) double = 0.01
    end

    arguments (Output)
        theta_est (1,1) double
        phi_est (1,1) double
    end

    theta_est = 0;
    phi_est = 0;
    % If position is (near) zero, best we can do is theta=0, phi -> 0
    if norm(p_target) < 1e-9
        theta_est = 0;
        phi_est   = 1e-6; % tiny positive to avoid division by zero downstream
        return;
    end

    % Step 1: Estimate theta using atan2
    theta_est = atan2(p_target(2), p_target(1));

    % Step 2: Estimate phi using least squares
    % options = optimoptions('Levenberg-Marquardt', ...
    %     'Display', 'off', ...
    %     'FunctionTolerance', 1e-12, ...
    %     'StepTolerance', 1e-12, ...
    %     'MaxIterations', 200, ...
    %     'MaxFunctionEvaluations', 400);

    options = optimoptions('fmincon','Display','off','Algorithm','sqp');

    fun = @(phi) modelErrorPhi(phi, p_target, theta_est, L, xi);

    % Solve for phi only
    phi_est = fmincon(fun, phi0, [], [], [], [], 1e-6, pi-1e-6, [], options);

    % Wrap results (optional)
    theta_est = wrapToPi(theta_est);

    % Step 3: Wrap phi estimate to ensure it is within the range [-pi, pi]
    phi_est = wrapToPi(phi_est);

    % Make sure that value is (0, pi)
    phi_est = min(max(phi_est, 1e-6), pi - 1e-6);
end

function err = modelErrorPhi(phi, p_target, theta, L, xi)


    % Avoid division by zero
    if abs(phi) < 1e-6
        phi = sign(phi) * 1e-6;
    end

    % Forward kinematics with fixed theta
    p = [cos(theta) * sin(xi * phi) * L(3) - cos(theta) * cos(xi * phi) * L(2) / phi + cos(theta) * L(2) / phi sin(theta) * sin(xi * phi) * L(3) - sin(theta) * cos(xi * phi) * L(2) / phi + sin(theta) * L(2) / phi cos(xi * phi) * L(3) + sin(xi * phi) * L(2) / phi + L(1)];
    p = p(:);

    % Residual
    err = norm(p - p_target);
end
