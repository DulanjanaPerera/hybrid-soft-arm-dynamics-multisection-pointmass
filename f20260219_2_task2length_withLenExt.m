function l = f20260219_2_task2length_withLenExt(p_target, L, xi, r, phi0)
% Estimate lengths from a 3D vector [x; y; z] for extended system:
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
%   l        : estimated length change [3x1] (m)


    arguments (Input)
        p_target (3,1) double
        L (3,1) double {mustBePositive} = [0.05; 0.1778; 0.05]
        xi (1,1) double {mustBePositive} = 1.0
        r (1,1) double {mustBePositive} = 0.013
        phi0 (1,1) double = 0.01
    end

    arguments (Output)
        l (3,1) double
    end

    [theta, phi] = f20260219_2_task2config_withLenExt(p_target, L, xi, phi0);

    l = [-r * cos(theta) * phi (0.5e0 * r * cos(theta) - sqrt(0.3e1) * r * sin(theta) / 0.2e1) * phi (0.5e0 * r * cos(theta) + sqrt(0.3e1) * r * sin(theta) / 0.2e1) * phi];
    l = l(:);
end
