% function dX  = DynamicModel(p, X)
%#codegen
% This file computes length velocity and accelration from pressures using the dynamic model.
%
% Inputs:
%   p   : pressures [3x1] (Pa)
%   X   : length derivatives [4x1] (m)
%         [l2; l3; dl2; dl3]
%
% Output:
%   dX  : length derivatives [4x1] (m)
%         [dl2; dl3; ddl2; ddl3]

function [dX, lengths_m, theta_rad, phi_rad, tau_N, ...
          pressureUsed_bar, thetaObservable] = forwardModel(pressureCommand_bar, X, stiff)
%#codegen
% pressureCommand_bar = [P1; P2; P3], nonnegative bar
% X = [l2; l3; dl2; dl3], with l1 = -l2-l3

L = 0.17411;
r = 0.013; % PMA radial offset from backbone (m)
mass = 0.1;
K = 650 * [2 1; 1 2];
% K = stiff * [2 1; 1 2];
D = 40 * eye(2);
gArm = [0; 0; -9.81];
mu = 2000;
lengthBounds = [-0.035; 0.035; 1e6];

pressureUsed_bar = min(max(pressureCommand_bar, 0), 3);
% pressureUsed_bar = max(min(max(pressureCommand_bar, 0), 3) - 0.8, 0);

pressureDiameter_m = 0.013; % Effective pressure-area diameter, independent of r
pressureRadius_m = pressureDiameter_m/2;
area = pi * pressureRadius_m^2;

% Positive P1, P2, P3 give the established -X, -60°, +60°
% nominal bending directions, respectively.
tau_N = area * 1e5 * [-1 1 0; -1 0 1] * pressureUsed_bar;

dX = armS_single_entry(0, X, L, r, mass, gArm, ...
                       K, D, tau_N, mu, lengthBounds);

l2 = X(1);
l3 = X(2);
lengths_m = [-l2-l3; l2; l3];

u = (l2+l3)/r;
v = (l3-l2)/(sqrt(3)*r);
phi_rad = hypot(u, v);
theta_rad = atan2(v, u);
thetaObservable = phi_rad > 0.5*pi/180;


end