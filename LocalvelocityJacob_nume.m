function [J_Omega, J_vel] = LocalvelocityJacob_nume(l, xi, L, r)
% Compute the angular and linear velocity jacobians. The angular velocity
% jacobian is skewsymmetric matrix of the angular velocity vector.
% 
% The derivations are referenced from following paper.
%   * Godage, Isuru S., Robert J. Webster, and Ian D. Walker. "Center-of-gravity-based 
%     approach for modeling dynamics of multisection continuum arms." IEEE transactions 
%     on robotics 35, no. 5 (2019): 1097-1108.
% 
% Inputs:
%   l   : length change of the PMA [3x1] [constant] (m)
%   xi  : selection factor of the backbone. xi=0 is the base and xi=1 is
%         the tip. [constant] {0,1}
%   L   : Length of the continuum arm [constant] (m)
%   r   : radial offset of the PMA [constant] (m)
% 
% Output:
%   J_Omega : The angular velocity Jacobian [3x6]
%   J_vel   : The linear velocity Jacobian [3x2]
%   

arguments (Input)
    l (1,3) double
    xi (1,1) double {mustBeBetween(xi, 0, 1, "closed")} = 1
    L (1,1) double {mustBeGreaterThanOrEqual(L, 0)} = 0.278
    r (1,1) double {mustBePositive} = 0.013
end

arguments (Output)
    J_Omega (3,6) double
    J_vel (3,2) double
end


[~, R, ~] = HTM_sym(l, xi, L, r);
[PosJ, RotJ] = LocalJacob_sym(l, xi, L, r);

J_Omega = R' * RotJ;
J_vel = R' * PosJ;

end