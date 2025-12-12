% The the local jacobians here has following dimensions
% 
% PosJ*     [3x2]
% RotJ*     [3x6]
% J_vel*    [3x2]
% J_Omega*  [3x6]


clear
N = 2;
L = 0.278;
r = 0.013;
cog_xi = [0.5, 0.5];
syms 'l' [N 3] 

J_Omegatip = zeros([3,6]);
J_veltip = zeros([3,2]);
for n=1:N 

    % get the current length variables
    length = [0, l(n,2), l(n,3)];

    % compute the tip and CoG frame of n-th section
    [Ttip, Rtip, Ptip] = HTM_sym(length, 1, L, r);
    [Tcog, Rcog, Pcog] = HTM_sym(length, cog_xi(n), L, r);

    % compute the jacobians of the above frames
    [PJtip, RJtip] = LocalJacob_sym(length, 1, L, r);
    [PJcog, RJcog] = LocalJacob_sym(length, cog_xi(n), L, r);
    
    % % only the tip velocity jacobians are enough
    % [J_Omegatip, J_veltip] = LocalvelocityJacob_sym(length, 1, L, r);

    % Accumulate the Jacobians for the n-th section
    % For the first iteration, standard velocity jacobian is applied. But
    % thereafter, the recursive method is applied.
    if n==1
        J_Omegatip = Rtip' * RJtip;
        J_veltip = Rtip' *  J_Omegatip;
    else
        % need to know how many block are there in the Jacobian (for
        % example, [3x6N]
        blocks = size(J_Omegatip,2)/3;
        
        % create a temporary matrices to compute block multiplication.
        temp_Omega_mat = zeros(3,3*blocks);
        temp_vel_mat = zeros(3, 1*blocks);

        for b=1:blocks
            % extract 3x3 block and multiply it with 3x3. Then assign it back to
            % the 3x3 block
            temp_Omega_mat(:,3*(b-1)+1: 3*(b-1)+3) = J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Rtip;

            % extract 3x3 block and multiply is with 3x1. Then assign it back to
            % the 3x1 block
            temp_vel_mat(:,2*(b-1)+1: 2*(b-1)+2) = J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Ptip;
        end

        J_Omegatip = Rtip' * [temp_Omega_mat, RJtip];
        J_veltip = Rtip' * [J_veltip + temp_vel_mat, PJtip];
    end

end