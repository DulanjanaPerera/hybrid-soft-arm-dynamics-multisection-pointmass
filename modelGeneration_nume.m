% The the local jacobians here has following dimensions
% 
% PJ*       [3x2]
% RJ*       [3x6]
% PJJ*      [6x2]
% RJJ*      [6x6]
% J_vel*    [3x2(n)]
% J_Omega*  [3x6(n)]
% H_Omega*  [6(n)x6(n)]
% H_vel*    [6(n)x2(n)]
% 
% Here * indicates the "tip" or "cog" in the naming convention. The n
% indicates the number of sections in the arm.



clear
N = 2;
L = 0.278;
r = 0.013;
cog_xi = 0.5*ones(1,N);
l = [
     0.0, 0.010, 0.010;
     0.0, 0.015, 0.015
    ];

% The tip velocity Jacobians and Hessian
J_Omegatip = zeros([3,6]);
J_veltip = zeros([3,2]);
H_Omegatip = zeros([6,6]);
H_veltip = zeros([6,2]);

% The COG velocity Jacobians
J_Omegacog = zeros([3,6]);
J_velcog = zeros([3,2]);
H_Omegacog = zeros([6,6]);
H_velcog = zeros([6,2]);

for n=1:N 
    fprintf("Section %d\n", n);
    % get the current length variables
    length = [0, l(n,2), l(n,3)];

    % compute the tip and CoG frame of n-th section
    [Ttip, Rtip, Ptip] = HTM_sym(length, 1, L, r);
    [Tcog, Rcog, Pcog] = HTM_sym(length, cog_xi(n), L, r);

    % compute the jacobians of the above frames
    [PJtip, RJtip, PJJtip, RJJtip] = LocalJacob_sym(length, 1, L, r);
    [PJcog, RJcog, PJJcog, RJJcog] = LocalJacob_sym(length, cog_xi(n), L, r);


    % compute the common block multiplications. This is the n-th section
    % Hessians pf Omega and Velocity. Here only bottom-right block's second
    % term is computed
    % R' * R_q_q  -- bottom-right block's second term in Hessian_Omega
    % R' * P_q_q  -- bottom-right block's second term in Hessian_vel
    syms temp_RRqq_mat [6,6] matrix
    syms temp_RPqq_mat [6,2] matrix
    syms temp_RRqq_mat_cog [6,6] matrix
    syms temp_RPqq_mat_cog [6,2] matrix

    for b=1:2 % always 2 blocks,
        
        % 
        % RJJtip = [dR/dl2l2 dR/dl2l3; 
        %           dR/dl3dl2 dR/dl3l3]
        % PJJtip = [dP/dl2l2 dP/dl2l3;
        %           dP/dl3l2 dP/dl3l3]

        % bottom-right block's second term in Hessian_Omega tip
        % extract [3x3] block and multiply is with [3x3] Rtip'. assign
        % it back to the [3x3]
        % Rtip.' * RJJtip
        % [3x3] * [6x6]
        temp_RRqq_mat(1:3, 3*(b-1)+1:3*(b-1)+3) = ...
            Rtip.' * RJJtip(1:3, 3*(b-1)+1:3*(b-1)+3);
        temp_RRqq_mat(4:6, 3*(b-1)+1:3*(b-1)+3) = ...
            Rtip.' * RJJtip(4:6, 3*(b-1)+1:3*(b-1)+3);
        
        % bottom-right block's second term in Hessian_Omega CoG
        % Rcog.' * RJJcog
        % [3x3] * [6x6]
        temp_RRqq_mat_cog(1:3, 3*(b-1)+1:3*(b-1)+3) = ...
            Rcog.' * RJJcog(1:3, 3*(b-1)+1:3*(b-1)+3);
        temp_RRqq_mat_cog(4:6, 3*(b-1)+1:3*(b-1)+3) = ...
            Rcog.' * RJJcog(4:6, 3*(b-1)+1:3*(b-1)+3);
        
        % bottom-right block's second term in Hessian_vel tip -  No need to
        % compute the first term because it is simply RJtip' * PJtip
        % Rtip.' * PJJtip
        % [3x3] * [6x2]
        temp_RPqq_mat(3*(b-1)+1:3*(b-1)+3, 1:2) = ...
            Rtip.' * PJJtip(3*(b-1)+1:3*(b-1)+3, 1:2);

        % bottom-right block's second term in Hessian_vel CoG - No need to
        % compute the first term because it is simply RJcog' * PJcog
        % Rcog.' * PJJcog
        % [3x3] * [6x2]
        temp_RPqq_mat_cog(3*(b-1)+1:3*(b-1)+3, 1:2) = ...
            Rcog.' * PJJcog(3*(b-1)+1:3*(b-1)+3, 1:2);
    end

    % Accumulate the Jacobians for the n-th section
    % For the first iteration, standard velocity jacobian is applied. But
    % thereafter, the recursive method is applied.
    if n==1

        J_Omegatip = Rtip.' * RJtip;
        J_veltip = Rtip.' *  PJtip;
        H_Omegatip = RJtip.' * RJtip + temp_RRqq_mat;
        H_veltip = RJtip.' * PJtip + temp_RPqq_mat;


        J_Omegacog = Rcog.' * RJcog;
        J_velcog = Rcog.' *  PJcog;
        H_Omegacog = RJcog.' * RJcog + temp_RRqq_mat_cog;
        H_velcog = RJcog.' * PJcog + temp_RPqq_mat_cog;

    else
        % need to know how many block are there in the Jacobian (for
        % example, [3x6N]
        blocks = size(J_Omegatip,2)/3;
        
        % create a temporary matrices to compute block multiplication.
        syms temp_JoR_mat [3,3*blocks] matrix;
        syms temp_JoP_mat [3, 1*blocks] matrix;

        syms temp_JoR_mat_cog [3,3*blocks] matrix;
        syms temp_JoP_mat_cog [3, 1*blocks] matrix;

        syms temp_RHoR_mat [3*blocks, 3*blocks] matrix % Rtip.' * H_Omegatip * Rtip
        syms temp_RqJoR_mat [6, 3*blocks] matrix % RJtip.' * temp_Omega_mat
        syms temp_RJoRq_mat [6, 3*blocks] matrix % RJtip.' * temp_Omega_mat
        syms temp_RHvHoP_mat [3*blocks, blocks] matrix % Rtip.' * (H_veltip + H_Omegatip*Ptip)
        syms temp_RqJvJoP_mat [6, blocks] matrix % Rq.' * (J_veltip + J_Omegatip * Ptip)
        syms temp_RJoPq_mat [6, blocks] matrix % Rtip.' J_Omega * PJtip

        syms temp_RHoR_mat_cog [3*blocks, 3*blocks] matrix % Rtip.' * H_Omegatip * Rtip
        syms temp_RqJoR_mat_cog [6, 3*blocks] matrix % RJtip.' * temp_Omega_mat
        syms temp_RJoRq_mat_cog [6, 3*blocks] matrix % RJtip.' * temp_Omega_mat
        syms temp_RHvHoP_mat_cog [3*blocks, blocks] matrix % Rtip.' * (H_veltip + H_Omegatip*Ptip)
        syms temp_RqJvJoP_mat_cog [6, blocks] matrix % Rq.' * (J_veltip + J_Omegatip * Ptip)
        syms temp_RJoPq_mat_cog [6, blocks] matrix % Rtip.' J_Omega * PJtip

        for b=1:blocks % column loop
            % left block of the Jacobian_Omega tip and CoG
            % extract 3x3 block and multiply it with 3x3. Then assign it back to
            % the 3x3 block
            temp_JoR_mat(:,3*(b-1)+1: 3*(b-1)+3) = ...
                J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Rtip;
            temp_JoR_mat_cog(:,3*(b-1)+1: 3*(b-1)+3) = ...
                J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Rcog;
            
            % 2nd term of left block of the Jacobian_velocity tip and CoG
            % extract 3x3 block and multiply is with 3x1. Then assign it back to
            % the 3x1 block
            temp_JoP_mat(:,b) = J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Ptip;
            temp_JoP_mat_cog(:,b) = J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Pcog;

            
            % This loop is mainly for top-left blocks where NxN is required
            for r=1:blocks % row loop
                
                % --------------- Tip -------------------------------------
                % Top-left block of Hessian_Omega tip
                % Rtip.' * Ho * Rtip
                % [3x3] * [3x3] * [3x3]
                temp_RHoR_mat(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) = ...
                    Rtip.' * H_Omegatip(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) * Rtip;
                
                % Top-left block of Hessian_velocity tip
                % Rtip.' * (Hvtip + HOtip*Ptip)
                % [3x3] * ([3x1(n)] + [3x3]*[3x1]) -- this is for one
                % iteration
                temp_RHvHoP_mat(3*(r-1)+1: 3*(r-1)+3, b) = ...
                    Rtip.' * (H_veltip(3*(r-1)+1: 3*(r-1)+3, b) + H_Omegatip(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) * Ptip);
            
                % --------------- CoG -------------------------------------
                % Top-left block of Hessian_Omega cog
                % Rcog.' * Ho * Rcog
                % [3x3] * [3x3] * [3x3]
                temp_RHoR_mat_cog(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) = ...
                    Rcog.' * H_Omegatip(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) * Rcog;
                
                % Top-left block of Hessian_velocity cog
                % Rcog.' * (Hvtip + HOtip*Pcog)
                % [3x3] * ([3x1(n)] + [3x3]*[3x1]) -- this is for one
                % iteration
                temp_RHvHoP_mat_cog(3*(r-1)+1: 3*(r-1)+3, b) = ...
                    Rcog.' * (H_veltip(3*(r-1)+1: 3*(r-1)+3, b) + H_Omegatip(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) * Pcog);

            end
            
            % this loop is maily for bottom -left blocks where 6xkN
            % is required. Here k={2, 6}
            for r=1:2
                % bottom-left block's first term in Hessian_Omega tip
                % RJtip.' * J_Omegatip * Rtip
                % [6x3] * [3x6(n)] * [3x3]
                temp_RqJoR_mat(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) = ...
                   RJtip(:,3*(r-1)+1: 3*(r-1)+3).' * temp_JoR_mat(:,3*(b-1)+1: 3*(b-1)+3);
            
                % bottom-left block's second term in Hessian_Omega tip
                % Rtip.' * J_Omegatip * RJtip
                % [3x3] * [3x6(n)] * [3x6]
                temp_RJoRq_mat(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) = ...
                    Rtip.' * J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * RJtip(:,3*(r-1)+1: 3*(r-1)+3);

                % bottom-left block's first term in Hessian_vel tip
                % RJtip.' * (J_veltip + J_Omegatip * Ptip)
                % [6x3] * ([3x2(n)] + [3x6(n)] * [3x1])
                temp_RqJvJoP_mat(3*(r-1)+1: 3*(r-1)+3, b) = ...
                    RJtip(:,3*(r-1)+1: 3*(r-1)+3).' * (J_veltip(:,b) + J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Ptip);

                % bottom-left block's second term in Hessian_vel tip
                % Rtip' * J_Omegatip * PJtip
                % [3x3] * [3x6(n)] * [3x2]
                temp_RJoPq_mat(3*(r-1)+1: 3*(r-1)+3, b) = ...
                    Rtip' * J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * PJtip(:,r);

                % -------------- CoG --------------------------------------
                % bottom-left block's first term in Hessian_Omega cog
                % RJcog.' * J_Omegatip * Rcog
                % [6x3] * [3x6(n)] * [3x3]
                temp_RqJoR_mat_cog(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) = ...
                   RJcog(:,3*(r-1)+1: 3*(r-1)+3).' * temp_JoR_mat_cog(:,3*(b-1)+1: 3*(b-1)+3);
            
                % bottom-left block's second term in Hessian_Omega cog
                % Rcog.' * J_Omegatip * RJcog
                % [3x3] * [3x6(n)] * [3x6]
                temp_RJoRq_mat_cog(3*(r-1)+1: 3*(r-1)+3, 3*(b-1)+1: 3*(b-1)+3) = ...
                    Rcog.' * J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * RJcog(:,3*(r-1)+1: 3*(r-1)+3);

                % bottom-left block's first term in Hessian_vel cog
                % RJcog.' * (J_veltip + J_Omegatip * Pcog)
                % [6x3] * ([3x2(n)] + [3x6(n)] * [3x1])
                temp_RqJvJoP_mat_cog(3*(r-1)+1: 3*(r-1)+3, b) = ...
                    RJcog(:,3*(r-1)+1: 3*(r-1)+3).' * (J_veltip(:,b) + J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Pcog);

                % bottom-left block's second term in Hessian_vel cog
                % Rcog' * J_Omegatip * PJcog
                % [3x3] * [3x6(n)] * [3x2]
                temp_RJoPq_mat_cog(3*(r-1)+1: 3*(r-1)+3, b) = ...
                    Rcog' * J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * PJcog(:,r);

            end
            
            

        end

        

        % IMPORTANT: update the CoG first and then tip. Because, if tip is
        % updated then the J_vel and J_Omega are n-th section, and not
        % (n-1)-th section
        J_velcog = Rcog.' * [J_veltip + temp_JoP_mat_cog, PJcog];
        J_veltip = Rtip.' * [J_veltip + temp_JoP_mat, PJtip];

        J_Omegacog = Rcog.' * [temp_JoR_mat_cog, RJcog];
        J_Omegatip = Rtip.' * [temp_JoR_mat, RJtip];

        % n-th section Hessian_Omega cog
        H_Omegacog = [temp_RHoR_mat_cog, zeros(size(temp_RHoR_mat_cog,1),6);
            temp_RqJoR_mat_cog + temp_RJoRq_mat_cog, RJcog.' * RJcog + temp_RRqq_mat_cog];
        % n-th section Hessian_vel cog
        H_velcog = [temp_RHvHoP_mat_cog, zeros(size(temp_RHvHoP_mat_cog,1), 2);
            temp_RqJvJoP_mat_cog + temp_RJoPq_mat_cog, RJcog.' * PJcog + temp_RPqq_mat_cog];

        % n-th section Hessian_Omega tip
        H_Omegatip = [temp_RHoR_mat, zeros(size(temp_RHoR_mat,1),6);
            temp_RqJoR_mat + temp_RJoRq_mat, RJtip.' * RJtip + temp_RRqq_mat];
        % n-th section Hessian_vel tip
        H_veltip = [temp_RHvHoP_mat, zeros(size(temp_RHvHoP_mat,1), 2);
            temp_RqJvJoP_mat + temp_RJoPq_mat, RJtip.' * PJtip + temp_RPqq_mat];


    end

end