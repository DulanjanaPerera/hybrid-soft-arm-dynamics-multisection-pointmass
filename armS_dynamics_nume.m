function dX = armS_dynamics_nume(t,X, params)
% The function to compute the EoM matrices. The model of N-section
% continuum arm is computed. The joint-space varibles (length changes) are
% considered for the generalized coordinates.
% 
% Inputs:
%   t : time
%   X : length and length velocity [4Nx1]
%       [l2, l3, ..., dl2, dl3, ...]
%   params : parameters for matrices (structure)


N = params.N;
L = params.L;
r = params.r;
cog_xi = params.cog_xi;
mi = params.mi;
g = params.g;
K = params.K;
D = params.D;
tau = params.tau;

l = reshape(X(1:2*N,:)',[2,N])';
dl = reshape(X(2*N+1:end,:)',[2,N])';


Rglob = eye(3);
Pglob = zeros(3,1);

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

% M, C, G matrices,
M = zeros(2, 2);
C = zeros(2, 2);
G = zeros(2, 1);


dM = cell(N,1);   % dM{i} is (2i x 2i x 2i)
for i = 1:N
    ni = 2*i;               % number of DoFs up to section i
    dM{i} = zeros(ni, ni, ni); % [[matrix size], # DoF]
end

for n=1:N 
    tic
    fprintf("Section %d\n", n);
    % get the current length variables
    length = [0, l(n,1), l(n,2)];

    % compute the tip and CoG frame of n-th section
    [Ttip, Rtip, Ptip] = HTM_nume(length, 1, L, r);
    [Tcog, Rcog, Pcog] = HTM_nume(length, cog_xi(n), L, r);

    % compute the jacobians of the above frames
    [PJtip, RJtip, PJJtip, RJJtip] = LocalJacob_nume(length, 1, L, r);
    [PJcog, RJcog, PJJcog, RJJcog] = LocalJacob_nume(length, cog_xi(n), L, r);
     
    % Compute the global R and P for the n-th section
    Rglob = Rglob * Rtip;
    Pglob = Pglob + Rglob * Ptip;
    
    % compute the common block multiplications. This is the n-th section
    % Hessians pf Omega and Velocity. Here only bottom-right block's second
    % term is computed
    % R' * R_q_q  -- bottom-right block's second term in Hessian_Omega
    % R' * P_q_q  -- bottom-right block's second term in Hessian_vel
    temp_RRqq_mat = zeros([6,6]);
    temp_RPqq_mat = zeros([6,2]);
    temp_RRqq_mat_cog = zeros([6,6]);
    temp_RPqq_mat_cog =zeros([6,2]);

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

        M = PJcog.' * PJcog;

        % compute (M),h 
        for h=1:2*n % Here h={l11, l12, l21, l22, ..., ln1, ln2, ...}
            dM{n}(:,:,h) = Mi_h(n, h, Pcog, PJcog, PJJcog, J_veltip, J_Omegatip, H_veltip, H_Omegatip);
        end

        G = ( (mi(n) * g' * Rglob * (Rglob * PJcog) ) + ( K(1:2*n,1:2*n) * length(2:3).' ).' ).';

    else
        % need to know how many block are there in the Jacobian (for
        % example, [3x6N] = [3x(3*b)]. For N=2|b=2, N=3|b=4 ==> (2(n-1))
        blocks = size(J_Omegatip,2)/3;
        
        % create a temporary matrices to compute block multiplication.
        temp_JoR_mat = zeros([3,3*blocks]);
        temp_JoP_mat = zeros([3, 1*blocks]);

        temp_JoR_mat_cog = zeros([3,3*blocks]);
        temp_JoP_mat_cog = zeros([3, 1*blocks]);

        temp_RHoR_mat = zeros([3*blocks, 3*blocks]); % Rtip.' * H_Omegatip * Rtip
        temp_RqJoR_mat = zeros([6, 3*blocks]); % RJtip.' * temp_Omega_mat
        temp_RJoRq_mat = zeros([6, 3*blocks]); % RJtip.' * temp_Omega_mat
        temp_RHvHoP_mat = zeros([3*blocks, blocks]); % Rtip.' * (H_veltip + H_Omegatip*Ptip)
        temp_RqJvJoP_mat = zeros([6, blocks]); % Rq.' * (J_veltip + J_Omegatip * Ptip)
        temp_RJoPq_mat = zeros([6, blocks] ); % Rtip.' J_Omega * PJtip

        temp_RHoR_mat_cog = zeros([3*blocks, 3*blocks]); % Rtip.' * H_Omegatip * Rtip
        temp_RqJoR_mat_cog = zeros([6, 3*blocks]); % RJtip.' * temp_Omega_mat
        temp_RJoRq_mat_cog = zeros([6, 3*blocks]); % RJtip.' * temp_Omega_mat
        temp_RHvHoP_mat_cog = zeros([3*blocks, blocks]); % Rtip.' * (H_veltip + H_Omegatip*Ptip)
        temp_RqJvJoP_mat_cog = zeros([6, blocks]); % Rq.' * (J_veltip + J_Omegatip * Ptip)
        temp_RJoPq_mat_cog = zeros([6, blocks]); % Rtip.' J_Omega * PJtip

        % MASS matrix
        temp_sigma_11 = zeros([blocks,blocks]);
        temp_sigma_12 = zeros([blocks, 2]);
        temp_sigma_22 = PJcog.' * PJcog;

        % (M),h matrix
        temp_neta_11 = zeros([blocks, blocks]);
        

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
            temp_JoP_mat_cog(:,b) = J_Omegatip(:,3*(b-1)+1: 3*(b-1)+3) * Pcog; % this is re-used in MASS matrix

       

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

                
                % sigma_11 of MASS matrix without JoP.'JoP
                % J_veltip' * (J_veltip + 2*J_Omegatip * Pcog)
                % [2(n)x3] * ( [3x2(n)] + [3x6(n)]*[3x1])
                % [2nx2n]
                temp_sigma_11(r, b) = J_veltip(:,r).' * (J_veltip(:,b) + 2 * temp_JoP_mat_cog(:,b));

               

            end
            
            % this loop is maily for bottom -left blocks where 6xkN
            % is required. Here k={2, 6}
            % Also it is used for generate EoM matrices
            for r=1:2
                % --------------- Tip -------------------------------------
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
                
                % ----------------- MASS matrix ---------------------------
                % sigma_12 of MASS matrix 
                % here this loop goes column wise and blocks are used to
                % row-wise
                % (J_veltip + J_Omega*Pcog)' * PJcog
                % ([3x2(n)] + [3x6(n)]*[3x1])' * [3x2]
                % [2nx2]
                temp_sigma_12(b,r) = (J_veltip(:,b) + temp_JoP_mat_cog(:,b)).' * PJcog(:,r);
            end
            
            

        end

        % Update the M, C, amd G matrices
        M = mi(n)*[M + temp_sigma_11 + (temp_JoP_mat_cog.' * temp_JoP_mat_cog), temp_sigma_12;
            temp_sigma_12.', PJcog.' * PJcog];
        
        % compute (M),h
        for h=1:2*n % Here h={l11, l12, l21, l22, ..., ln1, ln2, ...}
            dM{n}(:,:,h) = Mi_h(n, h, Pcog, PJcog, PJJcog, J_veltip, J_Omegatip, H_veltip, H_Omegatip);
            
        end

        % constructing G matrix
        Gp = mi(n) * g' * Rglob * ([J_veltip + temp_JoP_mat_cog, Rglob * PJcog]);
        Ge = [zeros(1,2*(n-1)), (K(2*(n-1)+1:2*(n-1)+2,2*(n-1)+1:2*(n-1)+2) * length(2:3).' ).'];
        G = [G; zeros(2,1)] + Gp.' + Ge.';
        
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
    toc
end % end of section loop

% constucting C matrices
for n=1:N
    if n==1
        l_d = reshape(dl(1:n,:)',[2*n,1]);
        C = C + christoffelSymbol(n, dM{n}, l_d);
    else
        l_d = reshape(dl(1:n,:)',[2*n,1]);
        Ci = christoffelSymbol(n, dM{n}, l_d);
        Ci(1:2*(n-1), 1:2*(n-1)) = Ci(1:2*(n-1), 1:2*(n-1)) + C;
        C = Ci;
    end
end
flat_dl = reshape(dl(1:N,:)',[2*N,1]);
ddl = M \ (tau - (C + D)*flat_dl - G);

dX = [flat_dl; ddl];

end