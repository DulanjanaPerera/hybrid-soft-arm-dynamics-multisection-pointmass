function [M, C, G] = armS_core_N3_mex(t, l, dl, L, r, cog_xi, mi, g, K)
%#codegen
% N fixed to 3
[N, ~] = size(cog_xi);
% N  = length(cog_xi);

% ---- paste your recursion here ----
% Replace:
%   params.N etc.
% with local variables above.
%
% IMPORTANT for Coder:
%   - Avoid cell arrays (replace dM{n} with numeric array if needed)
%   - Avoid growing arrays via concatenation
%   - Preallocate fixed-size matrices:


% You currently use:
%   dM = cell(N,1)
% Replace with fixed 4D array:
%   dM_all(:,:,h,n) if you need per-section
% or fixed 3D per-section buffers and accumulate.
%
% For N=3, easiest:
%   dM1 = zeros(2,2,2);
%   dM2 = zeros(4,4,4);
%   dM3 = zeros(6,6,6);
%
% Then call christoffelSymbol with those.

% ... your loop n=1:3, fill M, G, and fill dM1/dM2/dM3

Rglob = eye(3);
Pglob = zeros(3,1);

% The tip velocity Jacobians and Hessian
J_Omegatip = zeros([3,6*N]);
J_veltip = zeros([3,2*N]);
H_Omegatip = zeros([6*N,6*N]);
H_veltip = zeros([6*N,2*N]);

% M, C, G matrices,
M = zeros(2*N, 2*N);
C = zeros(2*N, 2*N);
G = zeros(2*N, 1);

r_offset = r;

dM1 = zeros(2,2,2);
dM2 = zeros(4,4,4);
dM3 = zeros(6,6,6);

for n=1:N 
    % get the current length variables
    length = [0, l(n,1), l(n,2)];

    % compute the tip and CoG frame of n-th section
    [~, Rtip, Ptip] = HTM_nume_mex(length, 1, L, r_offset);
    [~, Rcog, Pcog] = HTM_nume_mex(length, cog_xi(n), L, r_offset);

    % compute the jacobians of the above frames
    [PJtip, RJtip, PJJtip, RJJtip] = LocalJacob_nume_mex(length, 1, L, r_offset);
    [PJcog, ~, PJJcog, ~] = LocalJacob_nume_mex(length, cog_xi(n), L, r_offset);
     

    % compute the common blocks at Hessian.
    % bottom-right block's 1st term
    % RJtip.' * RJtip -- Hessian_Omega
    % RJtip.' * PJtip -- Hessian_vel
    % however, mutiplication is different
    %
    % RJtip.' * RJtip = [
    %                    (RJtip_1)'*RJtip_1, (RJtip_2)'*RJtip_1;
    %                    (RJtip_1)'*RJtip_2, (RJtip_2)'*RJtip_2
    %                   ]
    

    % compute the common block multiplications. This is the n-th section
    % Hessians pf Omega and Velocity. Here only bottom-right block's second
    % term is computed
    % R' * R_q_q  -- bottom-right block's second term in Hessian_Omega
    % R' * P_q_q  -- bottom-right block's second term in Hessian_vel
    temp_RqRq_mat = zeros(6,6);                                                 % CHECKED 2025/01/07
    temp_RqPq_mat = zeros([6,2]);                                               % CHECKED 2025/01/07
    temp_RRqq_mat = zeros([6,6]);                                               % CHECKED 2025/01/06 
    temp_RPqq_mat = zeros([6,2]);                                               % CHECKED 2025/01/07

    for c=1:2 % always 2 blocks, (column)
        for r=1:2 % row
        
            % 
            % RJJtip = [dR/dl2l2 dR/dl3l2; 
            %           dR/dl2dl3 dR/dl3l3]
            % PJJtip = [dP/dl2l2 dP/dl2l3;
            %           dP/dl3l2 dP/dl3l3]
            
            % bottom-right block's 1st term in Hessian_vel tip
            % The Hessian has followig format
            %
            % (H_vel)_{jk} = (J_vel)_k,qj
            % H_vel = [
            %          J1,q1  J2,q1  | J3,q1 ...
            %          J2,q2  J2,q2  | J3,q2 ...
            %         ]
            temp_RqPq_mat(3*(r-1)+1: 3*(r-1)+3, c) = RJtip(:, 3*(r-1)+1: 3*(r-1)+3).' * PJtip(:, c);            % CHECKED 2025/01/08
            
            % bottom-right block's 1st term in Hessian_Omega 
            temp_RqRq_mat(3*(r-1)+1: 3*(r-1)+3, 3*(c-1)+1: 3*(c-1)+3) = ...
                RJtip(:, 3*(r-1)+1: 3*(r-1)+3).' * RJtip(:, 3*(c-1)+1: 3*(c-1)+3);                              % CHECKED 2025/01/08

            % bottom-right block's second term in Hessian_Omega tip
            % extract [3x3] block and multiply is with [3x3] Rtip'. assign
            % it back to the [3x3]
            % Rtip.' * RJJtip
            % [3x3] * [6x6] = [6x6]
            temp_RRqq_mat(3*(r-1)+1: 3*(r-1)+3, 3*(c-1)+1: 3*(c-1)+3) = ...
                Rtip.' * RJJtip(3*(c-1)+1: 3*(c-1)+3, 3*(r-1)+1: 3*(r-1)+3);                                    % CHECKED 2025/01/08
 
            
            % bottom-right block's second term in Hessian_vel tip -  No need to
            % compute the first term because it is simply RJtip' * PJtip
            % Rtip.' * PJJtip
            % [3x3] * [6x2] = [6x2]
            temp_RPqq_mat(3*(r-1)+1:3*(r-1)+3, c) = ...
                Rtip.' * PJJtip(3*(c-1)+1:3*(c-1)+3, r);                                                        % CHECKED 2025/01/08
    
        end
                    
    end


    % Accumulate the Jacobians for the n-th section
    % For the first iteration, standard velocity jacobian is applied. But
    % thereafter, the recursive method is applied.
    if n==1

        J_Omegatip(:,1:6*(n)) = Rtip.' * RJtip;                                % CHECKED 2025/01/06                             
        J_veltip(:,1:2*(n)) = Rtip.' *  PJtip;                                 % CHECKED 2025/01/06
        H_Omegatip(1:6*(n),1:6*(n)) = temp_RqRq_mat + temp_RRqq_mat;                 % CHECKED 2025/01/06
        H_veltip(1:6*(n),1:2*(n)) = temp_RqPq_mat + temp_RPqq_mat;                   % CHECKED 2025/01/07

        M(1:2*n, 1:2*n) = mi(n) * (PJcog.' * PJcog);
        M(1:2*n, 1:2*n) = 0.5*(M(1:2*n, 1:2*n) + M(1:2*n, 1:2*n).');
        % compute (M),h 
        for h=1:2*n % Here h={l11, l12, l21, l22, ..., ln1, ln2, ...}
            dM1(:,:,h) = Mi_h(n, h, Pcog, PJcog, PJJcog, J_veltip(:,1:2*(n)), J_Omegatip(:,1:6*(n)), H_veltip(1:6*(n),1:2*(n)), H_Omegatip(1:6*(n),1:6*(n)));     
        end
        

        % G = ( (mi(n) * g' * Rglob * (Rglob * PJcog) ) + ( K(1:2*n,1:2*n) * length(2:3).' ).' ).';
        G(1:2*n,1) = ( ( mi(n) * g' * (Rglob * PJcog) ) ).';                      % CHECKED 2025/01/15 (removed the stiffness pot
        
        % Compute the global R and P for the n-th section 
        Pglob = Pglob + Rglob * Ptip;
        Rglob = Rglob * Rtip;


    else
        % need to know how many block are there in the Jacobian (for
        % example, [3x6N] = [3x(3*b)]. For N=2|b=2, N=3|b=4 ==> (2(n-1))
        blocks = 2*(n-1);
        
        % create a temporary matrices to compute block multiplication.
        temp_JoR_mat = zeros([3,3*blocks]);                                                             % CHECKED 2025/01/06                        
        temp_JoP_mat = zeros([3, blocks]);                                                              % CHECKED 2025/01/06 

        temp_JoR_mat_cog = zeros([3,3*blocks]);                                                         % CHECKED 2025/01/06 
        temp_JoP_mat_cog = zeros([3, blocks]);                                                          % CHECKED 2025/01/06 

        temp_RHoR_mat = zeros([3*blocks, 3*blocks]); % Rtip.' * H_Omegatip * Rtip                       % CHECKED 2025/01/06 
        temp_RqJoR_mat = zeros([6, 3*blocks]); % RJtip.' * temp_JoR_mat                                 % CHECKED 2025/01/06 
        temp_RJoRq_mat = zeros([6, 3*blocks]); % RJtip.' * temp_Omega_mat * RJtip                       % CHECKED 2025/01/07 
        temp_RHvHoP_mat = zeros([3*blocks, blocks]); % Rtip.' * (H_veltip + H_Omegatip*Ptip)            % CHECKED 2025/01/07
        temp_RqJvJoP_mat = zeros([6, blocks]); % Rq.' * (J_veltip + J_Omegatip * Ptip)                  % CHECKED 2025/01/07
        temp_RJoPq_mat = zeros([6, blocks] ); % Rtip.' J_Omega * PJtip                                  % CHECKED 2025/01/07 

        % MASS matrix (sigmas)
        temp_sigma_11 = zeros([blocks,blocks]);
        temp_sigma_12 = zeros([blocks, 2]);

        

        for c=1:blocks % column loop
            % left block of the Jacobian_Omega tip and CoG
            % extract 3x3 block and multiply it with 3x3. Then assign it back to
            % the 3x3 block
            temp_JoR_mat(:,3*(c-1)+1: 3*(c-1)+3) = ...
                J_Omegatip(:,3*(c-1)+1: 3*(c-1)+3) * Rtip;                              % CHECKED 2025/01/06            
            temp_JoR_mat_cog(:,3*(c-1)+1: 3*(c-1)+3) = ...
                J_Omegatip(:,3*(c-1)+1: 3*(c-1)+3) * Rcog;                              % CHECKED 2025/01/06 
            
            % 2nd term of left block of the Jacobian_velocity tip and CoG
            % extract 3x3 block and multiply is with 3x1. Then assign it back to
            % the 3x1 block
            temp_JoP_mat(:,c) = J_Omegatip(:,3*(c-1)+1: 3*(c-1)+3) * Ptip;              % CHECKED 2025/01/06 
            temp_JoP_mat_cog(:,c) = J_Omegatip(:,3*(c-1)+1: 3*(c-1)+3) * Pcog;          % CHECKED 2025/01/06 

       

            % This loop is mainly for top-left blocks where NxN is required
            for r=1:blocks % row loop
                
                % --------------- Tip -------------------------------------
                % Top-left block of Hessian_Omega tip
                % Rtip.' * Ho * Rtip
                % [3x3] * [3x3] * [3x3]
                temp_RHoR_mat(3*(r-1)+1: 3*(r-1)+3, 3*(c-1)+1: 3*(c-1)+3) = ...
                    Rtip.' * H_Omegatip(3*(r-1)+1: 3*(r-1)+3, 3*(c-1)+1: 3*(c-1)+3) * Rtip;                                         % CHECKED 2025/01/06                    
                
                % Top-left block of Hessian_velocity tip
                % Rtip.' * (Hvtip + HOtip*Ptip)
                % [3x3] * ([3x1(n)] + [3x3]*[3x1]) -- this is for one
                % iteration
                temp_RHvHoP_mat(3*(r-1)+1: 3*(r-1)+3, c) = ...
                    Rtip.' * (H_veltip(3*(r-1)+1: 3*(r-1)+3, c) + (H_Omegatip(3*(r-1)+1: 3*(r-1)+3, 3*(c-1)+1: 3*(c-1)+3) * Ptip));   % CHECKED 2025/01/07
                
                % sigma_11 of MASS matrix without JoP.'JoP
                % J_veltip' * (J_veltip + 2*J_Omegatip * Pcog)
                % [2(n)x3] * ( [3x2(n)] + [3x6(n)]*[3x1])
                % [2nx2n]
                temp_sigma_11(r, c) = J_veltip(:,r).' * ( J_veltip(:,c) + 2 * temp_JoP_mat_cog(:,c) );                              % CHECKED 2025/01/07

               

            end % row loop end
            
            % this loop is mainly for bottom -left blocks where 6xkN
            % is required. Here k={2, 6}
            % Also it is used for generate EoM matrices
            for r=1:2
                % --------------- Tip -------------------------------------
                % bottom-left block's first term in Hessian_Omega tip
                % RJtip.' * J_Omegatip * Rtip
                % [6x3] * [3x6(n)] * [3x3]
                % [6x3] * [3x6(n)] = [6x6(n)]
                temp_RqJoR_mat(3*(r-1)+1: 3*(r-1)+3, 3*(c-1)+1: 3*(c-1)+3) = ...
                   RJtip(:,3*(r-1)+1: 3*(r-1)+3).' * temp_JoR_mat(:,3*(c-1)+1: 3*(c-1)+3);                              % CHECKED 2025/01/07
            
                % bottom-left block's second term in Hessian_Omega tip
                % Rtip.' * J_Omegatip * RJtip
                % [3x3] * [3x6(n)] * [3x6]
                temp_RJoRq_mat(3*(r-1)+1: 3*(r-1)+3, 3*(c-1)+1: 3*(c-1)+3) = ...
                    Rtip.' * J_Omegatip(:,3*(c-1)+1: 3*(c-1)+3) * RJtip(:,3*(r-1)+1: 3*(r-1)+3);                        % CHECKED 2025/01/07       

                % bottom-left block's first term in Hessian_vel tip
                % RJtip.' * (J_veltip + J_Omegatip * Ptip)
                % [6x3] * ([3x2(n)] + [3x6(n)] * [3x1])
                temp_RqJvJoP_mat(3*(r-1)+1: 3*(r-1)+3, c) = ...
                    RJtip(:,3*(r-1)+1: 3*(r-1)+3).' * (J_veltip(:,c) + (J_Omegatip(:,3*(c-1)+1: 3*(c-1)+3) * Ptip));    % CHECKED 2025/01/07

                % bottom-left block's second term in Hessian_vel tip
                % Rtip' * J_Omegatip * PJtip
                % [3x3] * [3x6(n)] * [3x2]
                temp_RJoPq_mat(3*(r-1)+1: 3*(r-1)+3, c) = ...
                    Rtip.' * J_Omegatip(:,3*(c-1)+1: 3*(c-1)+3) * PJtip(:,r);                                            % CHECKED 2025/01/07 
                
                % ----------------- MASS matrix ---------------------------
                % sigma_12 of MASS matrix 
                % here this loop goes column wise and blocks are used to
                % row-wise
                % (J_veltip + J_Omega*Pcog)' * PJcog
                % ([3x2(n)] + [3x6(n)]*[3x1])' * [3x2]
                % [2nx2]
                temp_sigma_12(c,r) = ( J_veltip(:,c) + temp_JoP_mat_cog(:,c) ).' * PJcog(:,r);                          % CHECKED 2025/01/07
            end
            
            

        end
        

        % Update the M, C, amd G matrices
        M(1:2*n, 1:2*n) = [M(1:2*(n-1), 1:2*(n-1)) + mi(n) * (temp_sigma_11 + (temp_JoP_mat_cog.' * temp_JoP_mat_cog)), mi(n) * temp_sigma_12;
             mi(n) * temp_sigma_12.', mi(n) * (PJcog.' * PJcog)];                                                       % CHECKED 2025/01/07
        M(1:2*n, 1:2*n) = 0.5*(M(1:2*n, 1:2*n) + M(1:2*n, 1:2*n).');
        % compute (M),h
        if n==2
            for h=1:2*n % Here h={l11, l12, l21, l22, ..., ln1, ln2, ...}
                % dM2(:,:,h) = Mi_h(n, h, Pcog, PJcog, PJJcog, J_veltip, J_Omegatip, H_veltip, H_Omegatip);
                dM2(:,:,h) = Mi_h(n, h, Pcog, PJcog, PJJcog, J_veltip(:,1:2*(n-1)), J_Omegatip(:,1:6*(n-1)), H_veltip(1:6*(n-1),1:2*(n-1)), H_Omegatip(1:6*(n-1),1:6*(n-1)));
            end
        elseif n == 3
            for h=1:2*n % Here h={l11, l12, l21, l22, ..., ln1, ln2, ...}
                dM3(:,:,h) = Mi_h(n, h, Pcog, PJcog, PJJcog, J_veltip(:,1:2*(n-1)), J_Omegatip(:,1:6*(n-1)), H_veltip(1:6*(n-1),1:2*(n-1)), H_Omegatip(1:6*(n-1),1:6*(n-1))); 
            end
        end
        

        % constructing G matrix
        % Gp = mi(n) * g' * Rglob * ([J_veltip + temp_JoP_mat_cog, Rglob * PJcog]);
        Gp = mi(n) * (g' * Rglob * [ (J_veltip(:,1:2*(n-1)) + temp_JoP_mat_cog), PJcog]);
        % Ge = [zeros(1,2*(n-1)), (K(2*(n-1)+1:2*(n-1)+2,2*(n-1)+1:2*(n-1)+2) * length(2:3).' ).'];
        G(1:2*n, 1) = [G(1:2*(n-1), 1); zeros(2,1)] + Gp.';
        
        % Compute the global R and P for the n-th section 
        Pglob = Pglob + Rglob * Ptip;
        Rglob = Rglob * Rtip;

        % IMPORTANT: update the CoG first and then tip. Because, if tip is
        % updated then the J_vel and J_Omega are n-th section, and not
        % (n-1)-th section
        % J_velcog = Rcog.' * [J_veltip + temp_JoP_mat_cog, PJcog];
        J_veltip(:,1:2*(n)) = Rtip.' * [J_veltip(:,1:2*(n-1)) + temp_JoP_mat, PJtip];

        % J_Omegacog = Rcog.' * [temp_JoR_mat_cog, RJcog];
        J_Omegatip(:,1:6*(n)) = Rtip.' * [temp_JoR_mat, RJtip];


        H_Omegatip(1:6*(n),1:6*(n)) = [temp_RHoR_mat, zeros(size(temp_RHoR_mat,1),6);
            temp_RqJoR_mat + temp_RJoRq_mat, temp_RqRq_mat + temp_RRqq_mat];

        % n-th section Hessian_vel tip
        H_veltip(1:6*(n),1:2*(n)) = [temp_RHvHoP_mat, zeros(size(temp_RHvHoP_mat,1), 2);
            temp_RqJvJoP_mat + temp_RJoPq_mat, temp_RqPq_mat + temp_RPqq_mat];

    
    end
    
end % end of section loop

% constucting C matrices
for n=1:N
    if n==1
        l_d = reshape(dl(1:n,:)',[2*n,1]);
        C(1:2*n, 1:2*n) = C(1:2*n, 1:2*n) + christoffelSymbol(n, dM1, l_d);
    elseif n == 2
        l_d = reshape(dl(1:n,:)',[2*n,1]);
        Ci = christoffelSymbol(n, dM2, l_d);
        Ci(1:2*(n-1), 1:2*(n-1)) = Ci(1:2*(n-1), 1:2*(n-1)) + C(1:2*(n-1), 1:2*(n-1));
        C(1:2*n, 1:2*n) = Ci;
    elseif n == 3
        l_d = reshape(dl(1:n,:)',[2*n,1]);
        Ci = christoffelSymbol(n, dM3, l_d);
        Ci(1:2*(n-1), 1:2*(n-1)) = Ci(1:2*(n-1), 1:2*(n-1)) + C(1:2*(n-1), 1:2*(n-1));
        C(1:2*n, 1:2*n) = Ci;
    end
end


% ... build C from dM1/dM2/dM3 exactly as you already do ...

end
