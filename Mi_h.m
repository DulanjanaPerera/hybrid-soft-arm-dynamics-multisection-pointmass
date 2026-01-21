function Mih = Mi_h(n, h, p, p_q, p_qq, J_vel, J_omega, H_vel, H_omega)
%#codegen
% This function derives the partial derivatives of M (MASS) matrix for the
% given joint-space varible (h). The h here is an integer in {1,2,3,...,2n}
% For example, if n=2, then h can be {1,2,3,4}. However, h={1,2}
% corresponds to the previous section while, h={3,4} corresponds to the
% current section.
% The construction of (M_n),h uses common term for both eta and gamma to simplify the code.
% 
% The derivations are referenced from following paper.
%   * Godage, Isuru S., Robert J. Webster, and Ian D. Walker. "Center-of-gravity-based 
%     approach for modeling dynamics of multisection continuum arms." IEEE transactions 
%     on robotics 35, no. 5 (2019): 1097-1108.
% 
% Inputs:
%   n       : Current section [constant]
%   h       : current DoF [constant] (h can be one of the {1,2,3,4,...,2*n}
%   p       : CoG local position [3x1]
%   p_q     : dp/dq [3x2]
%   p_qq    : d2p/dq2 [6x2]
%   J_vel   : Velocity Jacobian [3x2(n-1)]
%   J_Omega : Angular velocity Jacobian [3x6(n-1)]
%   H_vel   : linear velocity hessian [6(n)x2(n)]
%   H_Omega : angular velocity hessian [6(n)x6(n)]
% 
% Output:
%   Mih     : The partial derivitives of M w.r.t. h [2nx2n]


d = 2*(n-1);

% For h < 2(n-1)
eta_11 = zeros(d,d);
eta_12 = zeros(d,2);
eta_22 = zeros(2,2);

% for 2(n-1)<h % using a common term for the both eta and gamma
% gamma_11 = zeros(d,d);
% gamma_12 = zeros(d,2);
% gamma_22 = zeros(2,2);

% beta_v1 = 1;
% beta_v2 = 1;
% beta_v3 = 1;
beta_v1 = 1.3155;
beta_v2 = 1.5053;
beta_v3 = 1.8011;

for r=1:d
    for c=1:d
        if(h <= 2*(n-1))
            % This is the eta calculation. when h <= 2*(n-1) that means h
            % is previous section's DoF

            % ci = 2 - mod(c,2); % this is to obtain the index 1 and 2 for odd and even numbers of column

            eta_11(r,c) = 2 * H_vel(3*(h-1)+1:3*(h-1)+3,r).' * (J_vel(:,c) + J_omega(:,3*(c-1)+1:3*(c-1)+3)*p) ...
                + 2 * (J_vel(:,r).' + beta_v1 * (J_omega(:,3*(r-1)+1:3*(r-1)+3)* p).') * (H_omega(3*(h-1)+1:3*(h-1)+3, 3*(c-1)+1:3*(c-1)+3) * p);                                  % CHECKED 2025/01/07
            
            if c<=2 % only two column exist
                eta_12(r,c) = (H_vel(3*(h-1)+1:3*(h-1)+3, r) ...
                    + beta_v2 * H_omega(3*(h-1)+1:3*(h-1)+3,3*(r-1)+1:3*(r-1)+3) * p).' ...
                    * p_q(:,c);                                                                             % CHECKED 2025/01/07
            end

        else
            % This is for gamma calculation. If h>2(n-1), then h is current
            % section's DoF.

            hi = 2 - mod(h,2); % this is to obtain the index 1 and 2 for odd and even numbers of h

            eta_11(r,c) = 2 * J_vel(:,r).' * (J_omega(:,3*(c-1)+1:3*(c-1)+3) * p_q(:,hi)) ...
                + 2 * beta_v1 *(J_omega(:,3*(r-1)+1:3*(r-1)+3) * p_q(:,hi)).' * (J_omega(:,3*(c-1)+1:3*(c-1)+3) * p); % CHECKED 2025/01/07

            if c<=2
                eta_12(r,c) = (J_vel(:,r) + beta_v2 * J_omega(:,3*(r-1)+1:3*(r-1)+3) * p).' ...
                    * p_qq(3*(c-1)+1:3*(c-1)+3, hi)...
                    + beta_v2 * (J_omega(:,3*(r-1)+1:3*(r-1)+3) * p_q(:,hi)).' * p_q(:,c);                             % CHECKED 2025/01/07 
                    % changed the row column in p_qq  % CHECKED 2025/01/15
                if r <=2
                    eta_22(r,c) = 2 * beta_v3 * p_qq(3*(r-1)+1:3*(r-1)+3, hi).' * p_q(:,c);                                % CHECKED 2025/01/07
                    % changed the row column in p_qq % CHECKED 2025/01/15
                end

            end
        end
    end

end

if d == 0
    eta_22 = zeros(2,2);
    for r=1:2
        for c=1:2
            eta_22(r,c) = 2 * beta_v3 * p_qq(3*(r-1)+1:3*(r-1)+3, h).' * p_q(:,c);                                        % CHECKED 2025/01/07
        end
    end
    eta_22 = 0.5*(eta_22 + eta_22.');
    % if section 1, then only gamma_22 is present
    Mih = eta_22;
else
    eta_11 = 0.5*(eta_11 + eta_11.');
    eta_22 = 0.5*(eta_22 + eta_22.');
    % from section 2 onwards, both eta and gamma exist based on the h
    Mih = [eta_11, eta_12;
          eta_12.', eta_22];
end
% miz = eye(size(Mih));
% Mih = miz;
% rcond(Mih)
% Mih(~isfinite(Mih)) = 0;

end