function Mih = Mi_h(n, h, p, p_q, p_qq, J_vel, J_omega, H_vel, H_omega)
% This function derives the partial derivatives of M (MASS) matrix for the
% given joint-space varible (h). The construction of (M_n),h uses common 
% term for both eta and gamma to simplify the code.
% 
% The derivations are referenced from following paper.
%   * Godage, Isuru S., Robert J. Webster, and Ian D. Walker. "Center-of-gravity-based 
%     approach for modeling dynamics of multisection continuum arms." IEEE transactions 
%     on robotics 35, no. 5 (2019): 1097-1108.
% 
% Inputs:
%   n       : Current section [constant]
%   h       : current derivative {constant]
%   p       : CoG local position [3x1]
%   p_q     : dp/dq [3x2]
%   p_qq    : d2p/dq2 [6x2]
%   J_vel   : Velocity Jacobian [3x2(n-1)]
%   J_Omega : Angular velocity Jacobian [3x6(n-1)]
%   H_vel   : linear velocity hessian [6(n-1)x2(n-1)]
%   H_Omega : angular velocity hessian [6(n-1)x6(n-1)]
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

for r=1:d
    for c=1:d
        if(h <= 2*(n-1))
            
            ci = 2 - mod(c,2); % this is to obtain the index 1 and 2 for odd and even numbers of column

            eta_11(r,c) = 2 * H_vel(3*(h-1)+1:3*(h-1)+3,r).' ...
                * (J_vel(:,c) + J_omega(:,3*(c-1)+1:3*(c-1)+3)*p) ...
                + 2 * (J_vel(:,r).' + (J_omega(:,3*(r-1)+1:3*(r-1)+3)* p).') ...
                * (H_omega(3*(h-1)+1:3*(h-1)+3, 3*(c-1)+1:3*(c-1)+3) * p);
            
            if c<=2 % only two column exist
                eta_12(r,c) = (H_vel(3*(h-1)+1:3*(h-1)+3, r) ...
                    + H_omega(3*(h-1)+1:3*(h-1)+3,3*(r-1)+1:3*(r-1)+3) * p).' ...
                    * p_q(:,ci);
            end

        else
            
            h = 2 - mod(h,2); % this is to obtain the index 1 and 2 for odd and even numbers of h
            ci = 2 - mod(c,2); % this is to obtain the index 1 and 2 for odd and even numbers of column
            ri = 2 - mod(r,2);

            eta_11(r,c) = 2 * J_vel(:,r).' ...
                * (J_omega(:,3*(c-1)+1:3*(c-1)+3) * p_q(:,h)) ...
                + (J_omega(:,3*(r-1)+1:3*(r-1)+3) * p_q(:,h)).' ...
                * (J_omega(:,3*(c-1)+1:3*(c-1)+3) * p_q(:,h));

            if c<=2
            eta_12(r,c) = (J_vel(:,r) + J_omega(:,3*(r-1)+1:3*(r-1)+3) * p).' ...
                * p_qq(3*(h-1)+1:3*(h-1)+3, ci)...
                + (J_omega(:,3*(r-1)+1:3*(r-1)+3) * p_q(:,h)).' * p_q(:,ci);

            if r <=2
                eta_22(r,c) = p_qq(3*(h-1)+1:3*(h-1)+3, ri).' * p_q(:,ci);
            end

            end
        end
    end

end

if d == 0
    eta_22 = zeros(2,2);
    for r=1:2
        for c=1:2
            eta_22(r,c) = p_qq(3*(h-1)+1:3*(h-1)+3, r).' * p_q(:,c);
        end
    end
    % if section 1, then only gamma_22 is present
    Mih = eta_22;
else
    % from section 2 onwards, both eta and gamma exist based on the h
    Mih = [eta_11, eta_12;
          eta_12.', eta_22];
end

if any(isnan(Mih), 'all') || any(isinf(Mih), 'all')
    Mih = eye(size(Mih));
end

end