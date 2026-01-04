function Mih = Mi_h(n, h, p, p_q, p_qq, J_vel, J_omega, H_vel, H_omega)
% This function derives the partial derivatives of M (MASS) matrix for the
% given joint-space varible (h). 

    d = 2*(n-1);

    % For h < 2(n-1)
    eta_11 = zeros(d,d);
    eta_12 = zeros(d,2);
    eta_22 = zeros(2,2);

    % for 2(n-1)<h
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

        Mih = eta_22;
    else

        Mih = [eta_11, eta_12;
              eta_12.', eta_22];
    end
end