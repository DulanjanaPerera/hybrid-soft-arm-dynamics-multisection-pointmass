function C = christoffelSymbol(n, M, dl)
% This function compute the C matrix using the christoffel symbol. The
% varible joint-space is given as dl. the dl is a vector concatinated of
% all the joint-space velocities. The input n represent the current section
% which is not used in here.
% 
% The M is a 3D matrix consisting of partial derivatives of Mass matrix
% w.r.t. the all the join-space varibles upto section n. The 3rd dimension
% represents the joint-space variable.
%   
%       [[dMass], variable]
% 
% Inputs:
%   n   : (just for the clarity) current section [constant] 
%   M   : 3D matrix with partial derivatives of M w.r.t. the joints-space
%         varibles upto the section [2n, 2n, h]
%   dl  : the respective joint-space varible velocity vector [2n,1] (m/s)
% 
% Output:
%   C   : current n section C matrix [2n,2n]

arguments (Input)
    n (1,1) double
    M (:,:,:) double 
    dl (:,1) double 
end

arguments (Output)
    C (:,:) double
end


[Row,Col,H] = size(M);
C = zeros(Row,Col);

for j=1:Row
    for k=1:Col
        for h=1:H
            C(j,k) = C(j,k) ...
                + (M(k,j,h) + M(k,h,j) - M(h,j,k)) * dl(h);
        end
        C(j,k) = 0.5 * C(j,k);
    end
end

end