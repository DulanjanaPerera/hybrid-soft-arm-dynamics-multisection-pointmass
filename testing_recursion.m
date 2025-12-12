clear
N = 2;
syms 'l' [N 3] 
% syms T(l)
T = eye(4); 
for n=1:N 
    [T_n, ~, ~] = HTM_sym([0, l(n,2), l(n,3)], 1, 0.278, 0.013);
    T = T * T_n;
end

l_num = zeros(N,3);   % initialize
l_num(1,2) = 0.01;
l_num(1,3) = 0.01;
l_num(2,2) = 0.01;
l_num(2,3) = 0.01;

% Evaluate numerically
T_num = double(subs(T, l, l_num));

TT=HTM(l_num(1,:))*HTM(l_num(2,:));

eT = T_num-TT;
