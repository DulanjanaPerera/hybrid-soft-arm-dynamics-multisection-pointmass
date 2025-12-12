clear
N = 2;
syms 'l' [N 3] 
syms T(l)
for n=1:N 
    T_n = HTM_sym([0; l(n,2); l(n,3)], 1, 0.278, 0.013);
    T = T * T_n;
end