function k_N_m = stiffnessFromPhiBlock(phi_rad)
%#codegen
% Scalar P1 effective stiffness from the 2026-10-02 fitted data (N/m).
% Input may be current phi or externally delayed phi, in radians.
% Endpoint extension: <=0.09425 rad -> 4132; >=2.05594 rad -> 953.333.
% Invalid input uses the near-straight endpoint (not a sensor-validity gate).
breaks = [0.094252336020426999 0.31719551885716302 ...
          0.84254655402063705 1.3628983100343199 ...
          1.7568402609821301 2.0559364288348001];
coefs = [-301.18295127845971 -90.779050026315701 -3508.295263988015 4132; ...
          8648.306884889811 -4672.8276757423591 -3593.6820922138418 3342; ...
          222.76060985869111 882.71017208803505 -1342.7968944369923 1418.3333333333301; ...
          341.14621515352195 107.09908322703643 -243.20935379879649 990; ...
          254.06580436465293 166.20780898136815 0 931.66666666666697];
if ~isfinite(phi_rad) || phi_rad <= breaks(1)
    k_N_m = 4132;
    return;
end
if phi_rad >= breaks(6)
    k_N_m = 953.333333333333;
    return;
end
segment = 1;
for j = 2:5
    if phi_rad >= breaks(j), segment = j; end
end
x = phi_rad-breaks(segment);
k_N_m = ((coefs(segment,1)*x+coefs(segment,2))*x ...
           +coefs(segment,3))*x+coefs(segment,4);
end
