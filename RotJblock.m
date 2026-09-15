function dR = RotJblock(RotJ, localIdx)
% Extract 3x3 block from RotJ (3x6) corresponding to local DoF 1 or 2.
% RotJ = [dR/dq1, dR/dq2] with each block 3 columns (stacked 3x3).
c1 = 3*(localIdx-1)+1;
c2 = c1 + 2;
dR = RotJ(:, c1:c2);
end
