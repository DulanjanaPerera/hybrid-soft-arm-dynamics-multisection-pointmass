function w = vee(S)
% Map a 3x3 skew-symmetric matrix S to vector w such that S = [w]x.
% Assumes S is (approximately) skew: S = -S^T.
w = [ S(3,2); S(1,3); S(2,1) ];
end
