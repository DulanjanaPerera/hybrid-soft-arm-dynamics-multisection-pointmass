function [M,C,G,dM] = armS_hybrid_core(q,dq,params)
% Distributed arm mass plus independent, physical backbone point masses.
% params.mi: existing distributed section masses (3x1).
% params.addedMass: added point masses, modules 1..3 (3x1, zeros allowed).
% params.addedXi: their backbone positions (3x1, normalized 0..1).
% The added masses use unit beta: ordinary translational kinetic energy.
assert(params.N==3 && numel(q)==6 && numel(dq)==6);
assert(isfield(params,'addedMass') && isfield(params,'addedXi'));
assert(isequal(size(params.addedMass),[3,1]) && ...
    all(isfinite(params.addedMass)) && all(params.addedMass>=0));
assert(isequal(size(params.addedXi),[3,1]) && ...
    all(isfinite(params.addedXi)) && ...
    all(params.addedXi>=0) && all(params.addedXi<=1));
[Md,Cd,Gd,dMd]=armS_standard_core(q,dq,params);
[Mp,Cp,Gp,dMp]=armS_core_N3_mex(0,reshape(q,2,3).', ...
    reshape(dq,2,3).',params.L,params.r,params.addedXi, ...
    params.addedMass,params.g,params.K,ones(3,3));
M=Md+Mp; C=Cd+Cp; G=Gd+Gp; dM=dMd+dMp;
end
