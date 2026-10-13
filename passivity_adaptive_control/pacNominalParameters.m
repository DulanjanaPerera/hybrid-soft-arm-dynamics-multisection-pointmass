function [params, lambda, dampingGain] = pacNominalParameters(gravity_arm)
%PACNOMINALPARAMETERS Initial nominal model; gains are offline starting values.
% lambda: 1/s; dampingGain: N*s/m. No payload or adaptation.
%#codegen
params.N = 1;
params.L = 0.17411;
params.r = 0.013;
params.mi = 0.1;
params.g = gravity_arm(:);
params.K = 650*[2 1;1 2];
params.D = 40*eye(2);
params.mu = 2000;
params.lKbounds = [-0.035;0.035;1e6];
lambda = [200;200];
dampingGain = [20;20];
end
