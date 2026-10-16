% addpath('BaseRotation')
% analyzeCommonPressureCircle('initial_test.mat');
%
% addpath('integral_term')
%
% analyzeCommonPressureCircle('initial_test.mat');

addpath('PAC_circle');
% results = analyzePacStaticSweep('PAC_circle/static_phi_lambda_200_damp_20_sweep.mat');

addpath('passivity_adaptive_control');
% % results = analyzePacAdaptive('passivity_adaptive_control/PAC_phi_0.2_static_reducedGains_engagedFalse.mat');
analyzePacAdaptive();
% analyzePacStaticSweep();
