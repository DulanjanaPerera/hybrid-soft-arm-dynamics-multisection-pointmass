% MODEL MATLAB process, while stopped; worker must run in OTHER instance.
sessionFolder=fileparts(mfilename('fullpath'));
repoFolder=fileparts(fileparts(fileparts(fileparts(sessionFolder))));
addpath(repoFolder,fullfile(repoFolder,'passivity_adaptive_control'), ...
    fullfile(repoFolder,'passivity_adaptive_control','autotuning'));
pacOnlineConnectionCheck(fullfile(sessionFolder,'pacOnlineSession.mat'));
