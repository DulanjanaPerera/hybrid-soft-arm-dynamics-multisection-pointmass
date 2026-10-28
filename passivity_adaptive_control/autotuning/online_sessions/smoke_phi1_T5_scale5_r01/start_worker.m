% SECOND MATLAB process. Starts only the fixed-gain smoke-test worker.
sessionFolder=fileparts(mfilename('fullpath'));
repoFolder=fileparts(fileparts(fileparts(fileparts(sessionFolder))));
addpath(repoFolder,fullfile(repoFolder,'passivity_adaptive_control'), ...
    fullfile(repoFolder,'passivity_adaptive_control','autotuning'));
pacOnlineWorker(fullfile(sessionFolder,'pacOnlineSession.mat'));
