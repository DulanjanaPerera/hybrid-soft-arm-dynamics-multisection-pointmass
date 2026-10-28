% Run in a SECOND MATLAB process. Never starts a model or operates hardware.
sessionFolder=fileparts(mfilename('fullpath'));
repoFolder=fileparts(fileparts(fileparts(fileparts(sessionFolder))));
addpath(repoFolder,fullfile(repoFolder,'passivity_adaptive_control'), ...
    fullfile(repoFolder,'passivity_adaptive_control','autotuning'));
pacOnlineWorker(fullfile(sessionFolder,'pacOnlineSession.mat'));
