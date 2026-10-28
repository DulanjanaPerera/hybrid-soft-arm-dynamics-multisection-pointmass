% SECOND MATLAB process. Use this r02 script, not the completed r01 attempt.
sessionFolder=fileparts(mfilename('fullpath'));
repoFolder=fileparts(fileparts(fileparts(fileparts(sessionFolder))));
addpath(repoFolder,fullfile(repoFolder,'passivity_adaptive_control'), ...
    fullfile(repoFolder,'passivity_adaptive_control','autotuning'));
pacOnlineWorker(fullfile(sessionFolder,'pacOnlineSession.mat'));
