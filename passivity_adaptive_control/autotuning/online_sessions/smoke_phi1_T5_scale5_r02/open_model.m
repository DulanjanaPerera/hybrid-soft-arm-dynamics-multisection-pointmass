% Open/review only. Does not enable pressure or start the model.
sessionFolder=fileparts(mfilename('fullpath'));
repoFolder=fileparts(fileparts(fileparts(fileparts(sessionFolder))));
addpath(repoFolder,fullfile(repoFolder,'passivity_adaptive_control'), ...
    fullfile(repoFolder,'passivity_adaptive_control','autotuning'));
pacOnlineConfigureBuild();
open_system(fullfile(repoFolder,'PAC_autotuning_online_smoke_r02.slx'));
