function buildRoot=pacOnlineConfigureBuild()
%PACONLINECONFIGUREBUILD Short, process-local Simulink build/cache directories.
% Call in the MODEL MATLAB process. No model is started or compiled.
% Per-process isolation prevents concurrent MATLAB builds sharing a cache.
buildRoot=fullfile(tempdir,'pac_bld',sprintf('p%d',feature('getpid')));
cache=fullfile(buildRoot,'cache');code=fullfile(buildRoot,'code');
assert(numel(buildRoot)<100,'Temporary path too long; select a shorter buildRoot.');
Simulink.fileGenControl('set','CacheFolder',cache,'CodeGenFolder',code,'createDir',true);
fprintf('PAC Simulink cache: %s\nPAC code generation: %s\n',cache,code);
end
