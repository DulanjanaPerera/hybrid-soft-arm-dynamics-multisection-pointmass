function build_yost
%BUILD_YOST Compile the Windows C++ MEX S-function for this MATLAB release.
root = fileparts(mfilename('fullpath'));
assert(ispc && strcmp(computer('arch'),'win64'), 'Windows 64-bit MATLAB is required.');
compiler = mex.getCompilerConfigurations('C++','Selected');
if isempty(compiler)
    error('Yost:Compiler','No selected C++ compiler. Run mex -setup C++ first.');
end
clear yost_imu_sfun
mex('-R2018a','-outdir',root,'-output','yost_imu_sfun', ...
    ['-I' root],['-I' fullfile(matlabroot,'simulink','include')], ...
    fullfile(root,'yost_imu_sfun.cpp'));
addpath(root);
fprintf('Built %s using %s.\n', fullfile(root,['yost_imu_sfun.' mexext]), compiler.Name);
end
