function build_armS_stationary_base_mex(buildFolder)
% Build the stationary-base wrapper and compare MATLAB and MEX at runtime rotations.
if nargin < 1, buildFolder = fullfile(tempdir,'armS_stationary_base_mex_build'); end
sourceFolder = fileparts(mfilename('fullpath'));
oldFolder = pwd; oldPath = path;
cleanup = onCleanup(@() restoreEnvironment(oldFolder,oldPath)); %#ok<NASGU>
if ~isfolder(buildFolder), mkdir(buildFolder); end
addpath(sourceFolder,'-begin');
clear armS_stationary_base_mex
cd(buildFolder);
args = {0,zeros(12,1),0.278,0.013,0.1*ones(3,1),eye(3), ...
    [0;0;-9.81],2200*eye(6),100*eye(6),zeros(6,1),2000,[-0.02;0.02;1e6]};
cfg = coder.config('mex');
cfg.GenerateReport = false; cfg.LaunchReport = false;
cfg.TargetLang = 'C++'; cfg.EnableOpenMP = false;
cfg.InlineBetweenUserFunctions = 'Never';
buildClock = tic;
codegen('-v','-config',cfg,'armS_stationary_base_entry', ...
    '-args',args,'-d',fullfile(pwd,'codegen'),'-o','armS_stationary_base_mex');
fprintf('Code generation and compilation finished in %.1f s\n',toc(buildClock));
rehash;
q = [-.008;.003;-.004;-.006;.002;-.005];
args{2} = [q;.01;-.02;.015;.003;-.005;.008];
args{5} = [.08;.10;.12]; args{10} = [.1;-.1;0;0;0;0];
rotations = stationaryBaseTestRotations();
worst = 0;
for k=1:numel(rotations)
    args{6}=rotations{k};
    if k==numel(rotations), args{7}=[1.1;-2.2;-8.7]; end
    expected=armS_stationary_base_entry(args{:});
    actual=armS_stationary_base_mex(args{:});
    assert(all(isfinite(actual)),'MEX returned nonfinite values.');
    err=max(abs(actual-expected)./(1+abs(expected)));
    worst=max(worst,err);
    fprintf('Orientation %d MEX scaled RHS error %.3e\n',k,err);
    assert(err<1e-8,'MATLAB and MEX disagree.');
end
binary=fullfile(pwd,['armS_stationary_base_mex.',mexext]);
assert(isfile(binary),'Compiled MEX was not found.');
clear armS_stationary_base_mex
copyfile(binary,fullfile(sourceFolder,['armS_stationary_base_mex.',mexext]),'f');
fprintf('Stationary-base MEX passed; maximum scaled error %.3e\n',worst);
end

function restoreEnvironment(oldFolder,oldPath)
cd(oldFolder); path(oldPath); rehash;
end
