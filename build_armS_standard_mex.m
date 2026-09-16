function build_armS_standard_mex(buildFolder)
% Build outside OneDrive, verify numerical agreement, copy MEX beside source.
% Run: build_armS_standard_mex
if nargin < 1, buildFolder = 'C:\MATLAB_build\armStd3'; end
sourceFolder = fileparts(mfilename('fullpath'));
oldFolder = pwd;
oldPath = path;
cleanup = onCleanup(@() restoreEnvironment(oldFolder,oldPath)); %#ok<NASGU>
if ~isfolder(buildFolder), mkdir(buildFolder); end
addpath(sourceFolder,'-begin');
clear armS_standard_mex
cd(buildFolder);
% Size/type examples only: these VALUES are not baked into the binary.
args = {0,zeros(12,1),0.278,0.013,0.1*ones(3,1),[0;0;-9.81], ...
    2200*eye(6),100*eye(6),zeros(6,1),2000,[-0.02;0.02;1e6]};
cfg = coder.config('mex');
cfg.GenerateReport = false;
cfg.LaunchReport = false;
cfg.TargetLang = 'C++';
cfg.EnableOpenMP = false;
cfg.InlineBetweenUserFunctions = 'Never';

fprintf('Starting code generation and compilation in %s\n', pwd);
drawnow;
buildClock = tic;

codegen('-v', '-config', cfg, 'armS_standard_entry', ...
    '-args', args, ...
    '-d', fullfile(pwd, 'codegen'), ...
    '-o', 'armS_standard_mex');

fprintf('Code generation and compilation finished in %.1f s\n', ...
    toc(buildClock));
rehash;
% Test representative poses with nonzero velocities and a changed mass/input.
poses = [zeros(6,1),[-.001;-.001;-.001;-.001;-1e-6;-1e-6], ...
    [-.008;.003;-.004;-.006;.002;-.005]];
worst = 0;
for k = 1:3
    testArgs = args;
    testArgs{2} = [poses(:,k); .01;-.02;.015;.003;-.005;.008];
    if k == 3
        testArgs{5} = [0.08;0.10;0.12];
        testArgs{9} = [0.1;-0.1;0;0;0;0];
    end
    expected = armS_standard_entry(testArgs{:});
    actual = armS_standard_mex(testArgs{:});
    assert(all(isfinite(actual)),'MEX returned nonfinite values.');
    err = max(abs(actual-expected)./(1+abs(expected)));
    worst = max(worst,err);
    fprintf('MEX check %d: scaled maximum error %.3e\n',k,err);
    assert(err < 1e-8,'MEX and MATLAB disagree; binary was not copied.');
end
binary = fullfile(pwd,['armS_standard_mex.',mexext]);
assert(isfile(binary),'Cannot find compiled binary: %s',binary);
clear armS_standard_mex
copyfile(binary,fullfile(sourceFolder,['armS_standard_mex.',mexext]),'f');
fprintf('Build and checks passed (maximum error %.3e).\n',worst);
fprintf('MEX copied to: %s\n',sourceFolder);
end

function restoreEnvironment(oldFolder,oldPath)
cd(oldFolder);
path(oldPath);
rehash;
end
