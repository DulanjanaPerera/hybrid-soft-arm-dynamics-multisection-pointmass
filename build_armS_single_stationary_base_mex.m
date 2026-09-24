function build_armS_single_stationary_base_mex(buildFolder)
% Build and verify the one-section stationary-orientation MEX.
if nargin<1
    buildFolder=fullfile(tempdir,'armS_single_stationary_base_mex_build');
end
sourceFolder=fileparts(mfilename('fullpath'));
oldFolder=pwd; oldPath=path;
cleanup=onCleanup(@() restoreEnvironment(oldFolder,oldPath)); %#ok<NASGU>
if ~isfolder(buildFolder), mkdir(buildFolder); end
addpath(sourceFolder,'-begin');
clear armS_single_stationary_base_mex
cd(buildFolder);
args={0,zeros(4,1),.278,.013,.1,eye(3),[0;0;9.81], ...
    2200*eye(2),600*eye(2),zeros(2,1),2000,[-.02;.02;1e6]};
cfg=coder.config('mex');
cfg.GenerateReport=false; cfg.LaunchReport=false;
cfg.TargetLang='C++'; cfg.EnableOpenMP=false;
cfg.InlineBetweenUserFunctions='Never';
buildClock=tic;
codegen('-v','-config',cfg,'armS_single_stationary_base_entry', ...
    '-args',args,'-d',fullfile(pwd,'codegen'), ...
    '-o','armS_single_stationary_base_mex');
fprintf('Single-module code generation and compilation: %.1f s\n', ...
    toc(buildClock));
rehash;

rotations=stationaryBaseTestRotations();
poses=[zeros(2,1),[-.001;-.001],[-.008;.003]];
worst=0;
for p=1:size(poses,2)
    args{2}=[poses(:,p);.01;-.02];
    args{5}=.08+.01*p;
    args{10}=[.1;-.1];
    for k=1:numel(rotations)
        args{6}=rotations{k};
        if k==numel(rotations), args{7}=[1.1;-2.2;-8.7]; end
        expected=armS_single_stationary_base_entry(args{:});
        actual=armS_single_stationary_base_mex(args{:});
        assert(all(isfinite(actual)),'MEX returned nonfinite values.');
        err=max(abs(actual-expected)./(1+abs(expected)));
        worst=max(worst,err);
        assert(err<1e-8,'Single-module MATLAB/MEX mismatch.');
    end
end
binary=fullfile(pwd,['armS_single_stationary_base_mex.',mexext]);
assert(isfile(binary),'Compiled single-module MEX was not found.');
clear armS_single_stationary_base_mex
copyfile(binary,fullfile(sourceFolder, ...
    ['armS_single_stationary_base_mex.',mexext]),'f');
fprintf('Single-module MEX passed %d checks; worst scaled RHS error %.3e\n', ...
    size(poses,2)*numel(rotations),worst);
end

function restoreEnvironment(oldFolder,oldPath)
cd(oldFolder); path(oldPath); rehash;
end
