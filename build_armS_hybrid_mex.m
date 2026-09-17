function build_armS_hybrid_mex(buildFolder)
% Build hybrid distributed-plus-added-point-mass RHS with runtime masses/xi.
% Builds a native MEX for the current OS; requires MATLAB Coder and a C++ compiler.
if nargin<1 || isempty(buildFolder)
    buildFolder=fullfile(tempdir,'armS_hybrid_mex_build');
end
assert(ischar(buildFolder) || (isstring(buildFolder) && isscalar(buildFolder)), ...
    'buildFolder must be a folder path.');
buildFolder=char(buildFolder);
if isempty(mex.getCompilerConfigurations('C++','Selected'))
    error('No C++ MEX compiler selected. Run mex -setup C++ in MATLAB.');
end
sourceFolder=fileparts(mfilename('fullpath'));
oldFolder=pwd; oldPath=path;
cleanup=onCleanup(@() restoreEnvironment(oldFolder,oldPath)); %#ok<NASGU>
if ~isfolder(buildFolder), mkdir(buildFolder); end
addpath(sourceFolder,'-begin');
clear armS_hybrid_mex
cd(buildFolder);
args={0,zeros(12,1),.278,.013,.1*ones(3,1),[0;0;-9.81], ...
    2200*eye(6),100*eye(6),zeros(6,1),2000,[-.02;.02;1e6], ...
    [.015;.025;.01],[.2;.5;.9]};
cfg=coder.config('mex'); cfg.GenerateReport=false; cfg.LaunchReport=false;
cfg.TargetLang='C++'; cfg.EnableOpenMP=false;
cfg.InlineBetweenUserFunctions='Never';
fprintf('Building hybrid MEX in %s\n',buildFolder); drawnow;
timer=tic;
codegen('-v','-config',cfg,'armS_hybrid_entry','-args',args, ...
    '-d',fullfile(pwd,'codegen'),'-o','armS_hybrid_mex');
fprintf('Code generation and compilation finished in %.1f s\n',toc(timer));
rehash;
poses=[zeros(6,1),[-.008;.003;-.004;-.006;.002;-.005]];
for k=1:2
    input=args; input{2}=[poses(:,k);.01;-.02;.015;.003;-.005;.008];
    if k==2
        input{12}=[0;.04;0]; input{13}=[.1;.8;1];
        input{9}=[.2;-.2;0;0;0;0];
    end
    expected=armS_hybrid_entry(input{:}); actual=armS_hybrid_mex(input{:});
    err=max(abs(actual-expected)./(1+abs(expected)));
    fprintf('Hybrid MEX runtime-input check %d: %.3e\n',k,err);
    assert(all(isfinite(actual)) && err<1e-8);
end
binary=fullfile(pwd,['armS_hybrid_mex.',mexext]);
assert(isfile(binary),'Native MEX was not created: %s',binary);
clear armS_hybrid_mex
copyfile(binary,fullfile(sourceFolder,['armS_hybrid_mex.',mexext]),'f');
fprintf('Copied validated hybrid MEX to %s\n',sourceFolder);
end

function restoreEnvironment(oldFolder,oldPath)
cd(oldFolder); path(oldPath); rehash;
end
