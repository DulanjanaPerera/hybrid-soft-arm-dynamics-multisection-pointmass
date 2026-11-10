function network=residualControllerNetwork(units)
%RESIDUALCONTROLLERNETWORK Fixed-size code-generation constants; run packer first.
% Set units64/128 before a fresh run. Never replace weights during a trial.
%#codegen
persistent n64 n128
if isempty(n64)
    a=coder.load('residualNetwork64_v1.mat','network');n64=a.network;
    b=coder.load('residualNetwork128_v1.mat','network');n128=b.network;
end
assert(units==64||units==128);
if units==64,network=n64;else,network=n128;end
end
