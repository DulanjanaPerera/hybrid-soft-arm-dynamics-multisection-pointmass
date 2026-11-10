function [shadowOptions,ResidualShadowDiagnosticsBus,ResidualShadowContextBus, ...
    shadowPressureOptions,shadowInverseOptions,ResidualShadowInverseDiagnosticsBus]=setupResidualShadow(root)
%SETUPRESIDUALSHADOW Paths, verified frozen packets and typed buses ONLY.
% Assign exact output names in the MATLAB BASE workspace. No SLX or device I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
addpath(root,fullfile(root,'residual_controller'),fullfile(root,'passivity_adaptive_control'), ...
    fullfile(root,'residual_dynamics_learning','force_model'), ...
    fullfile(root,'residual_dynamics_learning','force_model','streaming'));
% Verify the installed fixed packets/source hashes without rewriting weights.
info=readtable(fullfile(root,'residual_controller','NETWORKS.csv'),'TextType','string', ...
    'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve');
for k=1:height(info)
    packet=fullfile(root,'residual_controller',sprintf('residualNetwork%d_v1.mat',info.units(k)));
    assert(strcmp(sgHash(packet),info.packedSHA256(k)),'Frozen packet hash mismatch.');
end
shadowOptions=residualShadowOptions();shadowInverseOptions=residualControllerOptions();
% Diagnostic-only hypothetical allocator: raw training pressure map,0..3bar.
% Configure before testing. These are not changes to the running RL commands.
shadowPressureOptions=[1.6;0;0;0;3;2.5;2.5;6;0;pi*.0065^2*1e5];
ResidualShadowDiagnosticsBus=makeBus(residualShadowDiagnostics());
ResidualShadowContextBus=makeBus(residualShadowContext());
ResidualShadowInverseDiagnosticsBus=makeBus(residualShadowInverseDiagnostics());
fprintf('Shadow paths/packets/buses ready. No model opened, edited or run.\n');
end
function bus=makeBus(prototype)
names=fieldnames(prototype);elements=Simulink.BusElement.empty;
for j=1:numel(names)
    value=prototype.(names{j});e=Simulink.BusElement;e.Name=names{j};
    e.Dimensions=size(value);e.DataType=class(value);if islogical(value),e.DataType='boolean';end
    e.Complexity='real';elements(j)=e; %#ok<AGROW>
end
bus=Simulink.Bus;bus.Elements=elements;
end
