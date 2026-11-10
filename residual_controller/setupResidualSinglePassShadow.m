function [singlePassOptions,ResidualSinglePassDiagnosticsBus]=setupResidualSinglePassShadow(root)
%SETUPRESIDUALSINGLEPASSSHADOW Paths/packet verification and NEW bus only.
% Does not assign/reset the existing forward block's base-workspace options.
% No model/device I/O. Run alongside existing setupResidualShadow assignments.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
shadowOptions=setupResidualShadow(root);
singlePassOptions=[shadowOptions(1);.02];
prototype=residualSinglePassDiagnostics();names=fieldnames(prototype);
elements=Simulink.BusElement.empty;
for j=1:numel(names)
    v=prototype.(names{j});e=Simulink.BusElement;e.Name=names{j};
    e.Dimensions=size(v);e.DataType=class(v);if islogical(v),e.DataType='boolean';end
    e.Complexity='real';elements(j)=e; %#ok<AGROW>
end
ResidualSinglePassDiagnosticsBus=Simulink.Bus;ResidualSinglePassDiagnosticsBus.Elements=elements;
end
