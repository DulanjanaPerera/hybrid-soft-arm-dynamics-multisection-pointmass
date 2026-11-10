function [residualOptions,pdfbPressureOptions,pdfbSupervisorOptions,ResidualControllerDiagnosticsBus]=setupResidualControllers(root)
%SETUPRESIDUALCONTROLLERS Paths, frozen packets and logging bus only; no SLX/I/O.
% Assign all four outputs in the BASE workspace before compiling a model.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
addpath(root,fullfile(root,'residual_controller'),fullfile(root,'passivity_adaptive_control'));
prepareResidualControllerNetworks(root);
[pdfbPressureOptions,pdfbSupervisorOptions,residualOptions]=pdfbResidualSettings();
prototype=residualDiagnostics();names=fieldnames(prototype);elements=Simulink.BusElement.empty;
for k=1:numel(names)
    value=prototype.(names{k});element=Simulink.BusElement;
    element.Name=names{k};element.Dimensions=size(value);element.DataType=class(value);
    if islogical(value),element.DataType='boolean';end
    element.Complexity='real';elements(k)=element;
end
ResidualControllerDiagnosticsBus=Simulink.Bus;ResidualControllerDiagnosticsBus.Elements=elements;
fprintf('Residual paths/weights/bus prepared. No model opened or run. Defaults128units; keep enableResidual=false for initial wiring checks.\n');
end
