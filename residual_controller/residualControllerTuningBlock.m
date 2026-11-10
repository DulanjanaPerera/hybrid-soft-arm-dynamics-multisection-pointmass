function [enableResidual,residualOptions]=residualControllerTuningBlock()
%RESIDUALCONTROLLERTUNINGBLOCK Optional paste-ready source for the two new ports.
% Change only before a fresh run. true activates the frozen learned model.
%#codegen
enableResidual=false;
residualOptions=residualControllerOptions();
end
