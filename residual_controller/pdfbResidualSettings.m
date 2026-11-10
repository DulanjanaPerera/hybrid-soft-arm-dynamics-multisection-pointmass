function [pressureOptions,supervisorOptions,residualOptions]=pdfbResidualSettings()
%PDFBRESIDUALSETTINGS Read-only audit of saved FBL_commonP_circle_BRot.slx.
% Check live block settings before installing. Does not set gains/references.
% The saved inverse uses LINEAR pressure differences, so activation=0 here.
% Model-off/on in the NEW wrapper use the same PAC pressure supervisor/allocator.
%#codegen
pressureOptions=[1.6;-.0983981693363845;-.0819816714043049;.2;3;2.5;2.5;2;0;pi*.0065^2*1e5];
% Saved feedback starts15s/blends5s/stops220s, with no speed gate at anchor.
supervisorOptions=[10;15;5;.10;.075;1e9;220];
residualOptions=residualControllerOptions();
end
