function options=residualShadowOptions()
%RESIDUALSHADOWOPTIONS Configure before a fresh run; observer only.
% units; prospective horizon; velocity tau; velocity warmup; GRU warmup;
% max physical gap; max sim gap; max command age (sim); horizon tolerance;
% coordinate guard; maximum pressure. Times seconds; lengths m; pressure bar.
%#codegen
options=[128;.05;.20;.60;1.0;.10;.075;.075;.005;.05;3.0];
end
