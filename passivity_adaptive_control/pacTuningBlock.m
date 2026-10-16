function [adaptationGamma,parameterScale,parameterMaxRate, ...
    parameterLower,parameterUpper,parameterInitial,adaptationOptions, ...
    historyOptions,pressureOptions,supervisorOptions] = pacTuningBlock()
%#codegen
% One source block for the ten pacAdaptiveTunableBlock tuning inputs.
% Edit before a run. Initial estimates are used only at controller startup.
% Parameter order: k1 k2 k3 d2 d3 h1 h2 h3 b2 b3.

adaptationGamma = [0.2;0.2;0.2;0.0002;0.0002;0.5;0.5;0.5;1;1];
parameterScale = [1000;1000;1000;100;100;1;1;1;1;1];
parameterMaxRate = [200;200;200;10;10;0.2;0.2;0.2;0.2;0.2];
parameterLower = [50;50;0;1;1;0;0;0;-8;-8];
parameterUpper = [6000;6000;6000;250;250;10;10;10;8;8];
parameterInitial = [650;650;650;40;40;0;0;0;0;0];

adaptationOptions = [ ...
    0.005; ... % Filtered error s-norm deadband (m/s)
    0.02;  ... % Allocation residual tolerance (N)
    0.075];    % Maximum parameter update interval (s)

historyOptions = [ ...
    0.005; ... % Normalizing history length (m)
    1.0;   ... % History beta
    0.5;   ... % History gamma
    2.0];      % History state bound

pressureOptions = [ ...
    1.6; ...                  % Nominal common pressure (bar)
   -0.0983981693363845; ...   % Symmetric trim (bar)
   -0.0819816714043049; ...   % Antisymmetric trim (bar)
    0.2; ...                  % Minimum command pressure (bar)
    3.0; ...                  % Maximum command pressure (bar)
    2.5; ...                  % Symmetric correction cap (bar)
    2.5; ...                  % Antisymmetric correction cap (bar)
    2.0; ...                  % Command slew limit (bar/s)
    0.8; ...                  % Nominal activation threshold (bar)
    pi*0.0065^2*1e5];          % Force coefficient (N/bar)

supervisorOptions = [ ...
    10.0; ...   % Pressure enable start (s)
    15.0; ...   % Earliest feedback engagement (s)
    5.0; ...    % Feedback blend duration (s)
    0.10; ...   % Invalid measurement grace (s)
    0.075; ...  % Maximum measurement age (s)
    0.0005; ... % Maximum anchor speed (m/s)
    inf];      % Stop time (s): external enable controls stopping
end
