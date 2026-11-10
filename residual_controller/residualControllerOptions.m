function options=residualControllerOptions()
%RESIDUALCONTROLLEROPTIONS Structural settings; configure before a fresh run.
% units; max GN iterations; force tolerance N; valid warmup s; max gap s;
% pressure-projection finite difference N; max GN force step N; nominal period s.
%#codegen
options=[128;6;.02;1.0;.10;.001;20;.05];
end
