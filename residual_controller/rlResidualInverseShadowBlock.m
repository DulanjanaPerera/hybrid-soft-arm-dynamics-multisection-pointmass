function [residualPressurePrediction_bar,nominalPressurePrediction_bar,inverseDiagnostic]= ...
    rlResidualInverseShadowBlock(inverseContext,desiredAcceleration_mPerS2, ...
    shadowPressureOptions,variableCommon,shadowInverseOptions)
%RLRESIDUALINVERSESHADOWBLOCK Stateless diagnostic inverse. NO valve connection.
% Optional second block. History comes from forward observer; never commits it.
% Targets the fixed trained nominal LHS at measured q/v and supplied ddq.
% This is not PAC's adaptive/anchored control demand or a controller replacement.
%#codegen
residualPressurePrediction_bar=nan(3,1);nominalPressurePrediction_bar=nan(3,1);
inverseDiagnostic=residualShadowInverseDiagnostics();
if ~inverseContext.valid||~all(isfinite(desiredAcceleration_mPerS2)),return;end
network=residualControllerNetwork(inverseContext.units);residualValidateOptions(shadowInverseOptions);
assert(shadowInverseOptions(1)==inverseContext.units);
% Shared configuration validation only; no supervisor/anchor state is run.
cfg=residualPressureConfig(shadowPressureOptions,[0;0;1;.1;.1;1;Inf],variableCommon);
assert(abs(cfg.forcePerBar_N-network.forcePerBar)<1e-10*network.forcePerBar);
assert(cfg.activation_bar==0,'Shadow nominal reference is the training linear pressure map.');
physics=residualShadowPhysicsConfig();n=physics.nominal;
q=inverseContext.measured(1:2);v=inverseContext.measured(3:4);
P=forceResidualPhysics(q,v,n);M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
target=M*desiredAcceleration_mPerS2+(C+n.D)*v+n.K*q+ ...
    armS_single_boundForce(q,n.mu,n.lKbounds)+P(8:9);
step=cfg.slew_barPerS*inverseContext.horizon_s;
[nominal,~,~,nomLimited]=pacAllocatePressure(target,zeros(2,1),zeros(2,1),cfg);
nominalPressurePrediction_bar=min(max(nominal,inverseContext.pressure_bar-step),inverseContext.pressure_bar+step);
B=network.forcePerBar*[-1 1 0;-1 0 1];
[~,raw]=residualGruCell(network,[inverseContext.measured;inverseContext.pressure_bar;inverseContext.featureDt_s],inverseContext.historyPre);
[p,e,limited,iterations]=residualAllocatePressure(target,zeros(2,1),zeros(2,1),zeros(2,1), ...
    inverseContext.pressure_bar,target-raw,target,step,inverseContext.measured, ...
    inverseContext.featureDt_s,inverseContext.historyPre,network,cfg,zeros(2,1),shadowInverseOptions);
if ~all(isfinite([p;e;target])),return;end
residualPressurePrediction_bar=p;inverseDiagnostic.valid=true;
inverseDiagnostic.converged=norm(e)<=shadowInverseOptions(3);
inverseDiagnostic.limited=limited;
inverseDiagnostic.nominalLimited=nomLimited||any(abs(nominal-nominalPressurePrediction_bar)>1e-9);
inverseDiagnostic.iterations=iterations;inverseDiagnostic.targetForce_N=target;
inverseDiagnostic.predictedForce_N=target+e;inverseDiagnostic.inverseError_N=e;
inverseDiagnostic.nominalInverseError_N=B*nominalPressurePrediction_bar-target;
inverseDiagnostic.pressureDifference_bar=p-nominalPressurePrediction_bar;
inverseDiagnostic.sampleTime_host_s=inverseContext.sampleTime_host_s;
end
