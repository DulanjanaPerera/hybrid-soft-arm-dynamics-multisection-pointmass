function [singlePassPressurePrediction_bar,nominalPressurePrediction_bar,singlePassDiagnostic]= ...
    rlResidualSinglePassShadowBlock(inverseContext,desiredAcceleration_mPerS2, ...
    shadowPressureOptions,variableCommon,singlePassOptions)
%RLRESIDUALSINGLEPASSSHADOWBLOCK Paste into a separate discrete0.05s function.
% Observation ONLY. Outputs go to logging/Scopes, never valves/control/Kill.
% One GRU evaluation at ACTUAL current pressure. No inverse iteration or ODE.
% Uses measured q/v and supplied ddq, exactly as the iterative shadow block.
% singlePassOptions=[network width64/128; frozen-balance tolerance N].
%#codegen
singlePassPressurePrediction_bar=nan(3,1);nominalPressurePrediction_bar=nan(3,1);
singlePassDiagnostic=residualSinglePassDiagnostics();
if ~inverseContext.valid,return;end
assert(all(isfinite(singlePassOptions))&& ...
    (singlePassOptions(1)==64||singlePassOptions(1)==128)&&singlePassOptions(2)>0);
assert(singlePassOptions(1)==inverseContext.units);
if ~all(isfinite([inverseContext.measured;inverseContext.pressure_bar; ...
        inverseContext.historyPre;inverseContext.featureDt_s;inverseContext.horizon_s; ...
        inverseContext.sampleTime_host_s;desiredAcceleration_mPerS2])) || ...
        inverseContext.featureDt_s<=0 || inverseContext.horizon_s<=0
    return;
end
cfg=residualPressureConfig(shadowPressureOptions,[0;0;1;.1;.1;1;Inf],variableCommon);
network=residualControllerNetwork(singlePassOptions(1));
assert(abs(cfg.forcePerBar_N-network.forcePerBar)<1e-10*network.forcePerBar);
assert(cfg.activation_bar==0,'Single-pass shadow uses the as-trained linear pressure reference.');
if any(inverseContext.pressure_bar<cfg.minimum_bar-1e-9)|| ...
        any(inverseContext.pressure_bar>cfg.maximum_bar+1e-9)
    return;
end
physics=residualShadowPhysicsConfig();n=physics.nominal;
q=inverseContext.measured(1:2);v=inverseContext.measured(3:4);
P=forceResidualPhysics(q,v,n);M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
target=M*desiredAcceleration_mPerS2+(C+n.D)*v+n.K*q+ ...
    armS_single_boundForce(q,n.mu,n.lKbounds)+P(8:9);
% Evaluate once, using immutable pre-current history. Discard next hidden state;
% the forward observer alone owns the one-update-per-observation history.
[~,residual]=residualGruCell(network,[inverseContext.measured; ...
    inverseContext.pressure_bar;inverseContext.featureDt_s],inverseContext.historyPre);
singlePassDiagnostic.networkEvaluations=uint8(1);
request=target-residual;
[corrected,~,~,allocationLimited]=pacAllocatePressure(request,zeros(2,1),zeros(2,1),cfg);
[nominal,~,~,nominalLimited]=pacAllocatePressure(target,zeros(2,1),zeros(2,1),cfg);
step=cfg.slew_barPerS*inverseContext.horizon_s;
p=min(max(corrected,inverseContext.pressure_bar-step),inverseContext.pressure_bar+step);
pn=min(max(nominal,inverseContext.pressure_bar-step),inverseContext.pressure_bar+step);
if ~all(isfinite([target;residual;corrected;nominal;p;pn])),return;end
B=network.forcePerBar*[-1 1 0;-1 0 1];
singlePassPressurePrediction_bar=p;nominalPressurePrediction_bar=pn;
singlePassDiagnostic.valid=true;
singlePassDiagnostic.allocationLimited=allocationLimited;
singlePassDiagnostic.slewLimited=any(abs(p-corrected)>1e-9);
singlePassDiagnostic.nominalAllocationLimited=nominalLimited;
singlePassDiagnostic.nominalSlewLimited=any(abs(pn-nominal)>1e-9);
singlePassDiagnostic.targetForce_N=target;singlePassDiagnostic.residualForce_N=residual;
singlePassDiagnostic.correctedForceRequest_N=request;
singlePassDiagnostic.frozenTotalForce_N=B*p+residual;
singlePassDiagnostic.allocationError_N=B*corrected-request;
singlePassDiagnostic.slewError_N=B*(p-corrected);
singlePassDiagnostic.frozenBalanceError_N=B*p+residual-target;
singlePassDiagnostic.nominalBalanceError_N=B*pn-target;
singlePassDiagnostic.frozenBalanceSatisfied= ...
    norm(singlePassDiagnostic.frozenBalanceError_N)<=singlePassOptions(2);
singlePassDiagnostic.residualInputPressure_bar=inverseContext.pressure_bar;
singlePassDiagnostic.correctedPressureBeforeSlew_bar=corrected;
singlePassDiagnostic.nominalPressureBeforeSlew_bar=nominal;
singlePassDiagnostic.pressureDifferenceBeforeSlew_bar=corrected-nominal;
singlePassDiagnostic.pressureDifference_bar=p-pn;
singlePassDiagnostic.sampleTime_host_s=inverseContext.sampleTime_host_s;
singlePassDiagnostic.horizon_s=inverseContext.horizon_s;
end
