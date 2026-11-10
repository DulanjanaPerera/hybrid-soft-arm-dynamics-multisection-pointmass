function [singlePassPressurePrediction_bar,nominalPressurePrediction_bar,singlePassDiagnostic]= ...
    rlResidualSinglePassCollectionShadowBlock(inverseContext,desiredAcceleration_mPerS2, ...
    shadowPressureOptions,variableCommon,singlePassOptions) %#ok<INUSD>
%RLRESIDUALSINGLEPASSCOLLECTIONSHADOWBLOCK Paste into discrete0.05s function.
% Observation ONLY: outputs to logging/Scopes, never valves/control/Kill.
% Same5inputs/3outputs/bus as rlResidualSinglePassShadowBlock.
% variableCommon is a retained compatibility input; NOT used by this policy.
% Pressure options used:4minimum(0),5maximum,8slew,9activation(0),10force/bar.
% Common pressure/trim/s-a caps(1:3,6:7) are not part of collection allocation.
% One GRU evaluation at actual pressure; no inverse iterations/state commits.
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
assert(all(isfinite(shadowPressureOptions([4;5;8;9;10]))));
assert(shadowPressureOptions(4)==0,'Collection policy requires zero minimum pressure.');
assert(shadowPressureOptions(9)==0,'Collection policy uses the as-trained zero-activation map.');
maximum=shadowPressureOptions(5);slew=shadowPressureOptions(8);alpha=shadowPressureOptions(10);
assert(maximum>0&&slew>0&&alpha>0);
network=residualControllerNetwork(singlePassOptions(1));
assert(abs(alpha-network.forcePerBar)<1e-10*network.forcePerBar);
if any(inverseContext.pressure_bar<-1e-9)||any(inverseContext.pressure_bar>maximum+1e-9),return;end
physics=residualShadowPhysicsConfig();n=physics.nominal;
q=inverseContext.measured(1:2);v=inverseContext.measured(3:4);
P=forceResidualPhysics(q,v,n);M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
target=M*desiredAcceleration_mPerS2+(C+n.D)*v+n.K*q+ ...
    armS_single_boundForce(q,n.mu,n.lKbounds)+P(8:9);
[~,residual]=residualGruCell(network,[inverseContext.measured; ...
    inverseContext.pressure_bar;inverseContext.featureDt_s],inverseContext.historyPre);
singlePassDiagnostic.networkEvaluations=uint8(1);
request=target-residual;
% Same minimum-nonnegative conversion as the collection inverse. No preload.
[corrected,allocationLimited]=minimumPressure(request,alpha,maximum);
[nominal,nominalLimited]=minimumPressure(target,alpha,maximum);
step=slew*inverseContext.horizon_s;
p=min(max(corrected,inverseContext.pressure_bar-step),inverseContext.pressure_bar+step);
pn=min(max(nominal,inverseContext.pressure_bar-step),inverseContext.pressure_bar+step);
if ~all(isfinite([target;residual;corrected;nominal;p;pn])),return;end
B=alpha*[-1 1 0;-1 0 1];
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
singlePassDiagnostic.frozenBalanceSatisfied=norm(singlePassDiagnostic.frozenBalanceError_N)<=singlePassOptions(2);
singlePassDiagnostic.residualInputPressure_bar=inverseContext.pressure_bar;
singlePassDiagnostic.correctedPressureBeforeSlew_bar=corrected;
singlePassDiagnostic.nominalPressureBeforeSlew_bar=nominal;
singlePassDiagnostic.pressureDifferenceBeforeSlew_bar=corrected-nominal;
singlePassDiagnostic.pressureDifference_bar=p-pn;
singlePassDiagnostic.sampleTime_host_s=inverseContext.sampleTime_host_s;
singlePassDiagnostic.horizon_s=inverseContext.horizon_s;
end
function [pressure,limited]=minimumPressure(force,alpha,maximum)
forceBar=force/alpha;p1=max([0;-forceBar(1);-forceBar(2)]);
pressure=max([p1;p1+forceBar(1);p1+forceBar(2)],0);
peak=max(pressure);limited=peak>maximum;
if limited,pressure=pressure*(maximum/peak);end
pressure(pressure<=1e-8)=0;
end
