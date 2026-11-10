function [desiredPressurePrediction_bar,nominalPressurePrediction_bar,desiredPressureDiagnostic]= ...
    rlResidualDesiredPressureBlock(inverseContext,desiredMotion, ...
    shadowPressureOptions,variableCommon,singlePassOptions) %#ok<INUSD>
%RLRESIDUALDESIREDPRESSUREBLOCK Paste into a discrete 0.05 s MATLAB Function.
% desiredMotion = [l2d;l3d;dl2d;dl3d;ddl2d;ddl3d], double 6x1, SI units.
% Feedforward pressure for DESIRED motion, with a frozen observed residual.
% One GRU evaluation at actual measured q/v, actual command and causal history.
% No recurrent-state update, inverse iterations or integration in this block.
% This is a single-pass inverse approximation, not an exact nonlinear solve.
% Initial installation: Scopes/logging only; existing valve/Kill wiring stays.
% Diagnostic output type: Bus: ResidualSinglePassDiagnosticsBus (unchanged).
% variableCommon is retained for port compatibility and is not used.
% Pressure options: 4=min(0), 5=max, 8=slew, 9=activation(0), 10=N/bar.
%#codegen
desiredPressurePrediction_bar=nan(3,1);nominalPressurePrediction_bar=nan(3,1);
desiredPressureDiagnostic=residualSinglePassDiagnostics();
assert(isequal(size(desiredMotion),[6 1]),'desiredMotion must be a 6x1 column.');
if ~inverseContext.valid,return;end
assert(all(isfinite(singlePassOptions))&& ...
    (singlePassOptions(1)==64||singlePassOptions(1)==128)&&singlePassOptions(2)>0);
assert(singlePassOptions(1)==inverseContext.units);
if ~all(isfinite([inverseContext.measured;inverseContext.pressure_bar; ...
        inverseContext.historyPre;inverseContext.featureDt_s;inverseContext.horizon_s; ...
        inverseContext.sampleTime_host_s;desiredMotion])) || ...
        inverseContext.featureDt_s<=0 || inverseContext.horizon_s<=0
    return;
end
assert(all(isfinite(shadowPressureOptions([4;5;8;9;10]))));
assert(shadowPressureOptions(4)==0,'Collection policy requires zero minimum pressure.');
assert(shadowPressureOptions(9)==0,'Use the as-trained zero-activation pressure map.');
maximum=shadowPressureOptions(5);slew=shadowPressureOptions(8);alpha=shadowPressureOptions(10);
assert(maximum>0&&slew>0&&alpha>0);
network=residualControllerNetwork(singlePassOptions(1));
assert(abs(alpha-network.forcePerBar)<1e-10*network.forcePerBar);
if any(inverseContext.pressure_bar<-1e-9)||any(inverseContext.pressure_bar>maximum+1e-9),return;end
physics=residualShadowPhysicsConfig();n=physics.nominal;
% Both branches target precisely the SAME desired trajectory.
qd=desiredMotion(1:2);vd=desiredMotion(3:4);ad=desiredMotion(5:6);
P=forceResidualPhysics(qd,vd,n);M=[P(1),P(2);P(2),P(3)];C=[P(4),P(5);P(6),P(7)];
target=M*ad+(C+n.D)*vd+n.K*qd+ ...
    armS_single_boundForce(qd,n.mu,n.lKbounds)+P(8:9);
if ~all(isfinite(target)),return;end
% History and GRU features remain ACTUAL observations, not desired motion.
[~,residual]=residualGruCell(network,[inverseContext.measured; ...
    inverseContext.pressure_bar;inverseContext.featureDt_s],inverseContext.historyPre);
desiredPressureDiagnostic.networkEvaluations=uint8(1);
request=target-residual;
[corrected,allocationLimited]=minimumPressure(request,alpha,maximum);
[nominal,nominalLimited]=minimumPressure(target,alpha,maximum);
step=slew*inverseContext.horizon_s;
p=min(max(corrected,inverseContext.pressure_bar-step),inverseContext.pressure_bar+step);
pn=min(max(nominal,inverseContext.pressure_bar-step),inverseContext.pressure_bar+step);
if ~all(isfinite([target;residual;corrected;nominal;p;pn])),return;end
B=alpha*[-1 1 0;-1 0 1];
desiredPressurePrediction_bar=p;nominalPressurePrediction_bar=pn;
desiredPressureDiagnostic.valid=true;
desiredPressureDiagnostic.allocationLimited=allocationLimited;
desiredPressureDiagnostic.slewLimited=any(abs(p-corrected)>1e-9);
desiredPressureDiagnostic.nominalAllocationLimited=nominalLimited;
desiredPressureDiagnostic.nominalSlewLimited=any(abs(pn-nominal)>1e-9);
desiredPressureDiagnostic.targetForce_N=target;desiredPressureDiagnostic.residualForce_N=residual;
desiredPressureDiagnostic.correctedForceRequest_N=request;
desiredPressureDiagnostic.frozenTotalForce_N=B*p+residual;
desiredPressureDiagnostic.allocationError_N=B*corrected-request;
desiredPressureDiagnostic.slewError_N=B*(p-corrected);
desiredPressureDiagnostic.frozenBalanceError_N=B*p+residual-target;
desiredPressureDiagnostic.nominalBalanceError_N=B*pn-target;
desiredPressureDiagnostic.frozenBalanceSatisfied=norm(desiredPressureDiagnostic.frozenBalanceError_N)<=singlePassOptions(2);
desiredPressureDiagnostic.residualInputPressure_bar=inverseContext.pressure_bar;
desiredPressureDiagnostic.correctedPressureBeforeSlew_bar=corrected;
desiredPressureDiagnostic.nominalPressureBeforeSlew_bar=nominal;
desiredPressureDiagnostic.pressureDifferenceBeforeSlew_bar=corrected-nominal;
desiredPressureDiagnostic.pressureDifference_bar=p-pn;
desiredPressureDiagnostic.sampleTime_host_s=inverseContext.sampleTime_host_s;
desiredPressureDiagnostic.horizon_s=inverseContext.horizon_s;
end
function [pressure,limited]=minimumPressure(force,alpha,maximum)
forceBar=force/alpha;p1=max([0;-forceBar(1);-forceBar(2)]);
pressure=max([p1;p1+forceBar(1);p1+forceBar(2)],0);
peak=max(pressure);limited=peak>maximum;
if limited,pressure=pressure*(maximum/peak);end
pressure(pressure<=1e-8)=0;
end
