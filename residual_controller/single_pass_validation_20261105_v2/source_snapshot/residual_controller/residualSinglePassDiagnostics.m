function d=residualSinglePassDiagnostics()
%RESIDUALSINGLEPASSDIAGNOSTICS Frozen-residual algebra, NOT inverse convergence.
%#codegen
d.valid=false;d.networkEvaluations=uint8(0);
d.allocationLimited=false;d.slewLimited=false;
d.nominalAllocationLimited=false;d.nominalSlewLimited=false;
d.frozenBalanceSatisfied=false;
d.targetForce_N=nan(2,1);d.residualForce_N=nan(2,1);
d.correctedForceRequest_N=nan(2,1);d.frozenTotalForce_N=nan(2,1);
d.allocationError_N=nan(2,1);d.slewError_N=nan(2,1);
d.frozenBalanceError_N=nan(2,1);d.nominalBalanceError_N=nan(2,1);
d.residualInputPressure_bar=nan(3,1);
d.correctedPressureBeforeSlew_bar=nan(3,1);
d.nominalPressureBeforeSlew_bar=nan(3,1);
d.pressureDifferenceBeforeSlew_bar=nan(3,1);d.pressureDifference_bar=nan(3,1);
d.sampleTime_host_s=NaN;d.horizon_s=NaN;
end
