function d=residualShadowInverseDiagnostics()
%#codegen
d.valid=false;d.converged=false;d.limited=false;d.nominalLimited=false;d.iterations=uint8(0);
d.targetForce_N=nan(2,1);d.predictedForce_N=nan(2,1);d.inverseError_N=nan(2,1);
d.nominalInverseError_N=nan(2,1);d.pressureDifference_bar=nan(3,1);
d.sampleTime_host_s=NaN;
end
