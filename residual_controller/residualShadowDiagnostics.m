function d=residualShadowDiagnostics()
%RESIDUALSHADOWDIAGNOSTICS Fixed-size logging bus prototype.
%#codegen
d.status=uint8(0);d.accepted=false;d.advanced=false;d.ready=false;
d.historyReset=false;d.predictionValid=false;
d.scoreValid=false;d.retimedScoreValid=false;d.inverseContextValid=false;
d.qMeasured_m=nan(2,1);d.velocity_mPerS=nan(2,1);d.rawForce_N=nan(2,1);
d.errorAtMeasurement_m=nan(2,3);d.retimedError_m=nan(2,3);
d.sampleTime_host_s=NaN;d.predictionTarget_host_s=NaN;d.scoreOrigin_host_s=NaN;
d.dtPhysical_s=NaN;d.commandAge_sim_s=NaN;d.horizonMismatch_s=NaN;
d.pressureChangeSinceOrigin_bar=nan(3,1);d.velocityWarm_s=0;d.neuralWarm_s=0;
d.maximumTrainingZ=NaN;d.units=0;d.segment=uint32(0);
end
