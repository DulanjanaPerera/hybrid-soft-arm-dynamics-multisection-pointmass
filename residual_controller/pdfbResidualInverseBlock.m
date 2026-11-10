function [pressure_bar,tauDesired_N,tauAvailable_N,limited,modelValid, ...
    MassForce_N,CoriolisForce_N,GnStiff_N,status,residualDiagnostic] = ...
    pdfbResidualInverseBlock(measured,desiredAcceleration,orientationValid,t,gravity,baseOrientationValid, ...
    measurementTime_s,enablePressure,enableResidual,enableCommonPressureAdaptation,pressureOptions,supervisorOptions,residualOptions)
%PDFBRESIDUALINVERSEBLOCK Paste-ready model/inverse; outer PDFB/PID unchanged.
% First6 inputs/first8 outputs match saved FBL_commonP_circle_BRot inverse.
% Append physical measurement time/enable, residual switch/configuration.
% Fresh model-off trials in THIS wrapper are the paired nominal baseline.
% The old embedded bound/allocator is not claimed bit-identical to this pair.
%#codegen
network=residualControllerNetwork(residualOptions(1));residualValidateOptions(residualOptions);
cfg=residualPressureConfig(pressureOptions,supervisorOptions,enableCommonPressureAdaptation);
persistent ps ns
if isempty(ps),ps=residualPressureInitialState();ns=residualNeuralInitialState();end
tauDesired_N=zeros(2,1);MassForce_N=zeros(2,1);CoriolisForce_N=zeros(2,1);GnStiff_N=zeros(2,1);
good=orientationValid&&baseOrientationValid&&all(isfinite(measured(1:4)))&&all(isfinite(desiredAcceleration))&&all(isfinite(gravity));
modelValid=false;shift=zeros(2,1);
if good
    [params,~,~]=pacNominalParameters(gravity);q=measured(1:2);dq=measured(3:4);
    [M,C,G]=armS_single_core(q,dq,params);bound=armS_single_boundForce(q,params.mu,params.lKbounds);
    MassForce_N=M*desiredAcceleration;CoriolisForce_N=(C+params.D)*dq;GnStiff_N=G+params.K*q+bound;
    tauDesired_N=MassForce_N+CoriolisForce_N+GnStiff_N;
    modelValid=all(isfinite(tauDesired_N))&&all(isfinite(M(:)))&&M(1,1)>0&&det(M)>0;
    if enableResidual&&modelValid,shift=residualGravityShift(measured,gravity);end
end
if ~modelValid,tauDesired_N(:)=0;MassForce_N(:)=0;CoriolisForce_N(:)=0;GnStiff_N(:)=0;end
[pressure_bar,tauAvailable_N,~,limited,status,~,~,ps,ns,residualDiagnostic]=residualPressureStep( ...
    tauDesired_N,GnStiff_N,measured,modelValid,measurementTime_s,t,enablePressure,cfg,ps,ns,network,residualOptions,enableResidual,shift);
end
