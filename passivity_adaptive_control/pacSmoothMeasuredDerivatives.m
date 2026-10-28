function [l,dl,ddl]=pacSmoothMeasuredDerivatives(lengths3_m,orientationValid,pollStart_s)
%PACSMOOTHMEASUREDDERIVATIVES Paste-ready replacement for measured derivatives.
%#codegen
% Same 3 inputs/3 outputs as filteredDerivatives in the saved manual model.
% Candidate only: offline validated, not yet physically validated.
% Two equal continuous poles smooth raw velocity; exact ZOH state propagation.
% tauV=.15 gives approximately .30 s low-frequency derivative-filter delay.
% Acceleration is a smoothed derivative of velocity, not an independent sensor.
tauV_s=.15;
tauA_s=.35;
maxGap_s=.15;
persistent lastL lastT stage1 velocity acceleration
if isempty(lastT)
    lastL=nan(2,1);lastT=-1;
    stage1=zeros(2,1);velocity=zeros(2,1);acceleration=zeros(2,1);
end
l=lastL;dl=velocity;ddl=acceleration;
if ~orientationValid||~all(isfinite(lengths3_m))||~isfinite(pollStart_s),return;end
q=[lengths3_m(2);lengths3_m(3)];
if lastT<0
    lastL=q;lastT=pollStart_s;l=q;return
end
dt=pollStart_s-lastT;
if dt<=1e-6,return;end % Ignore repeated/backward timestamps.
if dt>maxGap_s
    % Re-anchor after a long gap; do not differentiate the missing interval.
    % Retain derivative states on this frame; supervisor still owns validity.
    lastL=q;lastT=pollStart_s;l=q;return
end
rawVelocity=(q-lastL)/dt;
h=dt/tauV_s;a=exp(-h);
oldStage1=stage1;oldVelocity=velocity;
stage1=a*oldStage1+(1-a)*rawVelocity;
velocity=a*oldVelocity+h*a*oldStage1+(1-a-h*a)*rawVelocity;
rawAcceleration=(velocity-oldVelocity)/dt;
alphaA=1-exp(-dt/tauA_s);
acceleration=acceleration+alphaA*(rawAcceleration-acceleration);
lastL=q;lastT=pollStart_s;
l=lastL;dl=velocity;ddl=acceleration;
end
