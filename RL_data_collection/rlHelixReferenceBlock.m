function [theta,phi,trialId,phaseId,cycleId,direction,done,testEnd_s] = rlHelixReferenceBlock(t)
%#codegen
% Requested theta/phi: outward then inward spiral/helix-like arm-tip path.
% Actual XYZ follows arm geometry; this is not a constant-pitch Cartesian helix.
peakPhi_rad=2.5;
circleTime_s=10.0;
turnsLoading=3;
turnsUnloading=3;
numTrials=4;
initialZero_s=5.0;
zeroHold_s=10.0;
alternateDirection=true; % CW/CCW alternate BETWEEN trials, at phi=0.
% phase: -1 invalid,0 baseline,4 zero hold,6 loading helix,
% 7 unloading helix,9 complete. Rotation does not pause at maximum phi.
assert(peakPhi_rad>=0 && peakPhi_rad<=pi && circleTime_s>0);
assert(turnsLoading>=1 && turnsLoading==floor(turnsLoading));
assert(turnsUnloading>=1 && turnsUnloading==floor(turnsUnloading));
assert(numTrials>=1 && numTrials==floor(numTrials) && initialZero_s>=0 && zeroHold_s>=0);
loading_s=turnsLoading*circleTime_s;
unloading_s=turnsUnloading*circleTime_s;
moving_s=loading_s+unloading_s;
trialDuration_s=moving_s+zeroHold_s;
turnsTotal=turnsLoading+turnsUnloading;
testEnd_s=initialZero_s+numTrials*trialDuration_s;
theta=0;phi=0;trialId=0;phaseId=0;cycleId=0;direction=0;done=false;
if ~isfinite(t) || t<0,phaseId=-1;return;end
if t<initialZero_s,return;end
if t>=testEnd_s
    turnsAccumulated=turnsTotal*numTrials;
    if alternateDirection,turnsAccumulated=turnsTotal*mod(numTrials,2);end
    theta=2*pi*turnsAccumulated;phaseId=9;done=true;return;
end
index=floor((t-initialZero_s)/trialDuration_s);
trialId=index+1;local=t-initialZero_s-index*trialDuration_s;
rotationDirection=1;turnsBefore=turnsTotal*index;
if alternateDirection
    rotationDirection=1-2*mod(index,2);
    turnsBefore=turnsTotal*mod(index,2);
end
thetaStart=2*pi*turnsBefore;
theta=thetaStart+rotationDirection*2*pi*min(local,moving_s)/circleTime_s;
if local<loading_s
    phaseId=6;direction=rotationDirection;
    cycleId=floor(local/circleTime_s)+1;
    u=local/loading_s;phi=peakPhi_rad*u*u*(3-2*u);
elseif local<moving_s
    phaseId=7;direction=rotationDirection;
    cycleId=floor(local/circleTime_s)+1;
    u=(local-loading_s)/unloading_s;
    phi=peakPhi_rad*(1-u*u*(3-2*u));
else
    phaseId=4;
end
end
