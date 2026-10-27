function [theta,phi,trialId,phaseId,cycleId,direction,done,testEnd_s] = rlCircleSweepReferenceBlock(t)
%#codegen
% Requested ANGLES, not pressures or measured arm motion.
% Five phi levels, five circles each, no pauses between revolutions.
phiLevels_rad=[0.2 0.8 1.0 1.2 1.5];
circleTime_s=20.0;
circlesPerTrial=5;
initialZero_s=5.0;
phiRise_s=5.0;
targetWait_s=3.0;
phiFall_s=5.0;
zeroHold_s=10.0;
rotationDirection=1; % +1 or -1; keep unwrapped theta for downstream rates.
% phase: -1 invalid,0 baseline,1 phi rise,2 target wait,3 phi fall,
% 4 zero hold,5 circles,9 complete.
assert(all(phiLevels_rad>=0 & phiLevels_rad<=pi));
assert(circleTime_s>0 && circlesPerTrial>=1 && circlesPerTrial==floor(circlesPerTrial));
assert(phiRise_s>0 && phiFall_s>0 && targetWait_s>=0 && zeroHold_s>=0 && initialZero_s>=0);
assert(abs(rotationDirection)==1);
circleDuration_s=circlesPerTrial*circleTime_s;
trialDuration_s=phiRise_s+targetWait_s+circleDuration_s+phiFall_s+zeroHold_s;
nTrials=numel(phiLevels_rad);
testEnd_s=initialZero_s+nTrials*trialDuration_s;
theta=0;phi=0;trialId=0;phaseId=0;cycleId=0;direction=0;done=false;
if ~isfinite(t) || t<0,phaseId=-1;return;end
if t<initialZero_s,return;end
if t>=testEnd_s
    theta=rotationDirection*2*pi*circlesPerTrial*nTrials;
    phaseId=9;done=true;return;
end
index=floor((t-initialZero_s)/trialDuration_s);
trialId=index+1;
local=t-initialZero_s-index*trialDuration_s;
targetPhi=phiLevels_rad(trialId);
thetaStart=rotationDirection*2*pi*circlesPerTrial*index;
theta=thetaStart;
circleStart=phiRise_s+targetWait_s;
circleEnd=circleStart+circleDuration_s;
if local<phiRise_s
    phaseId=1;u=local/phiRise_s;phi=targetPhi*u*u*(3-2*u);
elseif local<circleStart
    phaseId=2;phi=targetPhi;
elseif local<circleEnd
    phaseId=5;direction=rotationDirection;
    elapsed=local-circleStart;
    cycleId=floor(elapsed/circleTime_s)+1;
    theta=thetaStart+rotationDirection*2*pi*elapsed/circleTime_s;
    phi=targetPhi;
elseif local<circleEnd+phiFall_s
    phaseId=3;theta=thetaStart+rotationDirection*2*pi*circlesPerTrial;
    u=(local-circleEnd)/phiFall_s;phi=targetPhi*(1-u*u*(3-2*u));
else
    phaseId=4;theta=thetaStart+rotationDirection*2*pi*circlesPerTrial;
end
end
