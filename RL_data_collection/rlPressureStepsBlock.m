function [pressureRequest_bar,trialId,phaseId,channelId,cycleId,done,testEnd_s] = rlPressureStepsBlock(t)
%#codegen
% OPEN-LOOP REQUEST ONLY. Connect through the existing Kill/pressure limits.
% channelSelection: 1=P1, 2=P2, 3=P3, 0=all three tested separately.
channelSelection=0;
peakLevels_bar=[1.0 2.0 3.0];
repeats=2;
initialZero_s=5.0;
rise_s=3.0;
hold_s=10.0;
fall_s=3.0;
relax_s=10.0;
% phase: -1 invalid time, 0 initial zero, 1 rise, 2 hold, 3 fall,
% 4 zero relaxation, 9 complete. No deadzone subtraction.
assert(channelSelection>=0 && channelSelection<=3 && channelSelection==floor(channelSelection));
assert(repeats>=1 && repeats==floor(repeats));
assert(all(peakLevels_bar>=0 & peakLevels_bar<=3));
assert(initialZero_s>=0 && rise_s>0 && hold_s>0 && fall_s>0 && relax_s>0);
nChannels=1;if channelSelection==0,nChannels=3;end
nLevels=numel(peakLevels_bar);
trialDuration_s=rise_s+hold_s+fall_s+relax_s;
nTrials=nChannels*repeats*nLevels;
testEnd_s=initialZero_s+nTrials*trialDuration_s;
pressureRequest_bar=zeros(3,1);trialId=0;phaseId=0;channelId=0;cycleId=0;done=false;
if ~isfinite(t) || t<0,phaseId=-1;return;end
if t<initialZero_s,return;end
if t>=testEnd_s,phaseId=9;done=true;return;end
index=floor((t-initialZero_s)/trialDuration_s);
trialId=index+1;
channelId=channelSelection;
if channelSelection==0,channelId=floor(index/(repeats*nLevels))+1;end
cycleId=mod(floor(index/nLevels),repeats)+1;
level=mod(index,nLevels)+1;
local=t-initialZero_s-index*trialDuration_s;
peak=peakLevels_bar(level);
if local<rise_s
    phaseId=1;u=local/rise_s;p=peak*u*u*(3-2*u);
elseif local<rise_s+hold_s
    phaseId=2;p=peak;
elseif local<rise_s+hold_s+fall_s
    phaseId=3;u=(local-rise_s-hold_s)/fall_s;p=peak*(1-u*u*(3-2*u));
else
    phaseId=4;p=0;
end
pressureRequest_bar(channelId)=p;
end
