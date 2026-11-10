function [pressureRequest_bar,trialId,phaseId,channelId,cycleId,done,testEnd_s] = rlPressureSineBlock(t)
%#codegen
% OPEN-LOOP raw pressure request: peak/2*(1-cos(omega*time)).
% channelSelection: 1=P1, 2=P2, 3=P3, 0=all three separately.
channelSelection=0;
peak_bar=3.0;
omega_radPerS=0.8;
cyclesPerTrial=3;
repeats=2;
initialZero_s=5.0;
relax_s=10.0;
% phase: -1 invalid time, 0 initial zero, 4 zero relaxation,
% 5 periodic excitation, 9 complete. Cycles join at zero with zero slope.
assert(channelSelection>=0 && channelSelection<=3 && channelSelection==floor(channelSelection));
assert(peak_bar>=0 && peak_bar<=3 && omega_radPerS>0);
assert(cyclesPerTrial>=1 && cyclesPerTrial==floor(cyclesPerTrial));
assert(repeats>=1 && repeats==floor(repeats) && initialZero_s>=0 && relax_s>0);
nChannels=1;if channelSelection==0,nChannels=3;end
period_s=2*pi/omega_radPerS;
excitation_s=cyclesPerTrial*period_s;
trialDuration_s=excitation_s+relax_s;
nTrials=nChannels*repeats;
testEnd_s=initialZero_s+nTrials*trialDuration_s;
pressureRequest_bar=zeros(3,1);trialId=0;phaseId=0;channelId=0;cycleId=0;done=false;
if ~isfinite(t) || t<0,phaseId=-1;return;end
if t<initialZero_s,return;end
if t>=testEnd_s,phaseId=9;done=true;return;end
index=floor((t-initialZero_s)/trialDuration_s);
trialId=index+1;channelId=channelSelection;
if channelSelection==0,channelId=floor(index/repeats)+1;end
local=t-initialZero_s-index*trialDuration_s;
if local<excitation_s
    phaseId=5;cycleId=floor(local/period_s)+1;
    pressureRequest_bar(channelId)=peak_bar*.5*(1-cos(omega_radPerS*local));
else
    phaseId=4;
end
end
