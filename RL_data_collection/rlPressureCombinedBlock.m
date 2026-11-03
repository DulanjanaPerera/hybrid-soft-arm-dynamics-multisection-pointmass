function [pressureRequest_bar,trialId,phaseId,channelId,cycleId,done,testEnd_s] = rlPressureCombinedBlock(t)
%#codegen
% OPEN-LOOP common/differential pressure excitation, not angle tracking.
% trialId: 1 common pulse, 2 symmetric s, 3 antisymmetric a, 4 rotating s/a.
% channelId=0 means combined channels. No deadzone compensation.
commonPeak_bar=2.0;
bias_bar=1.5;
differentialAmplitude_bar=0.5;
omega_radPerS=0.5;
cyclesPerTrial=2;
initialZero_s=5.0;
preloadRise_s=5.0;
preloadFall_s=5.0;
relax_s=10.0;
% phase: -1 invalid, 0 initial zero, 1 preload rise, 3 preload fall,
% 4 zero relaxation, 5 periodic excitation, 9 complete.
assert(commonPeak_bar>=0 && commonPeak_bar<=3 && differentialAmplitude_bar>=0);
% Conservative bound for the rotating differential combination.
assert(bias_bar-4*differentialAmplitude_bar/3>=0 && bias_bar+4*differentialAmplitude_bar/3<=3);
assert(omega_radPerS>0 && cyclesPerTrial>=1 && cyclesPerTrial==floor(cyclesPerTrial));
assert(initialZero_s>=0 && preloadRise_s>0 && preloadFall_s>0 && relax_s>0);
period_s=2*pi/omega_radPerS;
excitation_s=cyclesPerTrial*period_s;
trialDuration_s=preloadRise_s+excitation_s+preloadFall_s+relax_s;
testEnd_s=initialZero_s+4*trialDuration_s;
pressureRequest_bar=zeros(3,1);trialId=0;phaseId=0;channelId=0;cycleId=0;done=false;
if ~isfinite(t) || t<0,phaseId=-1;return;end
if t<initialZero_s,return;end
if t>=testEnd_s,phaseId=9;done=true;return;end
index=floor((t-initialZero_s)/trialDuration_s);trialId=index+1;
local=t-initialZero_s-index*trialDuration_s;
% Each mode has matching endpoints; preload/fall remove inter-mode jumps.
endpoint=zeros(3,1);
if trialId==2 || trialId==3,endpoint(:)=bias_bar;end
if trialId==4
    endpoint=[bias_bar-2*differentialAmplitude_bar/3; ...
        bias_bar+differentialAmplitude_bar/3;bias_bar+differentialAmplitude_bar/3];
end
if local<preloadRise_s
    phaseId=1;u=local/preloadRise_s;pressureRequest_bar=endpoint*u*u*(3-2*u);
elseif local<preloadRise_s+excitation_s
    phaseId=5;tau=local-preloadRise_s;angle=omega_radPerS*tau;
    cycleId=floor(tau/period_s)+1;
    if trialId==1
        pressureRequest_bar(:)=commonPeak_bar*.5*(1-cos(angle));
    else
        s=0;a=0;
        if trialId==2,s=differentialAmplitude_bar*sin(angle);end
        if trialId==3,a=differentialAmplitude_bar*sin(angle);end
        if trialId==4
            s=differentialAmplitude_bar*cos(angle);
            a=differentialAmplitude_bar*sin(angle);
        end
        pressureRequest_bar=[bias_bar-2*s/3;bias_bar+s/3-a;bias_bar+s/3+a];
    end
elseif local<preloadRise_s+excitation_s+preloadFall_s
    phaseId=3;u=(local-preloadRise_s-excitation_s)/preloadFall_s;
    pressureRequest_bar=endpoint*(1-u*u*(3-2*u));
else
    phaseId=4;
end
end
