function [theta,phi,testId,phaseId] = pacStaticSweepReferenceBlock(t)
%#codegen
% Paste into reference MATLAB Function block. All settings are internal.
% 10 trials: three amplitudes in three directions, then repeat theta=0/.2.
% Stop Time 535 s. Start only if controller status=3 before 25 s.
start_s=25.0;
rise_s=5.0; hold_s=30.0; fall_s=5.0; zeroHold_s=10.0;
thetaList=[0 0 0 pi/2 pi/2 pi/2 pi pi pi 0];
phiList=[.1 .2 .3 .1 .2 .3 .1 .2 .3 .2];
theta=0.0;phi=0.0;testId=0.0;phaseId=0.0;
if ~isfinite(t) || t<start_s, return; end
trial_s=rise_s+hold_s+fall_s+zeroHold_s;
idx=floor((t-start_s)/trial_s)+1;
if idx>numel(phiList),return;end
testId=double(idx);theta=thetaList(idx);
local=t-start_s-(idx-1)*trial_s;
if local<rise_s
    phaseId=1.0;u=local/rise_s;phi=phiList(idx)*u*u*(3-2*u);
elseif local<rise_s+hold_s
    phaseId=2.0;phi=phiList(idx);
elseif local<rise_s+hold_s+fall_s
    phaseId=3.0;u=(local-rise_s-hold_s)/fall_s;
    phi=phiList(idx)*(1-u*u*(3-2*u));
else
    phaseId=4.0;
end
end
