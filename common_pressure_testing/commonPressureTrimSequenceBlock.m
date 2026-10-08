function [P1_bar, P2_bar, P3_bar, phaseID] = commonPressureTrimSequenceBlock(t)
%#codegen
% Paste this function into a MATLAB Function block; feed t from Clock.
% Connect P1_bar/P2_bar/P3_bar to the existing pressure-command Mux.
% Log phaseID. Set model Stop Time to 387 s.
%
% Open-loop neutral-trim calibration at mean commands 1.2, 1.6, 2.0 bar.
% Each level: baseline, symmetric +/-, antisymmetric +/-, returning to
% baseline after each pulse. The pulse size is 0.2 bar for symmetric
% differential and 0.1 bar for antisymmetric differential.
%
% p = [c-2*s/3; c+s/3-a; c+s/3+a], so mean(p) = c.
% All active-level commands are 1.067 to 2.133 bar. Initial and final
% ramps necessarily pass through the nominal 0.8 bar deadzone.
% No automatic pose-based adjustment or pressure clipping is performed.

P1_bar = 0.0;
P2_bar = 0.0;
P3_bar = 0.0;
phaseID = 0.0;
if ~isfinite(t) || t < 0.0
    return;
end

initialZero = 10.0;
groupSeconds = 119.0;
groupsEnd = initialZero + 3.0*groupSeconds;  % 367 s
returnSeconds = 5.0;

if t < initialZero
    return;
elseif t < groupsEnd
    level = floor((t - initialZero)/groupSeconds) + 1.0;
    localTime = t - initialZero - (level - 1.0)*groupSeconds;
    targetCommon = 0.8 + 0.4*level;
    previousCommon = 0.0;
    if level > 1.0
        previousCommon = targetCommon - 0.4;
    end

    s = 0.0;
    a = 0.0;
    if localTime < 5.0
        z = localTime/5.0;
        blend = z*z*(3.0 - 2.0*z);
        c = previousCommon + blend*(targetCommon - previousCommon);
        phaseID = 100.0*level + 1.0;  % baseline ramp
    elseif localTime < 15.0
        c = targetCommon;
        phaseID = 100.0*level + 2.0;  % initial baseline hold
    else
        c = targetCommon;
        pulseIndex = floor((localTime - 15.0)/26.0);  % 0,1,2,3
        pulseTime = localTime - 15.0 - 26.0*pulseIndex;
        phaseBase = 100.0*level + 10.0*(pulseIndex + 1.0);

        if pulseTime < 4.0
            z = pulseTime/4.0;
            weight = z*z*(3.0 - 2.0*z);
            phaseID = phaseBase + 1.0;  % ramp to perturbation
        elseif pulseTime < 14.0
            weight = 1.0;
            phaseID = phaseBase + 2.0;  % perturbation hold
        elseif pulseTime < 18.0
            z = (pulseTime - 14.0)/4.0;
            weight = 1.0 - z*z*(3.0 - 2.0*z);
            phaseID = phaseBase + 3.0;  % return ramp
        else
            weight = 0.0;
            phaseID = phaseBase + 4.0;  % return baseline hold
        end

        if pulseIndex == 0.0
            s =  0.2*weight;  % P2 and P3 rise relative to P1
        elseif pulseIndex == 1.0
            s = -0.2*weight;
        elseif pulseIndex == 2.0
            a =  0.1*weight;  % P3 rises relative to P2
        else
            a = -0.1*weight;
        end
    end
elseif t < groupsEnd + returnSeconds
    z = (t - groupsEnd)/returnSeconds;
    blend = z*z*(3.0 - 2.0*z);
    c = 2.0*(1.0 - blend);
    s = 0.0;
    a = 0.0;
    phaseID = 900.0;  % final ramp to zero
else
    phaseID = 901.0;  % zero hold after 372 s, including after Stop Time
    return;
end

P1_bar = c - (2.0/3.0)*s;
P2_bar = c + (1.0/3.0)*s - a;
P3_bar = c + (1.0/3.0)*s + a;
end
