function [P1_bar, P2_bar, P3_bar, phaseID] = commonPressureSequenceBlock(t)
%#codegen
% Paste this entire function into one Simulink MATLAB Function block.
% Input t is simulation time in seconds from a Clock block.
% Connect P1_bar, P2_bar, P3_bar to the existing P1/P2/P3 Mux inputs.
% Log phaseID alongside the measured and commanded signals.
%
% One 300 s open-loop test:
%   0 bar reference -> common 1.2, 1.6, 2.0 bar -> 0 bar.
%   At each nonzero common level, try +0.2 and -0.2 bar differential
%   pressure, returning to equal pressure between bends.
%   Every active-level command is between 1.067 and 2.067 bar.
%
% Columns: target common pressure (bar), target differential (bar),
%          ramp duration (s), hold duration (s).
% A positive differential makes P2=P3>P1 (nominal theta=0 direction).
schedule = [ ...
    0.0,  0.0, 0, 10; ...  % phase 1: zero reference
    1.2,  0.0, 5, 15; ...  % phases 2-6: common 1.2 bar
    1.2,  0.2, 5, 15; ...
    1.2,  0.0, 5, 10; ...
    1.2, -0.2, 5, 15; ...
    1.2,  0.0, 5, 10; ...
    1.6,  0.0, 5, 15; ...  % phases 7-11: common 1.6 bar
    1.6,  0.2, 5, 15; ...
    1.6,  0.0, 5, 10; ...
    1.6, -0.2, 5, 15; ...
    1.6,  0.0, 5, 10; ...
    2.0,  0.0, 5, 15; ...  % phases 12-16: common 2.0 bar
    2.0,  0.2, 5, 15; ...
    2.0,  0.0, 5, 10; ...
    2.0, -0.2, 5, 15; ...
    2.0,  0.0, 5, 10; ...
    0.0,  0.0, 5, 15];     % phase 17: return to zero

P1_bar = 0.0;
P2_bar = 0.0;
P3_bar = 0.0;
phaseID = 0.0;
if ~isfinite(t) || t < 0.0
    return;
end

elapsed = 0.0;
previousCommon = 0.0;
previousDifferential = 0.0;
common = 0.0;
differential = 0.0;
for k = 1:size(schedule, 1)
    targetCommon = schedule(k, 1);
    targetDifferential = schedule(k, 2);
    rampSeconds = schedule(k, 3);
    holdSeconds = schedule(k, 4);

    if rampSeconds > 0.0 && t < elapsed + rampSeconds
        fraction = (t - elapsed)/rampSeconds;
        blend = fraction*fraction*(3.0 - 2.0*fraction);
        common = previousCommon + blend*(targetCommon - previousCommon);
        differential = previousDifferential + ...
            blend*(targetDifferential - previousDifferential);
        phaseID = double(k);
        break;
    end
    elapsed = elapsed + rampSeconds;

    if t < elapsed + holdSeconds
        common = targetCommon;
        differential = targetDifferential;
        phaseID = double(k);
        break;
    end
    elapsed = elapsed + holdSeconds;
    previousCommon = targetCommon;
    previousDifferential = targetDifferential;
end

% This keeps the mean pressure equal to common and gives
% P2-P1 = P3-P1 = differential throughout each active test.
P1_bar = common - (2.0/3.0)*differential;
P2_bar = common + (1.0/3.0)*differential;
P3_bar = P2_bar;

% Fixed schedule never reaches these limits; keep outputs bounded if edited.
P1_bar = min(max(P1_bar, 0.0), 3.0);
P2_bar = min(max(P2_bar, 0.0), 3.0);
P3_bar = min(max(P3_bar, 0.0), 3.0);
end
