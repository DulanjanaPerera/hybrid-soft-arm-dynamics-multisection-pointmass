function [P1_bar, P2_bar, P3_bar, phaseID] = commonPressureTrimValidationBlock(t)
%#codegen
% Paste into a MATLAB Function block. Feed t from a Clock block.
% Connect P1_bar/P2_bar/P3_bar to the existing pressure-command Mux.
% Log phaseID and set Stop Time to 199 s.
%
% Open-loop validation of the estimated mean-1.6-bar trim.
% Route A: zero -> 1.2 -> equal 1.6 -> two intermediate holds -> trim.
% Route B: zero -> equal 1.6 -> trim directly.
% The target [1.766, 1.766, 1.268] bar is an extrapolated estimate;
% this sequence measures it but does not adjust pressure from NDI.
%
% Each row: P1, P2, P3 (bar), smooth ramp time (s), hold time (s).
% Phase ID is 10*row+1 during a ramp and 10*row+2 during a hold.
schedule = [ ...
    0.000, 0.000, 0.000, 0, 10; ...  % 12: initial zero reference
    1.200, 1.200, 1.200, 5, 10; ...  % 22: route A, 1.2 baseline
    1.600, 1.600, 1.600, 5, 15; ...  % 32: route A, 1.6 baseline
    1.650, 1.650, 1.500, 4, 10; ...  % 42: intermediate hold 1
    1.700, 1.700, 1.400, 4, 10; ...  % 52: intermediate hold 2
    1.766, 1.766, 1.268, 4, 20; ...  % 62: candidate, route A
    1.600, 1.600, 1.600, 5, 10; ...  % 72: return to baseline
    0.000, 0.000, 0.000, 5, 10; ...  % 82: zero between routes
    1.600, 1.600, 1.600, 7, 15; ...  % 92: route B, 1.6 baseline
    1.766, 1.766, 1.268, 8, 20; ...  % 102: candidate, route B
    0.000, 0.000, 0.000, 7, 15];     % 112: final zero reference

P1_bar = 0.0;
P2_bar = 0.0;
P3_bar = 0.0;
phaseID = 0.0;
if ~isfinite(t) || t < 0.0
    return;
end

elapsed = 0.0;
previous = [0.0, 0.0, 0.0];
pressure = previous;
phaseID = 112.0; % Also labels zero output after Stop Time.
for row = 1:size(schedule, 1)
    target = schedule(row, 1:3);
    rampSeconds = schedule(row, 4);
    holdSeconds = schedule(row, 5);

    if rampSeconds > 0.0 && t < elapsed + rampSeconds
        z = (t - elapsed)/rampSeconds;
        blend = z*z*(3.0 - 2.0*z);
        pressure = previous + blend*(target - previous);
        phaseID = 10.0*row + 1.0;
        break;
    end
    elapsed = elapsed + rampSeconds;

    if t < elapsed + holdSeconds
        pressure = target;
        phaseID = 10.0*row + 2.0;
        break;
    end
    elapsed = elapsed + holdSeconds;
    previous = target;
    pressure = target;
end

P1_bar = pressure(1);
P2_bar = pressure(2);
P3_bar = pressure(3);
end
