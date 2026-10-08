function [P1_bar, P2_bar, P3_bar, phaseID] = commonPressureLocalCalibrationBlock(t)
%#codegen
% Paste this entire function into one Simulink MATLAB Function block.
% Input t: simulation time (s) from a Clock block.
% Outputs: three scalar pressure targets (bar) and a phase label.
% Positive phaseID = hold; negative phaseID = transition to that hold.
% Stop time: 395 s. The output stays at zero after 395 s.
%
% This is an open-loop local test around the best manual trial-4 trim.
% It applies one differential setting at a time and returns to trim
% between tests. Keep the model's downstream Rate Limiter and Kill path.

common_bar = 1.6;
s0_bar = -0.0983981693363845;
a0_bar = -0.0819816714043049;
deltaS_bar = 0.04;
deltaA_bar = 0.04;

% Columns: pressure enabled, s target, a target, ramp seconds, hold seconds.
% phase 1: zero; 2: trim; 3/5: s+/s-; 7/9: a+/a-;
% 4/6/8/10: return to trim; 11: zero.
schedule = [ ...
    0, 0,                 0,                 0, 15; ...
    1, s0_bar,            a0_bar,            5, 40; ...
    1, s0_bar+deltaS_bar, a0_bar,            5, 35; ...
    1, s0_bar,            a0_bar,            5, 30; ...
    1, s0_bar-deltaS_bar, a0_bar,            5, 35; ...
    1, s0_bar,            a0_bar,            5, 30; ...
    1, s0_bar,            a0_bar+deltaA_bar, 5, 35; ...
    1, s0_bar,            a0_bar,            5, 30; ...
    1, s0_bar,            a0_bar-deltaA_bar, 5, 35; ...
    1, s0_bar,            a0_bar,            5, 40; ...
    0, 0,                 0,                 5, 20];

P1_bar = 0.0;
P2_bar = 0.0;
P3_bar = 0.0;
phaseID = 0.0;
if ~isfinite(t) || t < 0.0
    return;
end

elapsed = 0.0;
previous = [0.0, 0.0, 0.0];
command = previous;
phaseID = 11.0; % All times after the schedule remain at zero.
for k = 1:size(schedule,1)
    if schedule(k,1) > 0.5
        s = schedule(k,2);
        a = schedule(k,3);
        target = [common_bar-(2.0/3.0)*s, ...
                  common_bar+(1.0/3.0)*s-a, ...
                  common_bar+(1.0/3.0)*s+a];
    else
        target = [0.0, 0.0, 0.0];
    end

    ramp = schedule(k,4);
    holdTime = schedule(k,5);
    if ramp > 0.0 && t < elapsed+ramp
        u = (t-elapsed)/ramp;
        blend = u*u*(3.0-2.0*u); % Smooth zero-slope transition.
        command = previous+blend*(target-previous);
        phaseID = -double(k);
        break;
    end
    elapsed = elapsed+ramp;
    if t < elapsed+holdTime
        command = target;
        phaseID = double(k);
        break;
    end
    elapsed = elapsed+holdTime;
    previous = target;
end

% Guard against accidental parameter changes beyond the model range.
P1_bar = min(max(command(1),0.0),3.0);
P2_bar = min(max(command(2),0.0),3.0);
P3_bar = min(max(command(3),0.0),3.0);
end
