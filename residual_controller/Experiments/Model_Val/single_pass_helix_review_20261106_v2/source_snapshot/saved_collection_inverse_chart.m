function [pressure_bar, tauDesired_N, tauAvailable_N, limited, modelValid] = ...
    rlInverseModelBlock(l, orientationValid, stiff)
%#codegen
% Paste-ready replacement for the collection inverseModel MATLAB Function.
% Keeps the existing 3 inputs / 5 outputs, invalid hold and allocation policy.
% Requires the corrected armS_single_boundForce.m on the MATLAB path.
% Bound value remains +/-0.035 m; no SLX edited or hardware run.
% l = [l2;l3;dl2;dl3;ddl2;ddl3], SI units.
% Invalid input/model evaluation holds the last successful outputs.
% modelValid=false identifies held outputs.
% Persistent values reset when the block is initialized for a new run.
% stiff remains unused, matching the original implementation.

persistent lastPressure lastTauDesired lastTauAvailable lastLimited

if isempty(lastPressure)
    lastPressure = zeros(3,1);
    lastTauDesired = zeros(2,1);
    lastTauAvailable = zeros(2,1);
    lastLimited = false;
end

% Default to the last successful outputs.
pressure_bar = lastPressure;
tauDesired_N = lastTauDesired;
tauAvailable_N = lastTauAvailable;
limited = lastLimited;
modelValid = false;

if ~orientationValid || ~all(isfinite(l(:)))
    return;
end

q = l(1:2);
dq = l(3:4);
ddq = l(5:6);

params.N = 1;
params.L = 0.17411;
params.r = 0.013;
params.mi = 0.1;
params.g = [0;0;-9.81];
params.K = 650*[2 1;1 2];
% params.K = stiff * [2 1;1 2];
params.D = 40*eye(2);
params.mu = 2000;
params.lKbounds = [-0.035;0.035;1e6];

[M,C,G] = armS_single_core(q,dq,params);

if ~all(isfinite(M(:))) || ~all(isfinite(C(:))) || ...
        ~all(isfinite(G(:)))
    return;
end

% Equal soft bounds on all three actuator lengths, including l1=-l2-l3.
bound = armS_single_boundForce(q,params.mu,params.lKbounds);
tauCandidate = M*ddq + (C+params.D)*dq + G + params.K*q + bound;
if ~all(isfinite(tauCandidate(:)))
    return;
end

forcePerBar = pi*(params.r/2)^2*1e5;
forceBar = tauCandidate/forcePerBar;
p1 = max([0;-forceBar(1);-forceBar(2)]);
effective = max([p1;p1+forceBar(1);p1+forceBar(2)],0);

deadzone_bar = 0.0;
maxCommand_bar = 3.0;
maxEffective_bar = maxCommand_bar-deadzone_bar;

limitedCandidate = false;
peak = max(effective);
if peak > maxEffective_bar
    effective = effective*(maxEffective_bar/peak);
    limitedCandidate = true;
end

pressureCandidate = zeros(3,1);
for i = 1:3
    if effective(i) > 1e-8
        pressureCandidate(i) = effective(i)+deadzone_bar;
    end
end

availableCandidate = ...
    forcePerBar*[-1 1 0;-1 0 1]*effective;

if ~all(isfinite(pressureCandidate(:))) || ...
        ~all(isfinite(availableCandidate(:)))
    return;
end

% Commit only a successful calculation.
lastPressure = pressureCandidate;
lastTauDesired = tauCandidate;
lastTauAvailable = availableCandidate;
lastLimited = limitedCandidate;

pressure_bar = lastPressure;
tauDesired_N = lastTauDesired;
tauAvailable_N = lastTauAvailable;
limited = lastLimited;
modelValid = true;
end