function [pressure_bar, tauDesired_N, tauAvailable_N, limited, modelValid, ...
    MassForce_N, CoriolisForce_N, GnStiff_N] = ...
    commonPressureFeedbackInverseBlock(measured, desiredAcceleration, orientationValid, t)
%#codegen
% Paste into the existing inverseModel MATLAB Function block in a COPY of
% the feedback-linearization model. Add Clock time t as the fourth input.
% Outputs and their order match the existing block.
%
% Test sequence: zero 0-10 s, fixed trial-4 trim 10-50 s, bounded
% feedback 50-220 s, zero after 220 s. Set model Stop Time to 240 s.
% The block rate-limits each pressure command to +/-0.8 bar/s.

pressure_bar = zeros(3,1);
tauDesired_N = zeros(2,1);
tauAvailable_N = zeros(2,1);
MassForce_N = zeros(2,1);
CoriolisForce_N = zeros(2,1);
GnStiff_N = zeros(2,1);
limited = false;
modelValid = false;

persistent tauAnchor_N anchorReady lastPressure_bar lastTime_s
if isempty(anchorReady)
    tauAnchor_N = zeros(2,1);
    anchorReady = false;
    lastPressure_bar = zeros(3,1);
    lastTime_s = -1.0;
end

if ~isfinite(t) || t < 0.0 || ...
        ~orientationValid || ~all(isfinite(measured(1:4))) || ...
        ~all(isfinite(desiredAcceleration))
    % Keep the model's fail-to-zero behavior on invalid NDI/model input.
    lastPressure_bar = zeros(3,1);
    if isfinite(t), lastTime_s = t; end
    return;
end

q = measured(1:2);
dq = measured(3:4);
params.N = 1;
params.L = 0.17411;
params.r = 0.013; % Radial PMA offset, not pressure-area radius.
params.mi = 0.1;
params.g = [0;0;-9.81];
params.K = 650*[2 1;1 2];
params.D = 40*eye(2);
params.mu = 2000;
params.lKbounds = [-0.035;0.035;1e6];
[M,C,G] = armS_single_core(q,dq,params);
bound = armS_single_boundForce(q,params.mu,params.lKbounds);

MassForce_N = M*desiredAcceleration;
CoriolisForce_N = (C+params.D)*dq;
GnStiff_N = G+params.K*q+bound;
tauDesired_N = MassForce_N+CoriolisForce_N+GnStiff_N;
if ~all(isfinite(tauDesired_N)) || ~all(isfinite(GnStiff_N))
    tauDesired_N = zeros(2,1);
    lastPressure_bar = zeros(3,1);
    lastTime_s = t;
    return;
end
modelValid = true;

common_bar = 1.6;
s0_bar = -0.0983981693363845;
a0_bar = -0.0819816714043049;
maxDeltaS_bar = 0.04;
maxDeltaA_bar = 0.04;
pressureRate_barPerS = 0.8;
forcePerBar_N = pi*(0.013/2)^2*1e5;

target_bar = zeros(3,1);
if t >= 10.0 && t < 220.0
    ds_bar = 0.0;
    da_bar = 0.0;
    if t >= 50.0
        if ~anchorReady
            % Anchor only the nominal static model force at the settled
            % trim. Do not cancel the PD correction to the desired pose.
            tauAnchor_N = GnStiff_N;
            anchorReady = true;
        end
        deltaTau_N = tauDesired_N-tauAnchor_N;
        dsWanted_bar = (deltaTau_N(1)+deltaTau_N(2))/(2*forcePerBar_N);
        daWanted_bar = (deltaTau_N(2)-deltaTau_N(1))/(2*forcePerBar_N);
        ds_bar = min(max(dsWanted_bar,-maxDeltaS_bar),maxDeltaS_bar);
        da_bar = min(max(daWanted_bar,-maxDeltaA_bar),maxDeltaA_bar);
        limited = abs(ds_bar-dsWanted_bar)>1e-12 || ...
            abs(da_bar-daWanted_bar)>1e-12;
        % Gradually engage feedback after the static trim hold.
        u = min(max((t-50.0)/5.0,0.0),1.0);
        blend = u*u*(3.0-2.0*u);
        ds_bar = blend*ds_bar;
        da_bar = blend*da_bar;
    end
    s_bar = s0_bar+ds_bar;
    a_bar = a0_bar+da_bar;
    target_bar = [common_bar-(2.0/3.0)*s_bar; ...
                  common_bar+(1.0/3.0)*s_bar-a_bar; ...
                  common_bar+(1.0/3.0)*s_bar+a_bar];
    bounded = min(max(target_bar,1.2),3.0);
    if any(abs(bounded-target_bar)>1e-12), limited = true; end
    target_bar = bounded;
end

% Discrete state implements a per-channel slew limit using Clock time.
% The initial trim takes about 2.1 s to reach its target.
if lastTime_s < 0.0 || t < lastTime_s
    lastPressure_bar = zeros(3,1);
    lastTime_s = t;
end
dt = t-lastTime_s;
if dt > 0.0
    step_bar = pressureRate_barPerS*dt;
    pressure_bar = min(max(target_bar, ...
        lastPressure_bar-step_bar),lastPressure_bar+step_bar);
    lastPressure_bar = pressure_bar;
    lastTime_s = t;
else
    pressure_bar = lastPressure_bar;
end
if any(abs(pressure_bar-target_bar)>1e-9), limited = true; end

% All three active channels have the same nominal 0.8 bar deadzone, so
% their pressure differences give the ideal generalized force directly.
tauAvailable_N = forcePerBar_N * ...
    [pressure_bar(2)-pressure_bar(1); ...
     pressure_bar(3)-pressure_bar(1)];
end
