function [pressure_bar, tauDesired_N, tauAvailable_N, limited, modelValid] = inverseModel(l, orientationValid, stiff)
%#codegen
% measured = [l2;l3;dl2;dl3;ddl2;ddl3] in SI units.
% Nonnegative pressure mapping uses the smallest common P1 pressure.

pressure_bar = zeros(3,1);
tauDesired_N = zeros(2,1);
tauAvailable_N = zeros(2,1);
limited = false;
modelValid = false;
if ~orientationValid || ~all(isfinite(l))
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
% params.K = stiff * [2 1; 1 2];
params.D = 40*eye(2);
params.mu = 2000;
params.lKbounds = [-0.035;0.035;1e6];
[M,C,G] = armS_single_core(q,dq,params);
Keff = params.K;
for i = 1:2
    Keff(i,i) = Keff(i,i) + 0.5*params.lKbounds(3)*(2 + ...
        tanh(params.mu*(q(i)-params.lKbounds(2))) - ...
        tanh(params.mu*(q(i)-params.lKbounds(1))));
end
tauDesired_N = M*ddq + (C+params.D)*dq + G + Keff*q;
if ~all(isfinite(tauDesired_N))
    tauDesired_N = zeros(2,1);
    return;
end
forcePerBar = pi*(params.r/2)^2*1e5;
forceBar = tauDesired_N/forcePerBar;
p1 = max([0;-forceBar(1);-forceBar(2)]);
effective = max([p1;p1+forceBar(1);p1+forceBar(2)],0);
deadzone_bar = 0.0;
maxCommand_bar = 3.0;
maxEffective_bar = maxCommand_bar-deadzone_bar;
peak = max(effective);
if peak > maxEffective_bar
    effective = effective*(maxEffective_bar/peak);
    limited = true;
end
for i = 1:3
    if effective(i) > 1e-8
        pressure_bar(i) = effective(i)+deadzone_bar;
    end
end
tauAvailable_N = forcePerBar*[-1 1 0;-1 0 1]*effective;
modelValid = true;
end
