function result = compareKiAngleErrors(path1, path2)
%COMPAREKIANGLEERRORS Compare desired/measured angles in two saved runs.
%   path1 = 'C:\path\to\without_i_1rad.mat';
%   path2 = 'C:\path\to\withi_1rad.mat';
%   result = compareKiAngleErrors(path1, path2);
%
% Opens three figures: theta, phi, and their errors. No files are written.
% Theta traces and errors are wrapped to [-pi, pi]. Invalid NDI samples are
% omitted; measured theta is also omitted when its direction is unobservable.
% The plots use only the time interval covered by both recordings.

narginchk(2, 2);
run1 = readRun(path1);
run2 = readRun(path2);

tStart = max(run1.Time_s(1), run2.Time_s(1));
tEnd = min(run1.Time_s(end), run2.Time_s(end));
assert(tEnd > tStart, 'The recordings have no overlapping time interval.');
in1 = run1.Time_s >= tStart & run1.Time_s <= tEnd;
in2 = run2.Time_s >= tStart & run2.Time_s <= tEnd;

tCheck = run1.Time_s(in1);
thetaReferenceDifference = max(abs( ...
    interp1(run1.ReferenceThetaTime_s, run1.ReferenceTheta_rad, tCheck) - ...
    interp1(run2.ReferenceThetaTime_s, run2.ReferenceTheta_rad, tCheck)));
phiReferenceDifference = max(abs( ...
    interp1(run1.ReferencePhiTime_s, run1.ReferencePhi_rad, tCheck) - ...
    interp1(run2.ReferencePhiTime_s, run2.ReferencePhi_rad, tCheck)));
sameReference = thetaReferenceDifference < 1e-6 && ...
    phiReferenceDifference < 1e-6;
if ~sameReference
    warning('compareKiAngleErrors:DifferentReferences', ...
        'The desired trajectories differ. Both references will be shown.');
end

fprintf('Path 1: %s\n', run1.Path);
fprintf('  Kp %s; Kd %s; Ki %s\n', ...
    gainText(run1.Kp), gainText(run1.Kd), gainText(run1.Ki));
fprintf('Path 2: %s\n', run2.Path);
fprintf('  Kp %s; Kd %s; Ki %s\n', ...
    gainText(run2.Kp), gainText(run2.Kd), gainText(run2.Ki));
fprintf('Common comparison interval: %.2f to %.2f s\n', tStart, tEnd);
if ~isequal(unique(run1.Kp, 'rows'), unique(run2.Kp, 'rows')) || ...
        ~isequal(unique(run1.Kd, 'rows'), unique(run2.Kd, 'rows'))
    warning('compareKiAngleErrors:OtherGainsDiffer', ...
        ['Kp or Kd also differs between runs. Tracking differences ' ...
         'cannot be attributed to Ki alone.']);
end

blue = [0.000 0.447 0.741];
orange = [0.850 0.325 0.098];
referenceColor = [0.15 0.15 0.15];

figTheta = figure('Name', 'Ki comparison - theta', 'Color', 'w');
hold on;
plot(run1.Time_s(in1), breakWrap(wrapAngle(run1.DesiredTheta_rad(in1))), ...
    '--', 'Color', referenceColor, 'LineWidth', 1.4, ...
    'DisplayName', 'Desired theta');
if ~sameReference
    plot(run2.Time_s(in2), breakWrap(wrapAngle(run2.DesiredTheta_rad(in2))), ...
        ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
        'DisplayName', 'Desired theta, path 2');
end
plot(run1.Time_s(in1), breakWrap(wrapAngle(run1.MeasuredTheta_rad(in1))), ...
    'Color', blue, 'LineWidth', 1.2, 'DisplayName', run1.Label);
plot(run2.Time_s(in2), breakWrap(wrapAngle(run2.MeasuredTheta_rad(in2))), ...
    'Color', orange, 'LineWidth', 1.2, 'DisplayName', run2.Label);
hold off; grid on; xlim([tStart tEnd]); ylim([-pi pi]);
xlabel('Time (s)'); ylabel('Theta (rad)');
title('Desired and measured theta, wrapped to [-\pi, \pi]');
legend('Location', 'best');

figPhi = figure('Name', 'Ki comparison - phi', 'Color', 'w');
hold on;
plot(run1.Time_s(in1), run1.DesiredPhi_rad(in1), '--', ...
    'Color', referenceColor, 'LineWidth', 1.4, ...
    'DisplayName', 'Desired phi');
if ~sameReference
    plot(run2.Time_s(in2), run2.DesiredPhi_rad(in2), ':', ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
        'DisplayName', 'Desired phi, path 2');
end
plot(run1.Time_s(in1), run1.MeasuredPhi_rad(in1), ...
    'Color', blue, 'LineWidth', 1.2, 'DisplayName', run1.Label);
plot(run2.Time_s(in2), run2.MeasuredPhi_rad(in2), ...
    'Color', orange, 'LineWidth', 1.2, 'DisplayName', run2.Label);
hold off; grid on; xlim([tStart tEnd]);
xlabel('Time (s)'); ylabel('Phi (rad)');
title('Desired and measured phi');
legend('Location', 'best');

figError = figure('Name', 'Ki comparison - angle errors', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile; hold on;
plot(run1.Time_s(in1), breakWrap(run1.ThetaError_rad(in1)), ...
    'Color', blue, 'LineWidth', 1.2, 'DisplayName', run1.Label);
plot(run2.Time_s(in2), breakWrap(run2.ThetaError_rad(in2)), ...
    'Color', orange, 'LineWidth', 1.2, 'DisplayName', run2.Label);
yline(0, 'k:');
hold off; grid on; xlim([tStart tEnd]);
ylabel('Theta error (rad)');
title('Wrapped theta error: desired - measured');
legend('Location', 'best');
nexttile; hold on;
plot(run1.Time_s(in1), run1.PhiError_rad(in1), ...
    'Color', blue, 'LineWidth', 1.2, 'DisplayName', run1.Label);
plot(run2.Time_s(in2), run2.PhiError_rad(in2), ...
    'Color', orange, 'LineWidth', 1.2, 'DisplayName', run2.Label);
yline(0, 'k:');
hold off; grid on; xlim([tStart tEnd]);
xlabel('Time (s)'); ylabel('Phi error (rad)');
title('Phi error: desired - measured');
legend('Location', 'best');

result.Path1 = run1;
result.Path2 = run2;
result.CommonInterval_s = [tStart tEnd];
result.SameReference = sameReference;
result.Figures = [figTheta figPhi figError];
end

function run = readRun(path)
assert(ischar(path) || (isstring(path) && isscalar(path)), ...
    'Each path must be a character vector or string scalar.');
path = char(path);
if ~isfile(path)
    path = fullfile(fileparts(mfilename('fullpath')), path);
end
assert(isfile(path), 'MAT file not found: %s', path);
saved = load(path, 'out');
assert(isfield(saved, 'out') && isa(saved.out, 'Simulink.SimulationOutput'), ...
    'Expected a Simulink.SimulationOutput named out in %s.', path);
ds = saved.out.logsout;
assert(isa(ds, 'Simulink.SimulationData.Dataset'), ...
    'Expected out.logsout to be a Simulink Dataset in %s.', path);

[t, measuredTheta] = readSignal(ds, 'thetaOrientation');
[phiTime, measuredPhi] = readSignal(ds, 'phiOrientation');
assert(numel(t) == numel(phiTime) && max(abs(t - phiTime)) < 1e-8, ...
    'Theta and phi measurement times differ in %s.', path);
[thetaRefTime, thetaRef] = readSignal(ds, 'des_theta');
[phiRefTime, phiRef] = readSignal(ds, 'des_phi');
[orientationTime, orientationValid] = readSignal(ds, 'orientationValid');
[observableTime, thetaObservable] = readSignal(ds, 'thetaObservable');

desiredTheta = interp1(thetaRefTime, thetaRef, t, 'linear', NaN);
desiredPhi = interp1(phiRefTime, phiRef, t, 'linear', NaN);
validOrientation = interp1(orientationTime, orientationValid, ...
    t, 'previous', 0) > 0.5;
validTheta = validOrientation & interp1(observableTime, ...
    thetaObservable, t, 'previous', 0) > 0.5;
measuredPhi(~validOrientation) = NaN;
measuredTheta(~validTheta) = NaN;

[~, label] = fileparts(path);
run.Path = path;
run.Label = label;
run.Time_s = t;
run.ReferenceThetaTime_s = thetaRefTime;
run.ReferenceTheta_rad = thetaRef;
run.ReferencePhiTime_s = phiRefTime;
run.ReferencePhi_rad = phiRef;
run.DesiredTheta_rad = desiredTheta;
run.DesiredPhi_rad = desiredPhi;
run.MeasuredTheta_rad = measuredTheta;
run.MeasuredPhi_rad = measuredPhi;
run.ThetaError_rad = wrapAngle(desiredTheta - measuredTheta);
run.PhiError_rad = desiredPhi - measuredPhi;
run.Kp = readGain(ds, 'Kp');
run.Kd = readGain(ds, 'Kd');
run.Ki = readGain(ds, 'Ki');
end

function [t, data] = readSignal(ds, name)
names = ds.getElementNames;
index = find(strcmp(names, name));
assert(numel(index) == 1, 'Expected one logged signal named %s.', name);
values = ds.getElement(index).Values;
assert(isa(values, 'timeseries'), '%s must be a timeseries.', name);
t = double(values.Time(:));
raw = double(values.Data);
if ismatrix(raw) && size(raw, 1) == numel(t)
    data = raw;
else
    data = reshape(raw, [], numel(t)).';
end
assert(size(data, 2) == 1, '%s must be scalar.', name);
assert(all(diff(t) > 0), '%s times must increase strictly.', name);
end

function gain = readGain(ds, name)
names = ds.getElementNames;
index = find(strcmp(names, name));
if isempty(index)
    gain = NaN(1, 2);
    return;
end
assert(numel(index) == 1, 'Expected one logged gain signal named %s.', name);
values = ds.getElement(index).Values;
raw = double(values.Data);
if ismatrix(raw) && size(raw, 1) == numel(values.Time)
    gain = raw;
else
    gain = reshape(raw, [], numel(values.Time)).';
end
end

function s = gainText(gain)
s = mat2str(unique(gain, 'rows', 'stable'), 5);
end

function angle = wrapAngle(angle)
angle = atan2(sin(angle), cos(angle));
end

function y = breakWrap(y)
% Do not draw a vertical line across a +/-pi branch cut.
jump = find(abs(diff(y)) > pi);
y(jump + 1) = NaN;
end
