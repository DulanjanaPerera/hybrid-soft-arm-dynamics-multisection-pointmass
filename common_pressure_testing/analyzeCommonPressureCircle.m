function report = analyzeCommonPressureCircle(matPath)
%ANALYZECOMMONPRESSURECIRCLE Print a compact report from a saved circle run.
%   analyzeCommonPressureCircle('FB_circle_phi_0.2rad.mat')
%   report = analyzeCommonPressureCircle('C:\full\path\to\run.mat');
%
% The function reads out.logsout only. It does not open a Simulink model,
% operate NI/NDI hardware, create figures, or write analysis files.
% Copy the printed report into the discussion after each run.

folder = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(matPath)
    files = dir(fullfile(folder, 'FB_circle*.mat'));
    assert(~isempty(files), 'No FB_circle*.mat files found in %s.', folder);
    [~, newest] = max([files.datenum]);
    matPath = fullfile(files(newest).folder, files(newest).name);
else
    matPath = char(matPath);
    if ~isfile(matPath), matPath = fullfile(folder, matPath); end
end
assert(isfile(matPath), 'MAT file not found: %s', matPath);

saved = load(matPath, 'out');
assert(isfield(saved, 'out') && isa(saved.out, 'Simulink.SimulationOutput'), ...
    'Expected a Simulink.SimulationOutput variable named out.');
ds = saved.out.logsout;
assert(isa(ds, 'Simulink.SimulationData.Dataset'), ...
    'Expected out.logsout to be a Dataset.');

[t, pressure] = readSignal(ds, 'pressureCommand_bar');
assert(size(pressure, 2) == 3, 'pressureCommand_bar must have three columns.');
desiredPhi = atTime(ds, 'des_phi', t, 'linear');
measuredPhi = atTime(ds, 'phi', t, 'linear');
desiredTheta = atTime(ds, 'des_theta', t, 'linear');
measuredTheta = atTime(ds, 'theta', t, 'linear');
desiredLength = atTime(ds, 'sampled des_length', t, 'previous');
measuredLength = atTime(ds, 'meas_length', t, 'linear');
assert(size(desiredLength, 2) == 2 && size(measuredLength, 2) == 3, ...
    'Expected desired [l2,l3] and measured [l1,l2,l3] lengths.');
measuredLength = measuredLength(:, 2:3);

limited = atTime(ds, 'pressureLimited', t, 'previous') > 0.5;
modelValid = atTime(ds, 'inverseModelValid', t, 'previous') > 0.5;
orientationValid = atTime(ds, 'orientationValid', t, 'previous') > 0.5;
configurationValid = atTime(ds, 'configurationValid', t, 'previous') > 0.5;
thetaObservable = atTime(ds, 'thetaObservable', t, 'previous') > 0.5;
sensorValid = all(atTime(ds, 'sensorValid', t, 'previous') > 0.5, 2);
valid = modelValid & orientationValid & configurationValid & sensorValid ...
    & isfinite(measuredPhi) & isfinite(measuredTheta) ...
    & all(isfinite(measuredLength), 2);

% Detect the circle from the logged reference, including multiple laps.
[refTime, rawDesiredTheta] = readSignal(ds, 'des_theta');
assert(size(rawDesiredTheta, 2) == 1, 'des_theta must be scalar.');
up = find(diff(rawDesiredTheta) > 1e-12);
assert(~isempty(up), 'No increasing circle reference found in des_theta.');
circleStart = refTime(up(1));
circleEnd = refTime(up(end) + 1);
numCircles = round(max(rawDesiredTheta) / (2*pi));
assert(numCircles >= 1 && circleEnd > circleStart, ...
    'Could not identify complete circle timing.');
lapDuration = (circleEnd - circleStart) / numCircles;
angleMismatch = abs(max(rawDesiredTheta) - 2*pi*numCircles);

holdWindow = t >= max(t(1), circleStart - 5) & t < circleStart;
circleWindow = t >= circleStart & t < circleEnd;
holdGood = holdWindow & valid;
circleGood = circleWindow & valid;
thetaGood = circleGood & thetaObservable;
assert(nnz(holdGood) >= 5 && nnz(circleGood) >= 5, ...
    'Too few valid samples in the pre-circle hold or circle window.');

thetaError = atan2(sin(desiredTheta - measuredTheta), ...
    cos(desiredTheta - measuredTheta));
lengthError_mm = 1000 * (desiredLength - measuredLength);
niNames = {'P1_NI_DAQ', 'P2_NI_DAQ', 'P#_NI_DAQ'};
if ~hasSignal(ds, niNames{3}) && hasSignal(ds, 'P3_NI_DAQ')
    niNames{3} = 'P3_NI_DAQ';
end
niAvailable = all(cellfun(@(name) hasSignal(ds, name), niNames));
if niAvailable
    niPath = [atTime(ds, niNames{1}, t, 'previous'), ...
        atTime(ds, niNames{2}, t, 'previous'), ...
        atTime(ds, niNames{3}, t, 'previous')];
    niMismatch = max(abs(pressure - niPath), [], 'all');
end

lines = {sprintf('COMMON-PRESSURE CIRCLE REPORT: %s', matPath)};
if numCircles == 1, circleWord = 'circle'; else, circleWord = 'circles'; end
lines{end+1} = sprintf('Logged reference: phi peak %.4f rad; %d %s, %.2f-%.2f s, %.2f s/circle', ...
    max(desiredPhi), numCircles, circleWord, circleStart, circleEnd, lapDuration);
if hasSignal(ds, 'Kp') && hasSignal(ds, 'Kd')
    kp = gainDuring(ds, 'Kp', circleStart, circleEnd);
    kd = gainDuring(ds, 'Kd', circleStart, circleEnd);
    lines{end+1} = sprintf('Logged gains: Kp %s; Kd %s', ...
        gainSummary(kp), gainSummary(kd));
else
    lines{end+1} = 'Logged gains: Kp/Kd unavailable in this MAT file';
end
if angleMismatch > 0.02
    lines{end+1} = sprintf('CHECK: theta maximum differs from %d*2*pi by %.3f rad', ...
        numCircles, angleMismatch);
end
lines{end+1} = sprintf('Pre-circle hold %.2f-%.2f s: desired/actual phi %.4f/%.4f rad', ...
    circleStart - 5, circleStart, mean(desiredPhi(holdGood)), ...
    mean(measuredPhi(holdGood)));
lines{end+1} = sprintf('  desired l2,l3 %s mm; actual %s mm; error %s mm', ...
    vec(mean(1000*desiredLength(holdGood,:), 1)), ...
    vec(mean(1000*measuredLength(holdGood,:), 1)), ...
    vec(mean(lengthError_mm(holdGood,:), 1)));
holdPressure = mean(pressure(holdGood,:), 1);
lines{end+1} = sprintf('  pressure P1,P2,P3 %s bar; P2-P1/P3-P1 %s bar', ...
    vec(holdPressure), vec(holdPressure(2:3)-holdPressure(1)));
if hasSignal(ds, 'desired acceleration force') ...
        && hasSignal(ds, 'gravity and stiffness force')
    accelerationForce = atTime(ds, 'desired acceleration force', t, 'linear');
    staticForce = atTime(ds, 'gravity and stiffness force', t, 'linear');
    lines{end+1} = sprintf('  force mean: acceleration %s N; gravity+stiffness %s N', ...
        vec(mean(accelerationForce(holdGood,:), 1)), ...
        vec(mean(staticForce(holdGood,:), 1)));
    lines{end+1} = sprintf('  circle mean |force|: acceleration %s N; gravity+stiffness %s N', ...
        vec(mean(abs(accelerationForce(circleGood,:)), 1)), ...
        vec(mean(abs(staticForce(circleGood,:)), 1)));
end
if hasSignal(ds, 'Coriolist force')
    dampingForce = atTime(ds, 'Coriolist force', t, 'linear');
    lines{end+1} = sprintf('  circle mean |Coriolis+damping force| %s N', ...
        vec(mean(abs(dampingForce(circleGood,:)), 1)));
end

lines{end+1} = sprintf(['Circle quality: valid %d/%d; model %d/%d; NDI orientation ' ...
    '%d/%d; sensor pair %d/%d; configuration %d/%d; theta observable %d/%d; limited %d/%d'], ...
    nnz(circleGood), nnz(circleWindow), nnz(modelValid(circleWindow)), nnz(circleWindow), ...
    nnz(orientationValid(circleWindow)), nnz(circleWindow), ...
    nnz(sensorValid(circleWindow)), nnz(circleWindow), ...
    nnz(configurationValid(circleWindow)), nnz(circleWindow), ...
    nnz(thetaGood), nnz(circleWindow), nnz(limited(circleWindow)), nnz(circleWindow));
lines{end+1} = sprintf('  phi desired/actual mean %.4f/%.4f rad; actual range %.4f..%.4f; phi MAE %.4f rad', ...
    mean(desiredPhi(circleGood)), mean(measuredPhi(circleGood)), ...
    min(measuredPhi(circleGood)), max(measuredPhi(circleGood)), ...
    mean(abs(desiredPhi(circleGood)-measuredPhi(circleGood))));
if nnz(thetaGood) >= 2
    actualAngle = unwrap(measuredTheta(thetaGood));
    measuredTurns = (actualAngle(end) - actualAngle(1)) / (2*pi);
    lines{end+1} = sprintf('  theta wrapped MAE %.4f rad; measured angle change %.3f turns', ...
        mean(abs(thetaError(thetaGood))), measuredTurns);
end
lines{end+1} = sprintf('  length RMSE l2,l3 %s mm; pressure range P1,P2,P3 %s to %s bar', ...
    vec(sqrt(mean(lengthError_mm(circleGood,:).^2, 1))), ...
    vec(min(pressure(circleWindow,:), [], 1)), ...
    vec(max(pressure(circleWindow,:), [], 1)));
circlePressure = pressure(circleWindow,:);
commonPressure = mean(circlePressure, 2);
sPressure = (circlePressure(:,2)+circlePressure(:,3))/2-circlePressure(:,1);
aPressure = (circlePressure(:,3)-circlePressure(:,2))/2;
lines{end+1} = sprintf('  mean common/s/a %.4f/%.4f/%.4f bar; s range %.4f..%.4f; a range %.4f..%.4f', ...
    mean(commonPressure), mean(sPressure), mean(aPressure), ...
    min(sPressure), max(sPressure), min(aPressure), max(aPressure));
if niAvailable
    lines{end+1} = sprintf('  max controller-to-NI-path command difference %.6g bar', ...
        niMismatch);
else
    lines{end+1} = '  NI-path command logs incomplete; no command-path comparison';
end

% A single-lap row duplicates the overall circle summary; print lap rows
% only when they can show differences between repeated circles.
if numCircles > 1
for lap = 1:numCircles
    lapStart = circleStart + (lap-1)*lapDuration;
    lapEnd = circleStart + lap*lapDuration;
    inLap = t >= lapStart & t < lapEnd;
    good = inLap & valid;
    goodTheta = good & thetaObservable;
    if nnz(good) < 2
        lines{end+1} = sprintf('Lap %d (%.1f-%.1f s): insufficient valid data', ...
            lap, lapStart, lapEnd);
        continue;
    end
    if nnz(goodTheta) >= 2
        thetaMAE = mean(abs(thetaError(goodTheta)));
        angle = unwrap(measuredTheta(goodTheta));
        actualLapTurns = (angle(end)-angle(1))/(2*pi);
    else
        thetaMAE = NaN;
        actualLapTurns = NaN;
    end
    lines{end+1} = sprintf(['Lap %d (%.1f-%.1f s): actual phi %.4f rad; ' ...
        'phi MAE %.4f; theta MAE %.3f; measured turns %.3f; ' ...
        'l2/l3 RMSE %s mm; limited %d'], ...
        lap, lapStart, lapEnd, mean(measuredPhi(good)), ...
        mean(abs(desiredPhi(good)-measuredPhi(good))), thetaMAE, ...
        actualLapTurns, vec(sqrt(mean(lengthError_mm(good,:).^2, 1))), ...
        nnz(limited(inLap)));
end
end

lines{end+1} = ['Interpretation: NI-path values are logged commands, not ' ...
    'measured chamber pressure; lengths are inferred from NDI orientation.'];
report = strjoin(lines, newline);
fprintf('%s\n', report);
end

function tf = hasSignal(ds, name)
tf = nnz(strcmp(ds.getElementNames, name)) == 1;
end

function x = atTime(ds, name, t, method)
[tt, xx] = readSignal(ds, name);
assert(tt(1) <= t(1)+1e-6 && tt(end) >= t(end)-1e-6, ...
    'Incomplete coverage for logged signal %s.', name);
if isequal(tt, t)
    x = xx;
else
    x = interp1(tt, xx, t, method, 'extrap');
end
end

function [t, x] = readSignal(ds, name)
indices = find(strcmp(ds.getElementNames, name));
assert(numel(indices) == 1, 'Need exactly one logged signal named %s.', name);
values = ds.getElement(indices).Values;
t = double(values.Time(:));
data = double(values.Data);
if ndims(data) >= 3 && size(data, ndims(data)) == numel(t)
    x = reshape(data, [], numel(t)).';
elseif size(data, 1) == numel(t)
    x = reshape(data, numel(t), []);
else
    error('Unexpected data shape for %s: %s.', name, mat2str(size(data)));
end
assert(all(diff(t) > 0), 'Signal %s has duplicate/nonincreasing times.', name);
end

function s = vec(x)
s = ['[', sprintf(' %.4f', x), ' ]'];
end

function s = gainSummary(x)
finiteRows = all(isfinite(x), 2);
if ~any(finiteRows)
    s = 'nonfinite';
    return;
end

lo = min(x(finiteRows,:), [], 1);
hi = max(x(finiteRows,:), [], 1);
if all(abs(hi-lo) < 1e-9)
    s = vec(lo);
else
    s = sprintf('%s to %s (varied)', vec(lo), vec(hi));
end
end

function x = gainDuring(ds, name, windowStart, windowEnd)
% Gains may be logged only on changes, so carry the last pre-window value.
[t, values] = readSignal(ds, name);
preceding = find(t <= windowStart, 1, 'last');
assert(~isempty(preceding), 'No %s gain value at circle start.', name);
inWindow = t > windowStart & t <= windowEnd;
x = [values(preceding,:); values(inWindow,:)];
end
