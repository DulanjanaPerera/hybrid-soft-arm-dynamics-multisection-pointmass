function result = plotCircleTracking(source)
%PLOTCIRCLETRACKING Plot one feedback-linearization circle run offline.
%   result = plotCircleTracking(out) uses a Simulink.SimulationOutput.
%   result = plotCircleTracking('circle_multi_rounds_20s.mat') loads a MAT file.
%   With no argument, use the base-workspace out variable if present;
%   otherwise load the newest circle_*.mat beside this function.
%
%   This function only reads logged data and opens figures. It does not
%   load or run a Simulink model, operate NDI/NI hardware, or save files.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(source)
    if evalin('base', 'exist(''out'',''var'')')
        source = evalin('base', 'out');
    else
        files = dir(fullfile(here, 'circle_*.mat'));
        assert(~isempty(files), 'No circle_*.mat file found in %s.', here);
        [~, newest] = max([files.datenum]);
        source = fullfile(files(newest).folder, files(newest).name);
    end
end

if ischar(source) || (isstring(source) && isscalar(source))
    inputFile = char(source);
    if ~isfile(inputFile), inputFile = fullfile(here, inputFile); end
    assert(isfile(inputFile), 'MAT file not found: %s', inputFile);
    saved = load(inputFile, 'out');
    assert(isfield(saved, 'out'), 'MAT file must contain an out variable.');
    out = saved.out;
    sourceName = inputFile;
else
    out = source;
    sourceName = 'Simulink out variable';
end
assert(isa(out, 'Simulink.SimulationOutput'), ...
    'Input must be a Simulink.SimulationOutput or a MAT file containing out.');
ds = out.logsout;

% Use the 0.05 s desired-length samples as the common plotting timeline.
[t, qDesired23] = readSignal(ds, 'sampled des_length', 2);
thetaDesired = alignSignal(ds, 'des_theta', t, 'linear', 1);
phiDesired = alignSignal(ds, 'des_phi', t, 'linear', 1);
thetaMeasured = alignSignal(ds, 'thetaOrientation', t, 'linear', 1);
phiMeasured = alignSignal(ds, 'phiOrientation', t, 'linear', 1);
qMeasured23 = alignSignal(ds, 'fitered_length', t, 'linear', 2);
pressure = alignSignal(ds, 'pressureCommand_bar', t, 'previous', 3);
xyzMeasured = alignSignal(ds, 'position', t, 'linear', 3);
orientationValid = logical(alignSignal(ds, 'orientationValid', t, 'previous', 1));
thetaObservable = logical(alignSignal(ds, 'thetaObservable', t, 'previous', 1));
configurationValid = logical(alignSignal(ds, 'configurationValid', t, 'previous', 1));

% The controller's reduced lengths are [l2,l3]; l1=-l2-l3.
lengthDesired = [-sum(qDesired23, 2), qDesired23];
lengthMeasured = [-sum(qMeasured23, 2), qMeasured23];
lengthMeasured(~orientationValid, :) = NaN;
thetaMeasured(~orientationValid) = NaN;
phiMeasured(~orientationValid) = NaN;
lengthError = lengthDesired - lengthMeasured;

% Theta is periodic. Show an unwrapped trace for lap tracking, but use
% the shortest signed angular difference for the physical theta error.
thetaGood = orientationValid & thetaObservable & ...
    all(isfinite([thetaDesired, thetaMeasured, phiDesired, phiMeasured]), 2);
thetaError = atan2(sin(thetaDesired-thetaMeasured), ...
                   cos(thetaDesired-thetaMeasured));
thetaError(~thetaGood) = NaN;
phiError = phiDesired - phiMeasured;
thetaMeasuredUnwrapped = NaN(size(thetaMeasured));
goodIndex = find(thetaGood);
if ~isempty(goodIndex)
    unwrapped = unwrap(thetaMeasured(goodIndex));
    unwrapped = unwrapped + 2*pi*round( ...
        (thetaDesired(goodIndex(1))-unwrapped(1))/(2*pi));
    thetaMeasuredUnwrapped(goodIndex) = unwrapped;
end

% Desired tip-sensor XYZ from the same constant-curvature geometry used
% by f20260219_2_task2config_withLenExt. The logged position is the
% measured tip-sensor XYZ in the calibrated arm frame.
geometry_m = [0.055, 0.17411, 0.05]; % base offset, module, tip offset
xyzDesired = desiredTipPosition(thetaDesired, phiDesired, geometry_m);
xyzMeasured(~configurationValid, :) = NaN;
xyzError = xyzDesired - xyzMeasured;
xyzErrorNorm = sqrt(sum(xyzError.^2, 2));

result.Source = sourceName;
result.Time_s = t;
result.DesiredTheta_rad = thetaDesired;
result.MeasuredTheta_rad = thetaMeasured;
result.MeasuredThetaUnwrapped_rad = thetaMeasuredUnwrapped;
result.DesiredPhi_rad = phiDesired;
result.MeasuredPhi_rad = phiMeasured;
result.DesiredLengths_m = lengthDesired;
result.MeasuredLengths_m = lengthMeasured;
result.PressureCommand_bar = pressure;
result.DesiredXYZ_m = xyzDesired;
result.MeasuredXYZ_m = xyzMeasured;
result.ThetaErrorWrapped_rad = thetaError;
result.PhiError_rad = phiError;
result.LengthError_m = lengthError;
result.XYZError_m = xyzError;
result.XYZErrorNorm_m = xyzErrorNorm;
result.OrientationValid = orientationValid;
result.ConfigurationValid = configurationValid;

fig = gobjects(7, 1);
actualColor = [0.08, 0.38, 0.76];

fig(1) = figure('Name', '1 Angles', 'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
plot(t, thetaDesired, 'k--', t, thetaMeasuredUnwrapped, ...
    'Color', actualColor, 'LineWidth', 1.2);
ylabel('\theta (rad)'); grid on;
legend('Desired', 'Measured (unwrapped)', 'Location', 'best');
nexttile;
plot(t, phiDesired, 'k--', t, phiMeasured, ...
    'Color', actualColor, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('\phi (rad)'); grid on;
legend('Desired', 'Measured', 'Location', 'best');

fig(2) = figure('Name', '2 Lengths', 'NumberTitle', 'off');
tiledlayout(3, 1, 'TileSpacing', 'compact');
for i = 1:3
    nexttile;
    plot(t, 1000*lengthDesired(:,i), 'k--', ...
         t, 1000*lengthMeasured(:,i), 'Color', actualColor, ...
         'LineWidth', 1.2);
    ylabel(sprintf('l_%d (mm)', i)); grid on;
    if i == 1, legend('Desired', 'Measured', 'Location', 'best'); end
    if i == 3, xlabel('Time (s)'); end
end

fig(3) = figure('Name', '3 Pressure commands', 'NumberTitle', 'off');
plot(t, pressure, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Command (bar)'); grid on;
legend('P1', 'P2', 'P3', 'Location', 'best');

fig(4) = figure('Name', '4 Tip XYZ path', 'NumberTitle', 'off');
plot3(xyzDesired(:,1), xyzDesired(:,2), xyzDesired(:,3), ...
    'k--', 'LineWidth', 1.4); hold on;
plot3(xyzMeasured(:,1), xyzMeasured(:,2), xyzMeasured(:,3), ...
    'Color', actualColor, 'LineWidth', 1.2);
plot3(xyzDesired(1,1), xyzDesired(1,2), xyzDesired(1,3), ...
    'ko', 'MarkerFaceColor', 'k');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m, down)');
legend('Desired', 'Measured', 'Desired start', 'Location', 'best');
axis equal; grid on; view(3); set(gca, 'ZDir', 'reverse');

fig(5) = figure('Name', '5 Angle errors', 'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
plot(t, thetaError, 'LineWidth', 1.2);
ylabel('\theta error (rad)'); grid on;
nexttile;
plot(t, phiError, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('\phi error (rad)'); grid on;

fig(6) = figure('Name', '6 Length errors', 'NumberTitle', 'off');
tiledlayout(3, 1, 'TileSpacing', 'compact');
for i = 1:3
    nexttile;
    plot(t, 1000*lengthError(:,i), 'LineWidth', 1.2);
    ylabel(sprintf('e_{l%d} (mm)', i)); grid on;
    if i == 3, xlabel('Time (s)'); end
end

fig(7) = figure('Name', '7 XYZ errors', 'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
plot(t, 1000*xyzError, 'LineWidth', 1.2);
ylabel('Error (mm)'); grid on;
legend('X', 'Y', 'Z', 'Location', 'best');
nexttile;
plot(t, 1000*xyzErrorNorm, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('XYZ norm (mm)'); grid on;

result.Figures = fig;
end

function [t, x] = readSignal(ds, name, channels)
names = ds.getElementNames;
index = find(strcmp(names, name));
assert(numel(index) == 1, 'Expected one logged signal named %s.', name);
values = ds.get(index).Values;
assert(isa(values, 'timeseries'), 'Signal %s must contain a timeseries.', name);
t = double(values.Time(:));
raw = double(values.Data);
if isvector(raw) && numel(raw) == numel(t)
    x = raw(:);
elseif size(raw, 1) == numel(t)
    x = reshape(raw, numel(t), []);
elseif size(raw, ndims(raw)) == numel(t)
    x = reshape(raw, [], numel(t)).';
else
    error('Cannot align the dimensions of logged signal %s.', name);
end
assert(size(x, 2) == channels, ...
    'Logged signal %s has %d channels; expected %d.', ...
    name, size(x, 2), channels);
end

function x = alignSignal(ds, name, t, method, channels)
[ts, xs] = readSignal(ds, name, channels);
if isequal(ts, t)
    x = xs;
    return;
end
[ts, keep] = unique(ts, 'stable');
xs = xs(keep, :);
assert(numel(ts) > 1, 'Signal %s has too few samples.', name);
assert(t(1) >= ts(1)-0.05 && t(end) <= ts(end)+0.05, ...
    'Signal %s does not cover the desired-length timeline.', name);
x = interp1(ts, xs, t, method, 'extrap');
end

function xyz = desiredTipPosition(theta, phi, L)
ratio = ones(size(phi));
nonzero = abs(phi) > 1e-8;
ratio(nonzero) = sin(phi(nonzero))./phi(nonzero);
rho = (L(2)/2 + L(3))*phi;
rho(nonzero) = L(2)*2*sin(phi(nonzero)/2).^2./phi(nonzero) ...
             + L(3)*sin(phi(nonzero));
z = L(1) + L(2)*ratio + L(3)*cos(phi);
xyz = [rho.*cos(theta), rho.*sin(theta), z];
end
