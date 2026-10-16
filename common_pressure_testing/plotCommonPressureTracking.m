function figures = plotCommonPressureTracking(source)
%PLOTCOMMONPRESSURETRACKING Plot a saved common-pressure tracking run.
%   plotCommonPressureTracking(out) reads out.logsout.
%   plotCommonPressureTracking(out.logsout) reads the dataset directly.
%   plotCommonPressureTracking('run.mat') loads its saved out variable.
%   With no input, use out from the base workspace.
%
% Opens six figures and writes no files. It does not run Simulink or hardware.
% Pressures are logged commands, not measured chamber pressures. Measured
% lengths are inferred from NDI orientation, not independent strain sensors.

if nargin < 1 || isempty(source)
    assert(evalin('base', 'exist(''out'',''var'')') == 1, ...
        'Pass out, out.logsout, or a MAT path, or define out in the base workspace.');
    source = evalin('base', 'out');
end

if ischar(source) || (isstring(source) && isscalar(source))
    file = char(source);
    if ~isfile(file)
        file = fullfile(fileparts(mfilename('fullpath')), file);
    end
    assert(isfile(file), 'MAT file not found: %s', file);
    saved = load(file, 'out');
    assert(isfield(saved, 'out'), 'MAT file must contain out.');
    source = saved.out;
end

if isa(source, 'Simulink.SimulationOutput')
    ds = source.logsout;
else
    ds = source;
end
assert(isa(ds, 'Simulink.SimulationData.Dataset'), ...
    'Input must be out, out.logsout, or a MAT file containing out.');

% Use the desired reduced-length samples as the common 0.05 s timeline.
[t, desired23] = readSignal(ds, 'sampled des_length', 2);
thetaDesired = alignSignal(ds, 'des_theta', t, 'linear', 1);
phiDesired = alignSignal(ds, 'des_phi', t, 'linear', 1);
thetaMeasured = alignSignal(ds, 'theta', t, 'linear', 1);
phiMeasured = alignSignal(ds, 'phi', t, 'linear', 1);
lengthMeasured = alignSignal(ds, 'meas_length', t, 'linear', 3);
pressure = alignSignal(ds, 'pressureCommand_bar', t, 'previous', 3);
xyzMeasured = alignSignal(ds, 'position', t, 'linear', 3);
orientationValid = alignSignal(ds, 'orientationValid', t, 'previous', 1) > 0.5;
thetaObservable = alignSignal(ds, 'thetaObservable', t, 'previous', 1) > 0.5;
configurationValid = alignSignal(ds, 'configurationValid', t, 'previous', 1) > 0.5;
isPAC = any(strcmp(ds.getElementNames, 'tauRequested_N'));
if isPAC
    modelValid = alignSignal(ds, 'modelValid', t, 'previous', 1) > 0.5;
else
    modelValid = alignSignal(ds, 'inverseModelValid', t, 'previous', 1) > 0.5;
end

% The backbone constraint gives l1 = -l2-l3.
lengthDesired = [-sum(desired23, 2), desired23];
lengthMeasured(~orientationValid, :) = NaN;
phiMeasured(~orientationValid) = NaN;

% Theta is undefined close to straight. Wrap only the displayed traces;
% retain NaN gaps where NDI measurements were invalid or unobservable.
thetaGood = orientationValid & thetaObservable & ...
    isfinite(thetaDesired) & isfinite(thetaMeasured);
thetaDesiredPlot = wrappedForPlot(thetaDesired);
thetaMeasured(~thetaGood) = NaN;
thetaMeasuredPlot = wrappedForPlot(thetaMeasured);
thetaError = atan2(sin(thetaDesired-thetaMeasured), ...
    cos(thetaDesired-thetaMeasured));
thetaError(~thetaGood) = NaN;
phiError = phiDesired-phiMeasured;

% Desired tip-sensor XYZ uses the calibrated arm-frame geometry and the
% constant-curvature model used by the earlier circle plotting workflow.
geometry_m = [0.055, 0.17411, 0.05]; % base, module, tip sensor offsets
xyzDesired = desiredTipPosition(thetaDesired, phiDesired, geometry_m);
xyzMeasured(~configurationValid, :) = NaN;

if isPAC
    massForce = alignSignal(ds, 'inertialForce_N', t, 'linear', 2);
    velocityForce = alignSignal(ds, 'coriolisForce_N', t, 'linear', 2) + ...
        alignSignal(ds, 'dampingForce_N', t, 'linear', 2);
    staticForce = alignSignal(ds, 'staticForce_N', t, 'linear', 2);
    tauDesired = alignSignal(ds, 'tauRequested_N', t, 'linear', 2);
    forceLabels = {'M a', '(C+D) v', 'Static force', 'Requested total'};
else
    massForce = alignSignal(ds, 'desired acceleration force', t, 'linear', 2);
    velocityForce = alignSignal(ds, 'Coriolist force', t, 'linear', 2);
    staticForce = alignSignal(ds, 'gravity and stiffness force', t, 'linear', 2);
    tauDesired = alignSignal(ds, 'tauDesired_N', t, 'linear', 2);
    forceLabels = {'M v', '(C+D) dq', 'G+Kq', 'Requested total'};
end
massForce(~modelValid, :) = NaN;
velocityForce(~modelValid, :) = NaN;
staticForce(~modelValid, :) = NaN;
tauDesired(~modelValid, :) = NaN;

figures = gobjects(6, 1);
measuredColor = [0.08 0.38 0.76];

figures(1) = figure('Name', 'Desired and measured angles', ...
    'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
plot(t, thetaDesiredPlot, 'k--', 'LineWidth', 1.2); hold on;
plot(t, thetaMeasuredPlot, 'Color', measuredColor, 'LineWidth', 1.2);
ylabel('\theta (rad, wrapped)'); ylim([-pi pi]); grid on;
legend('Desired', 'Measured', 'Location', 'best');
nexttile;
plot(t, phiDesired, 'k--', 'LineWidth', 1.2); hold on;
plot(t, phiMeasured, 'Color', measuredColor, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('\phi (rad)'); grid on;
legend('Desired', 'Measured', 'Location', 'best');

figures(2) = figure('Name', 'Theta and phi errors', ...
    'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
plot(t, thetaError, 'LineWidth', 1.2);
ylabel('Wrapped \theta error (rad)'); grid on;
nexttile;
plot(t, phiError, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('\phi error (rad)'); grid on;

figures(3) = figure('Name', 'Desired and measured lengths', ...
    'NumberTitle', 'off');
tiledlayout(3, 1, 'TileSpacing', 'compact');
for k = 1:3
    nexttile;
    plot(t, 1000*lengthDesired(:,k), 'k--', 'LineWidth', 1.2); hold on;
    plot(t, 1000*lengthMeasured(:,k), 'Color', measuredColor, ...
        'LineWidth', 1.2);
    ylabel(sprintf('l_%d (mm)', k)); grid on;
    if k == 1
        legend('Desired', 'Measured', 'Location', 'best');
    end
    if k == 3, xlabel('Time (s)'); end
end

figures(4) = figure('Name', 'Inverse-dynamics forces', ...
    'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact');
for k = 1:2
    nexttile;
    plot(t, massForce(:,k), 'LineWidth', 1.1); hold on;
    plot(t, velocityForce(:,k), 'LineWidth', 1.1);
    plot(t, staticForce(:,k), 'LineWidth', 1.1);
    plot(t, tauDesired(:,k), 'k--', 'LineWidth', 1.2);
    ylabel(sprintf('q_%d force (N)', k+1)); grid on;
    if k == 1
        legend(forceLabels, 'Location', 'best');
    else
        xlabel('Time (s)');
    end
end

figures(5) = figure('Name', 'Pressure commands', ...
    'NumberTitle', 'off');
plot(t, pressure, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Command (bar)'); grid on;
legend('P1', 'P2', 'P3', 'Location', 'best');

figures(6) = figure('Name', 'Desired and measured tip XYZ', ...
    'NumberTitle', 'off');
plot3(xyzDesired(:,1), xyzDesired(:,2), xyzDesired(:,3), ...
    'k--', 'LineWidth', 1.4); hold on;
plot3(xyzMeasured(:,1), xyzMeasured(:,2), xyzMeasured(:,3), ...
    'Color', measuredColor, 'LineWidth', 1.2);
plot3(xyzDesired(1,1), xyzDesired(1,2), xyzDesired(1,3), ...
    'ko', 'MarkerFaceColor', 'k');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m, down)');
legend('Desired', 'Measured', 'Desired start', 'Location', 'best');
axis equal; grid on; view(3); set(gca, 'ZDir', 'reverse');
end

function y = wrappedForPlot(x)
y = atan2(sin(x), cos(x));
% A wrap crossing is a plot discontinuity, not a physical 2*pi jump.
y([false; abs(diff(y)) > pi]) = NaN;
end

function [t, x] = readSignal(ds, name, channels)
indices = find(strcmp(ds.getElementNames, name), 1, 'first');
assert(~isempty(indices), 'Need a logged signal named %s.', name);
values = ds.getElement(indices).Values;
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
    error('Unexpected data shape for signal %s.', name);
end
assert(size(x, 2) == channels, ...
    'Signal %s has %d channels; expected %d.', name, size(x,2), channels);
assert(all(diff(t) > 0), 'Signal %s has duplicate/nonincreasing times.', name);
end

function x = alignSignal(ds, name, t, method, channels)
[ts, values] = readSignal(ds, name, channels);
assert(ts(1) <= t(1)+1e-6 && ts(end) >= t(end)-1e-6, ...
    'Signal %s does not cover the desired-length timeline.', name);
if isequal(ts, t)
    x = values;
else
    x = interp1(ts, values, t, method, 'extrap');
end
end

function xyz = desiredTipPosition(theta, phi, L)
ratio = ones(size(phi));
rho = (L(2)/2 + L(3))*phi;
nonzero = abs(phi) > 1e-8;
ratio(nonzero) = sin(phi(nonzero))./phi(nonzero);
rho(nonzero) = L(2)*2*sin(phi(nonzero)/2).^2./phi(nonzero) ...
    + L(3)*sin(phi(nonzero));
z = L(1) + L(2)*ratio + L(3)*cos(phi);
xyz = [rho.*cos(theta), rho.*sin(theta), z];
end
