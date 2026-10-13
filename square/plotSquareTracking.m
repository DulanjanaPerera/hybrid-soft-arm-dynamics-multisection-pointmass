function figures = plotSquareTracking(source, rotationAngle_rad)
%PLOTSQUARETRACKING Plot a square-tracking run using logged XYZ positions.
%   figures = plotSquareTracking(out)
%   figures = plotSquareTracking(out.logsout)
%   figures = plotSquareTracking('square/initial_test.mat')
%   figures = plotSquareTracking() uses out in the base workspace.
%   figures = plotSquareTracking(source, 0.2903) rotates desired XYZ about
%   +Z by 0.2903 rad before plotting it against measured XYZ. A positive
%   angle rotates +X toward +Y; use a negative angle for the opposite way.
%
% Reuses the five angle, error, length, force, and pressure figures from
% plotCommonPressureTracking. Figure 6 uses logged des_position and
% meas_position directly; it does not reconstruct desired XYZ from angles.
% No files are written, and no model or hardware is run.

if nargin < 2 || isempty(rotationAngle_rad)
    rotationAngle_rad = 0;
end
assert(isnumeric(rotationAngle_rad) && isreal(rotationAngle_rad) && ...
    isscalar(rotationAngle_rad) && isfinite(rotationAngle_rad), ...
    'rotationAngle_rad must be a finite, real scalar in radians.');

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

[~, desiredXYZ] = readXYZ(ds, 'des_position');
[measuredTime, measuredXYZ] = readXYZ(ds, 'meas_position');
desiredXY = desiredXYZ(:,1:2);
cosAngle = cos(rotationAngle_rad);
sinAngle = sin(rotationAngle_rad);
desiredXYZ(:,1) = cosAngle*desiredXY(:,1) - sinAngle*desiredXY(:,2);
desiredXYZ(:,2) = sinAngle*desiredXY(:,1) + cosAngle*desiredXY(:,2);
validTime = ds.getElement('configurationValid').Values.Time;
validData = double(ds.getElement('configurationValid').Values.Data(:));
valid = interp1(validTime(:), validData, measuredTime, 'previous', 0) > 0.5;
measuredXYZ(~valid, :) = NaN;

% Find the sibling plotter without requiring the caller to change folders.
oldPath = path;
restorePath = onCleanup(@() path(oldPath)); %#ok<NASGU>
commonFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'common_pressure_testing');
addpath(commonFolder);
figures = plotCommonPressureTracking(ds);

clf(figures(6));
figure(figures(6));
desiredLabel = 'Desired';
plotTitle = 'Logged desired and measured tip XYZ';
if rotationAngle_rad ~= 0
    desiredLabel = 'Desired (rotated)';
    plotTitle = sprintf('%s (desired yaw %+0.4f rad)', ...
        plotTitle, rotationAngle_rad);
end
plot3(desiredXYZ(:,1), desiredXYZ(:,2), desiredXYZ(:,3), ...
    'k--', 'LineWidth', 1.4, 'DisplayName', desiredLabel);
hold on;
plot3(measuredXYZ(:,1), measuredXYZ(:,2), measuredXYZ(:,3), ...
    'Color', [0.08 0.38 0.76], 'LineWidth', 1.2, ...
    'DisplayName', 'Measured');
plot3(desiredXYZ(1,1), desiredXYZ(1,2), desiredXYZ(1,3), ...
    'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Desired start');
hold off;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m, down)');
title(plotTitle);
legend('Location', 'best');
axis equal; grid on; view(3); set(gca, 'ZDir', 'reverse');
rotate3d(figures(6), 'on');
end

function [t, xyz] = readXYZ(ds, name)
index = find(strcmp(ds.getElementNames, name));
assert(numel(index) == 1, 'Need exactly one logged signal named %s.', name);
values = ds.getElement(index).Values;
assert(isa(values, 'timeseries'), '%s must be a timeseries.', name);
t = double(values.Time(:));
raw = double(values.Data);
if size(raw, 1) == numel(t)
    xyz = reshape(raw, numel(t), []);
elseif size(raw, ndims(raw)) == numel(t)
    xyz = reshape(raw, [], numel(t)).';
else
    error('Unexpected data shape for %s.', name);
end
assert(size(xyz, 2) == 3, '%s must have three XYZ channels.', name);
assert(all(diff(t) > 0), '%s times must increase strictly.', name);
end
