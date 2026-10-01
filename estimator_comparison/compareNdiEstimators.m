function [summary, results] = compareNdiEstimators(cfg)
%COMPARENDIESTIMATORS Compare two causal NDI length-velocity estimators.
%   [SUMMARY, RESULTS] = COMPARENDIESTIMATORS() reads the six saved
%   forwardModelValidation/*sineWave.mat runs. It does not open COM12,
%   simulate a model, or command the valves.
%
%   Optional CFG fields: SourceDir, OutputDir, TimestampSignal, FilterTau_s,
%   MeasurementStd_m, AccelerationStd_mps2, InitialVelocityStd_mps,
%   ResetGap_s, LowCommandBar, LowCommandMax_s, ShowFigures,
%   MakeFigures, WriteResults, Verbose. Set the last three false for sweeps.
%
%   The NDI reader's lengthChange_m is inferred from orientation under
%   constant-curvature geometry. It is not a direct strain measurement.

if nargin < 1
    cfg = struct();
end
here = fileparts(mfilename('fullpath'));
projectRoot = fileparts(here);
cfg = defaultValue(cfg, 'SourceDir', fullfile(projectRoot, 'forwardModelValidation'));
cfg = defaultValue(cfg, 'OutputDir', fullfile(here, 'results'));
cfg = defaultValue(cfg, 'TimestampSignal', 'pollStart_s');
cfg = defaultValue(cfg, 'FilterTau_s', 0.20);
cfg = defaultValue(cfg, 'MeasurementStd_m', 5e-5);
cfg = defaultValue(cfg, 'AccelerationStd_mps2', 0.003);
cfg = defaultValue(cfg, 'InitialVelocityStd_mps', 0.005);
cfg = defaultValue(cfg, 'ResetGap_s', 0.5);
cfg = defaultValue(cfg, 'LowCommandBar', 0.3);
cfg = defaultValue(cfg, 'LowCommandMax_s', 3.0);
cfg = defaultValue(cfg, 'ShowFigures', false);
cfg = defaultValue(cfg, 'MakeFigures', true);
cfg = defaultValue(cfg, 'WriteResults', true);
cfg = defaultValue(cfg, 'Verbose', true);

assert(cfg.FilterTau_s > 0 && cfg.MeasurementStd_m > 0 && ...
    cfg.AccelerationStd_mps2 > 0 && cfg.ResetGap_s > 0, ...
    'Estimator parameters must be positive.');
trialNames = {'P1_sineWave','P2_sineWave','P3_sineWave', ...
    'P1P2_sineWave','P1P3_sineWave','P2P3_sineWave'};
if (cfg.WriteResults || cfg.MakeFigures) && ~isfolder(cfg.OutputDir)
    mkdir(cfg.OutputDir);
end

results = struct([]);
rows = cell(numel(trialNames), 10);
for j = 1:numel(trialNames)
    file = fullfile(cfg.SourceDir, [trialNames{j} '.mat']);
    assert(isfile(file), 'Recording not found: %s', file);
    loaded = load(file, 'out');
    assert(isfield(loaded, 'out') && isprop(loaded.out, 'logsout'), ...
        'Expected out.logsout in %s', file);
    ds = loaded.out.logsout;

    qSignal = findSignal(ds, {'lengthChnage_m','lengthChange_m'});
    q = lengthRows(qSignal.Values.Data);
    simTime_s = qSignal.Values.Time(:);
    hostSignal = findSignal(ds, {'hostTime'});
    hostTime_s = double(squeeze(hostSignal.Values.Data));
    hostTime_s = hostTime_s(:);
    pollSignal = findSignal(ds, {'pollStart_s'});
    pollStart_s = double(squeeze(pollSignal.Values.Data));
    pollStart_s = pollStart_s(:);
    if strcmp(cfg.TimestampSignal, 'pollStart_s')
        sampleTime_s = pollStart_s;
    elseif strcmp(cfg.TimestampSignal, 'hostTime')
        sampleTime_s = hostTime_s;
    else
        error('TimestampSignal must be pollStart_s or hostTime.');
    end
    validSignal = findSignal(ds, {'orientationValid'});
    valid = logical(squeeze(validSignal.Values.Data));
    valid = valid(:) & all(isfinite(q), 2) & isfinite(sampleTime_s);
    n = size(q, 1);
    assert(numel(simTime_s) == n && numel(hostTime_s) == n && ...
        numel(pollStart_s) == n && ...
        numel(valid) == n, 'NDI signal lengths disagree in %s', file);

    pressureSignal = findSignal(ds, {'des_pressure'});
    pressureBar = double(pressureSignal.Values.Data);
    assert(size(pressureBar, 2) == 3, ...
        'Expected three commanded-pressure channels in %s', file);
    maxPressureBar = interp1(pressureSignal.Values.Time(:), ...
        max(pressureBar, [], 2), simTime_s, 'previous', 'extrap');

    [rawVelocity, filteredVelocity, kalmanQ, kalmanVelocity] = ...
        estimateVelocities(q, sampleTime_s, valid, cfg);
    [fdPredictionMm, kfPredictionMm, predictionPairs] = ...
        predictionErrors(q, sampleTime_s, valid, filteredVelocity, ...
        kalmanQ, kalmanVelocity, cfg.ResetGap_s);

    lowCommand = valid & simTime_s <= cfg.LowCommandMax_s & ...
        maxPressureBar <= cfg.LowCommandBar;
    if nnz(lowCommand) < 20
        fdLowCommandMmS = NaN;
        kfLowCommandMmS = NaN;
    else
        fdLowCommandMmS = vectorRms(filteredVelocity(lowCommand, :))*1000;
        kfLowCommandMmS = vectorRms(kalmanVelocity(lowCommand, :))*1000;
    end

    intervals = diff(sampleTime_s);
    intervals = intervals(isfinite(intervals) & intervals > 0);
    results(j).trial = trialNames{j};
    results(j).sourceFile = file;
    results(j).simTime_s = simTime_s;
    results(j).hostTime_s = hostTime_s;
    results(j).pollStart_s = pollStart_s;
    results(j).sampleTime_s = sampleTime_s;
    results(j).orientationValid = valid;
    results(j).ndiLengths_m = q;
    results(j).rawVelocity_mps = rawVelocity;
    results(j).filteredVelocity_mps = filteredVelocity;
    results(j).kalmanLengths_m = kalmanQ;
    results(j).kalmanVelocity_mps = kalmanVelocity;
    results(j).maxCommandPressure_bar = maxPressureBar;
    results(j).lowCommandMask = lowCommand;

    rows(j, :) = {string(trialNames{j}), nnz(valid), n, ...
        median(intervals), max(intervals), nnz(lowCommand), ...
        fdLowCommandMmS, kfLowCommandMmS, ...
        fdPredictionMm, kfPredictionMm};
    if cfg.MakeFigures
        writeTrialFigure(results(j), cfg.OutputDir, cfg.ShowFigures);
    end
    if cfg.Verbose
        fprintf(['%s: valid %d/%d, %s dt median %.4f s, low-command velocity RMS ' ...
            'FD %.3f / KF %.3f mm/s, one-step error FD %.3f / KF %.3f mm (%d pairs)\n'], ...
            trialNames{j}, nnz(valid), n, cfg.TimestampSignal, median(intervals), ...
            fdLowCommandMmS, kfLowCommandMmS, ...
            fdPredictionMm, kfPredictionMm, predictionPairs);
    end
end

summary = cell2table(rows, 'VariableNames', { ...
    'Trial','ValidSamples','TotalSamples','MedianSampleInterval_s', ...
    'MaxSampleInterval_s','LowCommandSamples', ...
    'FilteredLowCommandVelocityRMS_mm_s', ...
    'KalmanLowCommandVelocityRMS_mm_s', ...
    'FilteredOneStepPredictionRMS_mm', ...
    'KalmanOneStepPredictionRMS_mm'});
if cfg.WriteResults
    writetable(summary, fullfile(cfg.OutputDir, 'estimator_summary.csv'));
    save(fullfile(cfg.OutputDir, 'estimator_comparison.mat'), ...
        'summary', 'results', 'cfg', '-v7.3');
end
if cfg.Verbose
    disp(summary);
    if cfg.WriteResults
        fprintf('Saved figures, CSV, and MAT results in %s\n', cfg.OutputDir);
    end
end
end

function cfg = defaultValue(cfg, name, value)
if ~isfield(cfg, name) || isempty(cfg.(name))
    cfg.(name) = value;
end
end

function sig = findSignal(ds, candidates)
names = getElementNames(ds);
for i = 1:numel(candidates)
    matches = find(strcmp(names, candidates{i}));
    if isscalar(matches)
        sig = ds{matches};
        return
    end
end
error('Expected exactly one signal matching: %s', strjoin(candidates, ', '));
end

function q = lengthRows(data)
data = double(squeeze(data));
if size(data, 1) == 3
    q = data(2:3, :).';
elseif size(data, 2) == 3
    q = data(:, 2:3);
else
    error('Expected three NDI length-change channels.');
end
end

function [rawV, filtV, kfQ, kfV] = estimateVelocities(q, t, valid, cfg)
n = size(q, 1);
rawV = nan(n, 2);
filtV = nan(n, 2);
kfQ = nan(n, 2);
kfV = nan(n, 2);
H = [eye(2), zeros(2)];
R = cfg.MeasurementStd_m^2 * eye(2);
initialP = diag([cfg.MeasurementStd_m^2, cfg.MeasurementStd_m^2, ...
    cfg.InitialVelocityStd_mps^2, cfg.InitialVelocityStd_mps^2]);
initialized = false;
previousValidQ = [NaN, NaN];
previousValidTime = NaN;
previousHostTime = NaN;
velocityFiltered = [0, 0];

for k = 1:n
    if ~initialized
        if ~valid(k)
            continue
        end
        state = [q(k, :).'; 0; 0];
        P = initialP;
        initialized = true;
        previousValidQ = q(k, :);
        previousValidTime = t(k);
        previousHostTime = t(k);
    else
        dt = t(k) - previousHostTime;
        if ~isfinite(dt) || dt <= 0
            filtV(k, :) = velocityFiltered;
            kfQ(k, :) = state(1:2).';
            kfV(k, :) = state(3:4).';
            continue
        end
        previousHostTime = t(k);
        if dt > cfg.ResetGap_s
            % A long software gap invalidates the old velocity estimate.
            state(3:4) = 0;
            if valid(k)
                state(1:2) = q(k, :).';
            end
            P = initialP;
            velocityFiltered = [0, 0];
        else
            F = [eye(2), dt*eye(2); zeros(2), eye(2)];
            G = [0.5*dt^2*eye(2); dt*eye(2)];
            Q = cfg.AccelerationStd_mps2^2 * (G*G.');
            state = F*state;
            P = F*P*F.' + Q;
        end

        if valid(k)
            dtValid = t(k) - previousValidTime;
            if isfinite(dtValid) && dtValid > 0 && dtValid <= cfg.ResetGap_s
                raw = (q(k, :) - previousValidQ)/dtValid;
                weight = dtValid/(cfg.FilterTau_s + dtValid);
                velocityFiltered = velocityFiltered + ...
                    weight*(raw - velocityFiltered);
                rawV(k, :) = raw;
            else
                velocityFiltered = [0, 0];
            end
            previousValidQ = q(k, :);
            previousValidTime = t(k);

            innovation = q(k, :).'-H*state;
            S = H*P*H.' + R;
            gain = (P*H.')/S;
            state = state + gain*innovation;
            I_KH = eye(4) - gain*H;
            P = I_KH*P*I_KH.' + gain*R*gain.';
            P = 0.5*(P + P.');
        end
    end
    filtV(k, :) = velocityFiltered;
    kfQ(k, :) = state(1:2).';
    kfV(k, :) = state(3:4).';
end
end

function [fdMm, kfMm, nPairs] = predictionErrors(q, t, valid, ...
    filtV, kfQ, kfV, maxGap)
fdError = nan(size(q, 1)-1, 2);
kfError = fdError;
for k = 1:size(q, 1)-1
    dt = t(k+1)-t(k);
    if valid(k) && valid(k+1) && dt > 0 && dt <= maxGap && ...
            all(isfinite([filtV(k, :), kfQ(k, :), kfV(k, :)]))
        fdError(k, :) = q(k+1, :) - (q(k, :) + dt*filtV(k, :));
        kfError(k, :) = q(k+1, :) - (kfQ(k, :) + dt*kfV(k, :));
    end
end
nPairs = nnz(all(isfinite(fdError), 2));
fdMm = 1000*vectorRms(fdError);
kfMm = 1000*vectorRms(kfError);
end

function value = vectorRms(x)
if isempty(x) || ~any(all(isfinite(x), 2))
    value = NaN;
    return
end
x = x(all(isfinite(x), 2), :);
value = sqrt(mean(sum(x.^2, 2)));
end

function writeTrialFigure(result, outputDir, showFigures)
visible = 'off';
if showFigures
    visible = 'on';
end
f = figure('Visible', visible, 'Color', 'w', ...
    'Position', [100, 100, 1200, 900]);
if ~showFigures
    cleanup = onCleanup(@() close(f));
end
t = result.sampleTime_s-result.sampleTime_s(find(isfinite(result.sampleTime_s), 1));
q = result.ndiLengths_m*1000;
qKF = result.kalmanLengths_m*1000;
vRaw = result.rawVelocity_mps*1000;
vFD = result.filteredVelocity_mps*1000;
vKF = result.kalmanVelocity_mps*1000;
tl = tiledlayout(f, 3, 1, 'TileSpacing', 'compact');
title(tl, strrep(result.trial, '_', '\_'));
nexttile(tl);
plot(t, q(:, 1), 'k-', t, qKF(:, 1), 'b-');
hold on
plot(t, q(:, 2), '-', 'Color', [0.6 0.6 0.6]);
plot(t, qKF(:, 2), 'r-');
ylabel('Length change (mm)');
legend('NDI l2','KF l2','NDI l3','KF l3','Location','best');
grid on
nexttile(tl);
plot(t, vRaw(:, 1), 'Color', [0.7 0.7 0.7]);
hold on
plot(t, vFD(:, 1), 'b-', t, vKF(:, 1), 'r-');
ylabel('dl2/dt (mm/s)');
legend('Raw difference','Filtered difference','Kalman','Location','best');
grid on
nexttile(tl);
plot(t, vRaw(:, 2), 'Color', [0.7 0.7 0.7]);
hold on
plot(t, vFD(:, 2), 'b-', t, vKF(:, 2), 'r-');
ylabel('dl3/dt (mm/s)');
xlabel('NDI sample timestamp since first poll (s)');
legend('Raw difference','Filtered difference','Kalman','Location','best');
grid on
exportgraphics(f, fullfile(outputDir, [result.trial '.png']), ...
    'Resolution', 160);
end
