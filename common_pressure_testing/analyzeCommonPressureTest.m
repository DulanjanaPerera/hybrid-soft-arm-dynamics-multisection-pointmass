function results = analyzeCommonPressureTest(matPath, outputDir)
%ANALYZECOMMONPRESSURETEST Offline analysis of the common-pressure test.
%   R = analyzeCommonPressureTest() reads test_300s.mat next to this file.
%   R = analyzeCommonPressureTest(MATPATH, OUTPUTDIR) selects other paths.
%   This function loads saved logs only; it does not open or run Simulink.

folder = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(matPath)
    matPath = fullfile(folder, 'test_300s.mat');
end
if nargin < 2 || isempty(outputDir)
    [~, runName] = fileparts(matPath);
    outputDir = fullfile(folder, 'results', runName);
end
assert(isfile(matPath), 'Missing MAT-file: %s', matPath);
if ~isfolder(outputDir)
    mkdir(outputDir);
end

loaded = load(matPath, 'out');
assert(isfield(loaded, 'out') && isa(loaded.out, 'Simulink.SimulationOutput'), ...
    'Expected a Simulink.SimulationOutput variable named out.');
ds = loaded.out.logsout;
assert(isa(ds, 'Simulink.SimulationData.Dataset'), ...
    'Expected out.logsout to be a Dataset.');

[t, measuredLength] = getSignal(ds, 'meas_length');
assert(size(measuredLength, 2) == 3, 'meas_length must have three columns.');
data.t_s = t;
data.measuredLength_m = measuredLength;
data.modelLength_m = atTime(ds, 'comp_lengths', t, 'linear', '/MATLAB Function1');
data.command_bar = atTime(ds, 'des_pressure', t, 'linear');
data.effectiveModel_bar = atTime(ds, 'pressureUsed_bar', t, 'linear');
data.regulatorPath_bar = [ ...
    atTime(ds, 'P1 Regulator', t, 'linear'), ...
    atTime(ds, 'P2 Regulator', t, 'linear'), ...
    atTime(ds, 'P3 Regulator', t, 'linear')];
data.phaseID = round(atTime(ds, 'phaseID', t, 'previous'));
data.measuredPhi_rad = atTime(ds, 'phi', t, 'linear');
data.measuredTheta_rad = atTime(ds, 'theta', t, 'linear');
data.modelPhi_rad = atTime(ds, 'comp_phi', t, 'linear', '/MATLAB Function1');
data.position_m = atTime(ds, 'position', t, 'linear');
data.orientationValid = atTime(ds, 'orientationValid', t, 'previous', '/Sensor') > 0.5;
data.sensorValid = all(atTime(ds, 'sensorValid', t, 'previous') > 0.5, 2);
data.configurationValid = atTime(ds, 'configurationValid', t, 'previous') > 0.5;
data.pollStart_s = atTime(ds, 'pollStart_s', t, 'linear');
data.hostTime_s = atTime(ds, 'hostTime_s', t, 'linear');
[~, killState] = getSignal(ds, 'kitt_switch');
data.loggedKillState = killState(1); % Raw dashboard state; not a pressure measurement.

assert(size(data.modelLength_m, 2) == 3 && size(data.command_bar, 2) == 3 ...
    && size(data.position_m, 2) == 3, 'Unexpected vector signal width.');
data.valid = data.orientationValid & data.sensorValid & data.configurationValid ...
    & all(isfinite(data.measuredLength_m), 2) ...
    & isfinite(data.measuredPhi_rad) ...
    & all(isfinite(data.position_m), 2);

% This is the same signed bending projection used by the forward block.
% Unlike phi, it preserves the sign when the P1-versus-P2/P3 bend reverses.
radius_m = 0.013; % PMA radial offset from backbone; not pressure-area radius.
forcePerBar_N = pi*(0.013/2)^2*1e5;
data.measuredU_rad = (data.measuredLength_m(:, 2) ...
    + data.measuredLength_m(:, 3))/radius_m;
data.measuredV_rad = (data.measuredLength_m(:, 3) ...
    - data.measuredLength_m(:, 2))/(sqrt(3)*radius_m);
data.modelU_rad = (data.modelLength_m(:, 2) ...
    + data.modelLength_m(:, 3))/radius_m;
data.modelV_rad = (data.modelLength_m(:, 3) ...
    - data.modelLength_m(:, 2))/(sqrt(3)*radius_m);

phaseSummary = summarizePhases(data, 5.0);
conditionSummary = summarizeConditions(phaseSummary, radius_m, forcePerBar_N);
results.sourceMat = char(matPath);
results.outputDir = char(outputDir);
results.assumptions.radius_m = radius_m;
results.assumptions.forcePerBar_N = forcePerBar_N;
results.assumptions.tailWindow_s = 5.0;
results.data = data;
results.phaseSummary = phaseSummary;
results.conditionSummary = conditionSummary;

writetable(phaseSummary, fullfile(outputDir, 'phase_summary.csv'));
writetable(conditionSummary, fullfile(outputDir, 'condition_summary.csv'));
writeReport(results, fullfile(outputDir, 'report.txt'));
makePlots(results, outputDir);
save(fullfile(outputDir, 'analysis_results.mat'), 'results', '-v7');
end

function [t, x] = getSignal(ds, name, blockSuffix)
if nargin < 3
    blockSuffix = '';
end
names = ds.getElementNames;
candidates = find(strcmp(names, name));
matches = zeros(0, 1);
for i = candidates(:)'
    element = ds.get(i);
    path = char(element.BlockPath.getBlock(1));
    if isempty(blockSuffix) || endsWith(path, blockSuffix)
        matches(end+1, 1) = i; %#ok<AGROW>
    end
end
assert(numel(matches) == 1, ...
    'Expected one logged %s signal matching block suffix "%s"; found %d.', ...
    name, blockSuffix, numel(matches));
values = ds.get(matches).Values;
assert(isa(values, 'timeseries'), 'Expected timeseries for %s.', name);
t = double(values.Time(:));
raw = double(values.Data);
if size(raw, 1) == numel(t) && ndims(raw) <= 2
    x = reshape(raw, numel(t), []);
else
    x = reshape(raw, [], numel(t)).';
end
assert(size(x, 1) == numel(t), 'Time/data size mismatch for %s.', name);
end

function x = atTime(ds, name, t, method, blockSuffix)
if nargin < 5
    blockSuffix = '';
end
[sourceTime, source] = getSignal(ds, name, blockSuffix);
[sourceTime, keep] = unique(sourceTime, 'last');
source = source(keep, :);
if numel(sourceTime) == 1
    x = repmat(source, numel(t), 1);
else
    x = interp1(sourceTime, source, t, method, NaN);
end
end

function tableOut = summarizePhases(data, tailWindow_s)
ids = unique(data.phaseID(isfinite(data.phaseID) & data.phaseID > 0));
n = numel(ids);
phaseID = zeros(n, 1);
kind = strings(n, 1);
start_s = NaN(n, 1); end_s = NaN(n, 1);
validFraction = NaN(n, 1); nTail = zeros(n, 1);
common_bar = NaN(n, 1); delta2_bar = NaN(n, 1); delta3_bar = NaN(n, 1);
regulatorError_bar = NaN(n, 1);
measuredU_rad = NaN(n, 1); measuredV_rad = NaN(n, 1);
measuredPhi_rad = NaN(n, 1); modelU_rad = NaN(n, 1);
modelV_rad = NaN(n, 1); modelPhi_rad = NaN(n, 1);
tailURange_rad = NaN(n, 1); tailUSlope_radPerS = NaN(n, 1);
measuredLength_mm = NaN(n, 3); modelLength_mm = NaN(n, 3);
position_mm = NaN(n, 3);

for j = 1:n
    phaseID(j) = ids(j);
    inPhase = data.phaseID == ids(j);
    start_s(j) = min(data.t_s(inPhase));
    end_s(j) = max(data.t_s(inPhase));
    validFraction(j) = mean(data.valid(inPhase));
    tail = inPhase & data.t_s >= end_s(j) - tailWindow_s;
    validTail = tail & data.valid;
    nTail(j) = nnz(validTail);
    command = median(data.command_bar(tail, :), 1, 'omitnan');
    common_bar(j) = mean(command);
    delta2_bar(j) = command(2) - command(1);
    delta3_bar(j) = command(3) - command(1);
    regulatorError_bar(j) = max(abs(data.regulatorPath_bar(tail, :) ...
        - data.command_bar(tail, :)), [], 'all');
    if common_bar(j) < 0.1
        kind(j) = "zero";
    elseif abs(delta2_bar(j)) < 0.05 && abs(delta3_bar(j)) < 0.05
        kind(j) = "neutral";
    elseif delta2_bar(j) > 0.05 && delta3_bar(j) > 0.05
        kind(j) = "positive";
    elseif delta2_bar(j) < -0.05 && delta3_bar(j) < -0.05
        kind(j) = "negative";
    else
        kind(j) = "mixed";
    end
    if nTail(j) > 0
        measuredU_rad(j) = mean(data.measuredU_rad(validTail), 'omitnan');
        tailU = data.measuredU_rad(validTail);
        tailURange_rad(j) = max(tailU) - min(tailU);
        if nTail(j) >= 3
            tailTime = data.t_s(validTail);
            slopeFit = polyfit(tailTime - mean(tailTime), tailU, 1);
            tailUSlope_radPerS(j) = slopeFit(1);
        end
        measuredV_rad(j) = mean(data.measuredV_rad(validTail), 'omitnan');
        measuredPhi_rad(j) = mean(data.measuredPhi_rad(validTail), 'omitnan');
        modelU_rad(j) = mean(data.modelU_rad(validTail), 'omitnan');
        modelV_rad(j) = mean(data.modelV_rad(validTail), 'omitnan');
        modelPhi_rad(j) = mean(data.modelPhi_rad(validTail), 'omitnan');
        measuredLength_mm(j, :) = mean(data.measuredLength_m(validTail, :), 1, 'omitnan')*1e3;
        modelLength_mm(j, :) = mean(data.modelLength_m(validTail, :), 1, 'omitnan')*1e3;
        position_mm(j, :) = mean(data.position_m(validTail, :), 1, 'omitnan')*1e3;
    end
end

tableOut = table(phaseID, kind, start_s, end_s, validFraction, nTail, ...
    common_bar, delta2_bar, delta3_bar, regulatorError_bar, ...
    tailURange_rad, tailUSlope_radPerS, ...
    measuredU_rad, measuredV_rad, measuredPhi_rad, modelU_rad, ...
    modelV_rad, modelPhi_rad, ...
    measuredLength_mm(:,1), measuredLength_mm(:,2), measuredLength_mm(:,3), ...
    modelLength_mm(:,1), modelLength_mm(:,2), modelLength_mm(:,3), ...
    position_mm(:,1), position_mm(:,2), position_mm(:,3), ...
    'VariableNames', {'phaseID','kind','start_s','end_s','validFraction', ...
    'nTail','common_bar','delta2_bar','delta3_bar','regulatorError_bar', ...
    'tailURange_rad','tailUSlope_radPerS', ...
    'measuredU_rad','measuredV_rad','measuredPhi_rad','modelU_rad', ...
    'modelV_rad','modelPhi_rad','l1_meas_mm','l2_meas_mm','l3_meas_mm', ...
    'l1_model_mm','l2_model_mm','l3_model_mm','x_mm','y_mm','z_mm'});
end

function tableOut = summarizeConditions(phases, radius_m, forcePerBar_N)
eligible = phases.common_bar > 0.5 & phases.kind ~= "mixed";
levels = unique(round(phases.common_bar(eligible)*10)/10);
n = numel(levels);
common_bar = levels(:); nNeutral = zeros(n,1);
neutralU_rad = NaN(n,1); positiveU_rad = NaN(n,1);
negativeU_rad = NaN(n,1); neutralDrift_rad = NaN(n,1);
neutralV_rad = NaN(n,1); positiveV_rad = NaN(n,1);
negativeV_rad = NaN(n,1); measuredVGain_radPerBar = NaN(n,1);
returnAfterPositive_rad = NaN(n,1); returnAfterNegative_rad = NaN(n,1);
measuredGain_radPerBar = NaN(n,1); modelGain_radPerBar = NaN(n,1);
apparentScalarK_NperM = NaN(n,1); modelScalarK_NperM = NaN(n,1);
neutralPhi_rad = NaN(n,1); positivePhi_rad = NaN(n,1);
negativePhi_rad = NaN(n,1);

for j = 1:n
    atLevel = abs(phases.common_bar - levels(j)) < 0.05;
    neutral = find(atLevel & phases.kind == "neutral");
    positive = find(atLevel & phases.kind == "positive");
    negative = find(atLevel & phases.kind == "negative");
    nNeutral(j) = numel(neutral);
    if ~isempty(neutral)
        neutralU_rad(j) = mean(phases.measuredU_rad(neutral), 'omitnan');
        neutralV_rad(j) = mean(phases.measuredV_rad(neutral), 'omitnan');
        neutralPhi_rad(j) = mean(phases.measuredPhi_rad(neutral), 'omitnan');
        neutralDrift_rad(j) = phases.measuredU_rad(neutral(end)) ...
            - phases.measuredU_rad(neutral(1));
        if numel(neutral) >= 3
            returnAfterPositive_rad(j) = phases.measuredU_rad(neutral(2)) ...
                - phases.measuredU_rad(neutral(1));
            returnAfterNegative_rad(j) = phases.measuredU_rad(neutral(3)) ...
                - phases.measuredU_rad(neutral(2));
        end
    end
    if isempty(positive) || isempty(negative)
        continue;
    end
    positiveU_rad(j) = mean(phases.measuredU_rad(positive), 'omitnan');
    negativeU_rad(j) = mean(phases.measuredU_rad(negative), 'omitnan');
    positiveV_rad(j) = mean(phases.measuredV_rad(positive), 'omitnan');
    negativeV_rad(j) = mean(phases.measuredV_rad(negative), 'omitnan');
    positivePhi_rad(j) = mean(phases.measuredPhi_rad(positive), 'omitnan');
    negativePhi_rad(j) = mean(phases.measuredPhi_rad(negative), 'omitnan');
    plusDelta = mean((phases.delta2_bar(positive) ...
        + phases.delta3_bar(positive))/2, 'omitnan');
    minusDelta = mean((phases.delta2_bar(negative) ...
        + phases.delta3_bar(negative))/2, 'omitnan');
    deltaSpan = plusDelta - minusDelta;
    if deltaSpan <= 0
        continue;
    end
    measuredGain_radPerBar(j) = (positiveU_rad(j) - negativeU_rad(j))/deltaSpan;
    measuredVGain_radPerBar(j) = (positiveV_rad(j) - negativeV_rad(j))/deltaSpan;
    modelGain_radPerBar(j) = (mean(phases.modelU_rad(positive), 'omitnan') ...
        - mean(phases.modelU_rad(negative), 'omitnan'))/deltaSpan;
    if measuredGain_radPerBar(j) > 0
        apparentScalarK_NperM(j) = 2*forcePerBar_N ...
            /(3*radius_m*measuredGain_radPerBar(j));
    end
    if modelGain_radPerBar(j) > 0
        modelScalarK_NperM(j) = 2*forcePerBar_N ...
            /(3*radius_m*modelGain_radPerBar(j));
    end
end

tableOut = table(common_bar,nNeutral,neutralU_rad,positiveU_rad,negativeU_rad, ...
    neutralV_rad,positiveV_rad,negativeV_rad,neutralDrift_rad, ...
    returnAfterPositive_rad,returnAfterNegative_rad, ...
    measuredGain_radPerBar,measuredVGain_radPerBar,modelGain_radPerBar, ...
    apparentScalarK_NperM,modelScalarK_NperM,neutralPhi_rad, ...
    positivePhi_rad,negativePhi_rad);
end

function writeReport(results, path)
d = results.data; p = results.phaseSummary; c = results.conditionSummary;
file = fopen(path, 'w');
assert(file >= 0, 'Could not write %s.', path);
cleaner = onCleanup(@()fclose(file)); %#ok<NASGU>
fprintf(file, 'Common-pressure experiment: %s\n', results.sourceMat);
fprintf(file, 'Simulation range: %.2f to %.2f s; NDI samples: %d\n', ...
    d.t_s(1), d.t_s(end), numel(d.t_s));
fprintf(file, 'Valid orientation/paired-sensor/configuration samples: %d/%d (%.1f%%)\n', ...
    nnz(d.valid), numel(d.valid), 100*mean(d.valid));
pollIntervals = diff(d.pollStart_s);
pollIntervals = pollIntervals(isfinite(pollIntervals) & pollIntervals > 0);
fprintf(file, 'Physical poll interval median: %.4f s; observed poll span: %.2f s\n', ...
    median(pollIntervals), d.pollStart_s(end)-d.pollStart_s(1));
fprintf(file, 'Logged raw Kill dashboard state: %.0f\n', d.loggedKillState);
fprintf(file, 'Maximum regulator-path versus commanded signal difference: %.5f bar\n', ...
    max(abs(d.regulatorPath_bar-d.command_bar), [], 'all'));
fprintf(file, 'These are logged commands, not measured chamber pressures.\n');
fprintf(file, 'U=(l2+l3)/r is signed bending toward theta=0; V=(l3-l2)/(sqrt(3)*r).\n');
fprintf(file, 'Each phase row averages its final %.1f simulation seconds.\n\n', ...
    results.assumptions.tailWindow_s);
fprintf(file, 'Largest 5 s within-phase U range: %.5f rad; largest absolute slope: %.5f rad/s\n\n', ...
    max(p.tailURange_rad,[],'omitnan'), max(abs(p.tailUSlope_radPerS),[],'omitnan'));
fprintf(file, 'Common  neutral-U  +delta-U  -delta-U  neutral drift  meas gain  model gain  apparent K\n');
fprintf(file, ' bar      rad        rad       rad          rad       rad/bar    rad/bar      N/m\n');
for j = 1:height(c)
    fprintf(file, '%5.1f  %+9.4f  %+8.4f  %+8.4f  %+12.4f  %+9.4f  %+10.4f  %10.1f\n', ...
        c.common_bar(j), c.neutralU_rad(j), c.positiveU_rad(j), ...
        c.negativeU_rad(j), c.neutralDrift_rad(j), ...
        c.measuredGain_radPerBar(j), c.modelGain_radPerBar(j), ...
        c.apparentScalarK_NperM(j));
end
fprintf(file, '\nCommon  neutral-V  differential V gain  return after +  return after -\n');
fprintf(file, ' bar       rad           rad/bar              rad             rad\n');
for j = 1:height(c)
    fprintf(file, '%5.1f  %+9.4f       %+10.4f        %+10.4f      %+10.4f\n', ...
        c.common_bar(j), c.neutralV_rad(j), c.measuredVGain_radPerBar(j), ...
        c.returnAfterPositive_rad(j), c.returnAfterNegative_rad(j));
end
if height(p) >= 2
    fprintf(file, '\nZero-pressure U drift, final minus initial: %+.5f rad\n', ...
        p.measuredU_rad(end)-p.measuredU_rad(1));
end
fprintf(file, '\nApparent K assumes nominal equal pressure-force gains, the\n');
fprintf(file, 'symmetric K=k*[2 1;1 2] model, equilibrium, and repeatable valves.\n');
fprintf(file, 'It is a command-to-motion diagnostic, not identified material stiffness.\n');
fprintf(file, 'The forward model predicts zero bend at equal commands; measured\n');
fprintf(file, 'equal-command bend can reflect asymmetry, pressure dynamics, or preload.\n');
fprintf(file, 'Even the model-derived apparent K can differ from its configured\n');
fprintf(file, 'K=650 because these are finite-time phase-tail responses.\n');
fprintf(file, 'Cross-plane V and within-hold motion further limit that interpretation.\n');
fprintf(file, 'Positive and negative steps were run in one order at each level;\n');
fprintf(file, 'time drift and hysteresis are not separately identifiable here.\n');
end

function makePlots(results, outputDir)
d = results.data; p = results.phaseSummary; c = results.conditionSummary;
t = d.t_s;
boundaries = p.start_s(ismember(p.phaseID, [2 7 12 17]));
colors = lines(3);

f = figure('Visible','off','Position',[80 80 1300 850]);
tiles = tiledlayout(f,3,1,'TileSpacing','compact');
for j = 1:3
    ax = nexttile(tiles);
    plot(ax,t,d.command_bar(:,j),'Color',colors(j,:),'LineWidth',1.2); hold(ax,'on');
    plot(ax,t,d.regulatorPath_bar(:,j),'k--','LineWidth',0.8);
    yline(ax,0.8,':','Nominal deadzone');
    addBoundaries(ax,boundaries);
    ylabel(ax,sprintf('P%d (bar)',j)); ylim(ax,[0 3.1]); grid(ax,'on');
    if j == 1
        legend(ax,{'Command','Regulator-path command','0.8 bar'},'Location','northwest');
    end
end
xlabel(nexttile(tiles,3),'Simulation time (s)');
title(tiles,'Common-pressure commands and valve-path signals');
exportgraphics(f,fullfile(outputDir,'pressure_trajectories.png'),'Resolution',160); close(f);

f = figure('Visible','off','Position',[80 80 1300 850]);
tiles = tiledlayout(f,3,1,'TileSpacing','compact');
labels = {'l_1','l_2','l_3'};
for j = 1:3
    ax = nexttile(tiles);
    plot(ax,t,1e3*d.measuredLength_m(:,j),'b','LineWidth',1.0); hold(ax,'on');
    plot(ax,t,1e3*d.modelLength_m(:,j),'r--','LineWidth',1.0);
    addBoundaries(ax,boundaries);
    ylabel(ax,sprintf('%s (mm)',labels{j})); grid(ax,'on');
    if j == 1
        legend(ax,{'NDI orientation-derived','Forward model'},'Location','best');
    end
end
xlabel(nexttile(tiles,3),'Simulation time (s)');
title(tiles,'Actuator length changes');
exportgraphics(f,fullfile(outputDir,'lengths.png'),'Resolution',160); close(f);

f = figure('Visible','off','Position',[80 80 1300 900]);
tiles = tiledlayout(f,3,1,'TileSpacing','compact');
measured = {d.measuredU_rad,d.measuredV_rad,d.measuredPhi_rad};
modeled = {d.modelU_rad,d.modelV_rad,d.modelPhi_rad};
axisLabels = {'Signed U (rad)','Cross-plane V (rad)','Bend magnitude phi (rad)'};
for j = 1:3
    ax = nexttile(tiles);
    plot(ax,t,measured{j},'b','LineWidth',1.0); hold(ax,'on');
    plot(ax,t,modeled{j},'r--','LineWidth',1.0);
    addBoundaries(ax,boundaries);
    ylabel(ax,axisLabels{j}); grid(ax,'on');
    if j == 1
        legend(ax,{'NDI','Forward model'},'Location','best');
    end
end
xlabel(nexttile(tiles,3),'Simulation time (s)');
title(tiles,'Signed bending separates reversal from bend magnitude');
exportgraphics(f,fullfile(outputDir,'bending.png'),'Resolution',160); close(f);

f = figure('Visible','off','Position',[80 80 1500 850]);
tiles = tiledlayout(f,2,3,'TileSpacing','compact');
ax = nexttile(tiles);
plot(ax,c.common_bar,c.neutralU_rad,'ko-','LineWidth',1.3); hold(ax,'on');
plot(ax,c.common_bar,c.positiveU_rad,'bo-','LineWidth',1.3);
plot(ax,c.common_bar,c.negativeU_rad,'ro-','LineWidth',1.3);
ylabel(ax,'Settled measured U (rad)'); grid(ax,'on');
legend(ax,{'Equal pressure','Positive differential','Negative differential'},'Location','best');
ax = nexttile(tiles);
plot(ax,c.common_bar,c.measuredGain_radPerBar,'bo-','LineWidth',1.3); hold(ax,'on');
plot(ax,c.common_bar,c.modelGain_radPerBar,'rs--','LineWidth',1.3);
ylabel(ax,'Signed gain U / differential (rad/bar)'); grid(ax,'on');
legend(ax,{'NDI','Forward model'},'Location','best');
ax = nexttile(tiles);
plot(ax,c.common_bar,c.neutralV_rad,'ko-','LineWidth',1.3);
ylabel(ax,'Equal-pressure V (rad)'); grid(ax,'on');
ax = nexttile(tiles);
plot(ax,c.common_bar,c.measuredVGain_radPerBar,'mo-','LineWidth',1.3);
ylabel(ax,'Differential V gain (rad/bar)'); grid(ax,'on');
ax = nexttile(tiles);
plot(ax,c.common_bar,c.returnAfterPositive_rad,'bo-','LineWidth',1.3); hold(ax,'on');
plot(ax,c.common_bar,c.returnAfterNegative_rad,'ro-','LineWidth',1.3);
ylabel(ax,'Return-to-neutral U shift (rad)'); grid(ax,'on');
legend(ax,{'After + step','After - step'},'Location','best');
ax = nexttile(tiles);
plot(ax,c.common_bar,c.apparentScalarK_NperM,'bo-','LineWidth',1.3); hold(ax,'on');
plot(ax,c.common_bar,c.modelScalarK_NperM,'rs--','LineWidth',1.3);
ylabel(ax,'Conditional apparent scalar K (N/m)'); grid(ax,'on');
for ax = findall(f,'Type','Axes')'
    xlabel(ax,'Common command (bar)');
end
title(tiles,'Settled common-pressure comparison');
exportgraphics(f,fullfile(outputDir,'condition_comparison.png'),'Resolution',160); close(f);

baseline = dataBaseline(d);
f = figure('Visible','off','Position',[80 80 1300 850]);
tiles = tiledlayout(f,3,1,'TileSpacing','compact');
for j = 1:3
    ax = nexttile(tiles);
    plot(ax,t,1e3*(d.position_m(:,j)-baseline(j)),'Color',colors(j,:),'LineWidth',1);
    addBoundaries(ax,boundaries);
    axisNames = 'XYZ';
    ylabel(ax,sprintf('Delta %c (mm)',axisNames(j))); grid(ax,'on');
end
xlabel(nexttile(tiles,3),'Simulation time (s)');
title(tiles,'NDI tip-position change from the initial zero-pressure hold');
exportgraphics(f,fullfile(outputDir,'tip_position.png'),'Resolution',160); close(f);
end

function baseline = dataBaseline(data)
mask = data.phaseID == 1 & data.valid & data.t_s >= 5;
assert(any(mask), 'No valid zero-pressure baseline samples.');
baseline = mean(data.position_m(mask,:),1,'omitnan');
end

function addBoundaries(ax, boundaries)
for boundary = boundaries(:)'
    xline(ax,boundary,':','Color',[0.65 0.65 0.65]);
end
end
