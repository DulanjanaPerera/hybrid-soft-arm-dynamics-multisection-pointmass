# PAC candidate1/5/8/9 presentation comparison

Created2026-11-06 from four saved no-payload, requestedphi1rad,10s-circle runs.
No hardware, Simulink model, controller, acquisition or pressure execution;
no raw data, tuning settings, filters or model files changed. No commit/push.

## Single figure

`autotuning_candidates_1_5_8_9.png` is the six-panel comparison:

- Top: measured/reference theta and phi.
- Middle: signed circular theta error and signed phi error, radians.
- Bottom: signed l2 and l3 errors, millimetres (the two mechanical coordinates).

All errors are **desired minus measured**. Theta error uses atan2(sin,cos).
Displayed theta is lifted to the desired winding at each valid sample and both
reference/measurement have the initial reference theta subtracted; no fitted
phase/time shift or across-gap unwrap is applied. The top theta curves nearly
overlap at this scale; the azimuth-error panel exposes their differences.

Time0 is simulation20s. All2400samples of the common [20,140)s tracking window
(12circles) are retained; ramp-down is excluded. The reference period10s and
phi1rad were checked from logged values, and aligned reference trajectories
match across all four runs. No smoothing, downsampling or missing-value filling.

## Exact identity mapping

| Display label | Source MAT |
|---|---|
| Candidate1 (baseline) | baseline_00g_phi1_T10_12cycles_r01.mat |
| Candidate5 | candidate05_00g_phi1_T10_12cycles_r01.mat |
| Candidate8 | candidate08_00g_phi1_T10_12cycles_r01.mat |
| Candidate9 (smooth) | smooth_candidate09_00g_phi1_T10_12cycles_r01.mat |

No candidate01 MAT exists. Candidate1 is explicitly interpreted as the first
baseline observation, not baseline r02 or the short smoke run. Candidate9 belongs
to the later broader, smoothness-primary campaign. Its label is not evidence
these four were matched independent optimization repeats or that candidate9
improved the objective. The reported comparison preserves these identities.

Logged lambda=[18,16.5],Ks=[15,10] for allfour; no within-tracking gamma/lambda/Ks
changes. GammaK/D/H differ by candidate; recordings.csv contains exact values
and raw SHA256s. The physical reset follows prior operator attestation, not
independent physical-state verification.

## Validity and checks

Native logs are read through pacTuneReadTrial, with unambiguous signal names
and paths retained by the reader. Alignment is previous sample with maxage.075s,
no extrapolation. Only scalar flags exactly1 are accepted: orientationValid,
configurationValid,baseOreintationValid,modelValid,status3, plus nonnegative
measurementAge <= the recorded supervisor maximumage. Theta additionally needs
thetaObservable. Nonfinite angles/lengths are excluded from their own plots.
Held invalid measurements are masked with NaN. Validtracking rows are2380,
2380,2376,2376 in displayed candidate order. Lengths are logged inferred length
changes, not independent strain measurements.

All four full120s and late30s coordinateRMSE/angularMAE values independently
reproduce the archived campaign evaluation MATs within1e-9 (mm/rad). See
metrics.csv and verification.json. Candidate8 has the smallest late coordinate
errors here: [1.1725,1.3029]mm; baseline[1.5684,1.8887],candidate5[1.8360,2.0686],
candidate9[2.5279,2.7520]. These are separate experimental runs, not an isolated
causal test of a single gain. Roughness/acceptance is not scored by this figure.

## Formats / retained data

- PNG300dpi, vectorPDF/SVG and editableFIG of the same figure.
- Four candidate_XX_plot_data.csv files: simulation/relative time, plotted
  signals, reference/errors and validity masks; plot_data.mat retains arrays.
- recordings.csv (identities/gains/hashes), metrics.csv and verification.json.
- Source/dependency snapshots, input/source and exported-artifact SHA256 tables.

Reproduce offline from the project root, choosing a new output directory:

```matlab
addpath(fullfile(pwd,'passivity_adaptive_control'));
plotPacCandidateComparison(pwd,fullfile(pwd,'tmp','new_candidate_comparison'));
```
