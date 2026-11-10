# Sine-wave length prediction: presentation figures

Created2026-11-05 from sinewave_test.mat and sinewave_test_100g.mat in the
parent folder. Both saved operator runs use residual_FWD_INV_validation1.51,
0..212s, GRU128 and logged0.05s prediction horizon. Payload labels follow the
user's description/filenames; payload attachment/calibration timing is not
embedded in these out-only recordings. No fitting, controller/SLX changes,
hardware/model execution, commit or push during figure preparation.

## Ready-to-use figures

| Figure stem | Content |
|---|---|
| sinewave_00g_full | Full no-payload run, l2 and l3 |
| sinewave_100g_full | Full100g run, l2 and l3 |
| sinewave_00g_cycle2_by_channel | No payload: one cycle for each of the three pressure channels |
| sinewave_100g_cycle2_by_channel | 100g: matching cycles/channels |
| sinewave_payload_comparison_channel2_cycle2 | No payload/100g side-by-side, matching channel2 cycle2 |

Every stem has **PNG3840x2160 at300dpi, vector PDF/SVG, editable MATLAB FIG**.
Canvas/page aspect16:9. Black=measured NDI, orange dashed=nominal,
blue=residual+nominal. Length changes are in millimetres. Shared full-run
coordinate scales across payloads; close-up scales also match corresponding
channels/coordinates across payloads. Only the three requested curves appear.

## Prediction semantics and validity

Predictions use logged predictionAtMeasurement_m: column1=residual+nominal,
column2=nominal. Both were issued at the preceding observation for a0.05s
horizon. Nominal therefore means the matched one-step physics prediction,
not the separate freely integrated nominal-simulator X trace. Measurements are
independent orientation-inferred NDI l2/l3, not desired reference lengths.

All scored predictions were checked against their previously issued output;
origin identity and independent logged-error recalculation differences are0.
Diagnostic measured lengths match raw NDI exactly; context pressure commands
match the final NI command logs exactly. The two models share the same scoring
mask. Invalid/cold-history predictions remain gaps; no interpolation or
replacement with measured values is used. Valid NDI remains visible through
observer warmup. These are forward one-step prediction comparisons, not
closed-loop tracking, learned-inverse validation or free-run rollouts.

No payload:4229/4241physically valid rows,3866scored(91.16%).
100g:4231/4241valid,3984scored(93.94%). The files do not include an independent
force/pressure measurement or a complete saved shadowOptions vector.

## Close-up selection

Both recordings contain the same six trial IDs: two repeats of three sine
cycles for each pressure channel. Actual logged phase5=periodic excitation;
phase4=zero-pressure relaxation. Selection is independent of prediction error:
cycle2 in the first trial for each channel, identically in both recordings.

| Channel | Trial | Cycle | Simulation window (s) |
|---|---:|---:|---|
| 1 | 1 | 2 | 12.90–20.70 |
| 2 | 3 | 2 | 80.00–87.80 |
| 3 | 5 | 2 | 147.15–154.95 |

Close-ups use time relative to each window start. Channel3 contains prediction
gaps in both recordings and is retained rather than selecting a better cycle.
See selected_windows.csv and trial_cycle_windows.csv for exact coverage.

## Supporting metrics

RMSE in mm, ordered[l2,l3]. Prediction metrics are computed only on common
valid scoring rows; all-run and active-sine subsets are separated.

| Payload | Subset | Nominal RMSE | Residual+nominal RMSE |
|---|---|---|---|
| No payload | All scored | [3.696,2.897] | [0.482,0.429] |
| No payload | Active sine,phase5 | [4.582,3.586] | [0.594,0.529] |
| 100g | All scored | [4.505,3.733] | [0.729,0.809] |
| 100g | Active sine,phase5 | [5.556,4.587] | [0.896,0.994] |

This export contains the requested nominal/measurement/residual comparison.
Constant velocity is omitted from the figures as requested; earlier forward
audits show why it remains an important separate short-horizon baseline. These
plots alone do not establish improved feedback control or isolate a causal
payload effect across independently calibrated experiments.

## Files and reproduction

prediction_metrics.csv: RMSE/MAE/95th absolute errors and scored counts.
recording_summary.csv: source model versions, coverage and signal checks.
sinewave_00g_plot_data.csv and sinewave_100g_plot_data.csv: times, validity,
trial metadata, measured/predicted lengths in mm and final pressure commands.
plotting_data.mat preserves the plotting arrays and selected/scaled windows.
signal_inventory.csv records source paths/duplicate signal names.
source_snapshot and ARTIFACTS.csv preserve source/raw-file/export SHA256 hashes.

Reproduce from project root with addpath('residual_controller');
plotResidualSinePresentation(pwd,freshOutputFolder). It reads saved MAT files
only and deliberately refuses to overwrite an existing output directory.
Draft v1 stopped before plotting because phase4 was initially selected;
v2 uses actual verified phase5 labels. Original MAT files remain unchanged.
