# SineWave4 prediction plots: requested predictor grouping

2026-11-03. User requested sine prediction figures, separating constant velocity
and GRU from the other previously listed predictors. Used only those existing
five frozen predictors. No retraining, new predictor, force-residual learning,
controller/Simulink/hardware operation, commit or push.

## Figures

Each figure has rows P1/P2/P3 excitation and columns l2,l3,phi. Black is measured
orientation-derived response. Lengths are in mm and phi in rad. There is no
desired phi trajectory: this is sine pressure excitation. It is not tracking
a desired sinusoidal bend. Pressure commands remain setpoint proxies.

- sine_free_running_cv_gru: constant velocity (orange) and nominal+GRU (cyan).
- sine_free_running_others: nominal (blue), nominal+ridge (green),
  kinematic+ridge (purple).
- sine_one_step_cv_gru and sine_one_step_others: same grouping for the complete
  exported sine record. Measured history is supplied again at every sample;
  these are explicitly NOT autonomous whole-record propagations. Visible gaps
  retain excluded windows; no interpolation across acquisition gaps. Spikes
  and unfavorable predictions remain visible, without smoothing/clipping.

All four are retained as editable FIG, vector PDF and300dpiPNG. Vector PDF
generation issued MATLAB's content-complexity/performance warning; exports
completed and PNGs were visually checked for labels/legends and all traces.

## Free-running window selection

SineWave4.mat from November2, with its previously verified export and prepared
physical-time sequences. Period inferred from actual logged cycleId start
spacings: median31.425s simulation (sampled31.4/31.45s intervals), rather than
assuming current pressure-source code matches the historical recording.

Choose the FIRST chronological mature, continuous one-period window for EACH
channel. Entire initial history and predicted period must remain in the same
validity segment and periodic phase. No gap guards were relaxed to obtain a
window. This makes starts/repeats different across channels:

| Excited channel | Trial | Simulation interval, s | Physical duration, s |
|---|---:|---|---:|
| P1 | 1 | 40.80-72.25 | 31.880 |
| P2 | 3 | 214.60-246.05 | 31.860 |
| P3 | 6 | 527.35-558.80 | 31.854 |

630samples/window. Different initial sine phases and repeats mean this is not
a controlled channel-to-channel comparison. Selection is independent of model
errors. selected_windows.csv and logged_periods.csv retain the exact choices.

All models start from the same10frame measured q/v history for a given window.
Thereafter q and filtered v are predicted only. Recorded pressure/time history
is known; future measured states are used only for scoring/plotting. Original
legacy-boundmode0 learned context, physical dt, velocitytau.20s and sensor-row
pressure ZOH are unchanged. There is no per-sample true-state reset. The GRU
still uses10frame independent sequence inference, not a newly stateful model.

The original locally frozen affine nominal transition is used, not full
nonlinear ODE integration. Numerical chatter/spikes in predicted traces are
not measured arm oscillations and do not independently diagnose hardware or
physical instability. This plotting request does not investigate their cause.
All15 model/channel propagations complete without the inherited |q|>.05m
numerical guard; this is not a physical bound or stability guarantee.

## Selected-window results

Path vector RMSE over inferred l2/l3, mm, relative to measured response:

| Predictor | P1 excitation | P2 excitation | P3 excitation |
|---|---:|---:|---:|
| Constant velocity | 9.604 | 16.748 | 20.685 |
| Nominal + GRU | .951 | 3.899 | 2.401 |
| Nominal | 7.018 | 10.551 | 9.433 |
| Nominal + ridge | 3.376 | 2.724 | 4.714 |
| Kinematic + ridge | 2.954 | 3.092 | 4.851 |

GRU improves markedly on constant velocity in these windows, but has visible
spikes, especially P2. Both ridge models beat GRU on P2 in this example; the
requested visual grouping does not imply every other model is always worse.
One selected window per channel does not establish an architecture winner,
controller tracking improvement or readiness for inverse pressure control.
These GRU outputs remain position corrections, not generalized-force residuals.

## Verification and artifacts

Logged-period selection and same-segment/mature-history assertions pass.
q-to-phi mapping agrees with logged measured phi within1e-10rad.
First predicted step reproduces the EXISTING frozen one-step predictions for
all15 model/channel pairs, maximum difference2.52e-10m (float/batch arithmetic).
Independent CSV recomputation verifies every length-vector path RMSE. Raw MAT,
export, prepared sequences, normalizer, weights and inference-source hashes
are unchanged. See sources.csv and verification.txt.

Full free-running prediction/measured length,phi,pressure/time traces and
scores are in CSV/MAT. Saved one-step predictions retain source-row provenance.
Plotting/prediction source snapshot is included; raw MAT remains authoritative.
No measured pressure/force or independent length truth is invented.

Reproduce from repository root into a NEW directory:

    addpath(pwd,fullfile(pwd,'residual_dynamics_learning'));
    plotResidualSinePredictors(pwd,newFolder);
