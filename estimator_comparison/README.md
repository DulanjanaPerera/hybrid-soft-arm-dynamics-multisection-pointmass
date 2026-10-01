# Offline NDI velocity-estimator comparison

Run this from the project root in MATLAB R2025a:

```matlab
addpath(fullfile(pwd, 'estimator_comparison'));
[summary, results] = compareNdiEstimators();
```

The function reads the six `forwardModelValidation/*sineWave.mat` recordings.
It only loads saved data; it does not open COM12, run Simulink, or command NI
outputs. Results go into this folder's `results/` directory: one plot per
trial, `estimator_summary.csv`, and `estimator_comparison.mat` containing the
summary, numeric traces, and exact settings.

The comparison uses the NDI reader's orientation-derived length changes
`[l2,l3]`, not independently measured muscle strains. The source logging name
is currently misspelled `lengthChnage_m`; the script also accepts the corrected
`lengthChange_m`. Invalid orientation frames do not update either estimator.

Two causal estimators are compared:

1. Timestamp-based backward difference followed by a first-order low-pass
   filter. The current default filter time constant is 0.20 s.
2. A four-state, constant-velocity Kalman filter with state
   `[l2,l3,dl2,dl3]` and measurement `[l2,l3]`. Initial measurement standard
   deviation is 0.05 mm; per-step acceleration standard deviation is
   0.003 m/s^2. It predicts through invalid frames and corrects only on valid
   frames. Both estimators reset velocity after a gap longer than 0.5 s.

The default timestamp is `pollStart_s`. The saved `hostTime` occurs after the
BX reply and includes variable serial latency. `pollStart_s` is a practical
timing proxy, not a synchronized NDI frame timestamp. To compare using reply
completion time instead, run:

```matlab
cfg = struct('TimestampSignal', 'hostTime', ...
    'OutputDir', fullfile(pwd, 'estimator_comparison', 'results_hostTime'));
compareNdiEstimators(cfg);
```

The CSV reports two *proxies*, not velocity-estimation accuracy: velocity RMS
in the first three simulated seconds while commanded pressure is at most
0.3 bar, and one-step prediction RMS against the next NDI length sample. The
arm is not guaranteed motionless in the low-command window, and the next NDI
sample has measurement noise. Review the plots and repeat with an independent
motion reference before claiming absolute velocity accuracy.

Tuning values are optional `cfg` fields. For example:

```matlab
cfg = struct('FilterTau_s', 0.15, ...
    'MeasurementStd_m', 1e-4, ...
    'AccelerationStd_mps2', 0.003);
compareNdiEstimators(cfg);
```

The current findings and data provenance are in `../PROJECT_LOG.md`. The
script and figures are analysis artifacts; neither estimator is yet connected
to feedback linearization.
