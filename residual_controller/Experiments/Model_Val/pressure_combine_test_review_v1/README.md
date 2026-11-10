# Operator pressure-combine shadow-model test

Reviewed2026-11-04 from ../pressure_combine_test.mat, saved out only. Runtime
model residual_FWD_INV_validation, version1.50, normal simulation,0..191s.
ReachedStopTime with no saved simulation error. No hardware, model execution,
SLX modification, training, gain change, commit or push during this review.

## Result

The residual-augmented force/physics model predicts the next measured lengths
much better than the nominal physics model on this run. The reported values
below use **previously issued0.05s predictions**, aligned to the next accepted
NDI observation, with scoreValid and independent orientation/base validity.
These are prediction errors, not reference-tracking or physical force errors.

| Predictor | l2 RMSE, all scored (mm) | l3 RMSE, all scored (mm) | l2 RMSE, moving (mm) | l3 RMSE, moving (mm) |
|---|---:|---:|---:|---:|
| Residual + physics | 0.234 | 0.198 | 0.346 | 0.300 |
| Nominal physics | 1.007 | 1.569 | 1.700 | 2.491 |
| Constant velocity | 0.046 | 0.047 | 0.074 | 0.074 |

All-scored3560rows; moving1043rows. Motion is defined before inspecting errors
as norm(observer causal velocity)>0.001m/s. The pressure-present subset
(max final command>0.2bar),2174rows, gives residual RMSE[.296,.248]mm versus
nominal[1.284,1.994]mm and CV[.059,.060]mm. Full metrics include MAE,95th
percentile, maxima and each trial separately; quiet portions are not the only
reported subset. Residual improves each coordinate's all-row RMSE over nominal
by approximately76.8%/87.4%. Its95th absolute errors are[.521,.441]mm and maxima
[1.098,.808]mm on scored observations.

Constant velocity is clearly better at this short horizon, including motion.
It uses the current q and estimated velocity and needs no pressure model; it
does not provide an inverse pressure/force model. This baseline result means
visual agreement/low one-step error alone cannot establish learned dynamic
accuracy, longer rollout quality, accurate pressure inversion or better PAC/PDFB
control. Keep the CV result alongside the nominal/residual comparison.

## Independent checks of wiring, causality and model identity

- Diagnostic qMeasured and context measured q exactly match independent raw
  orientation-inferred NDI l2/l3. There is no accidental reference-length match.
- For every scored row, predictionAtMeasurement equals predictionLength issued
  exactly **one observer tick earlier**, with a matching logged origin time.
  No current measurement is substituted into the displayed prediction.
- Logged prediction errors exactly equal independently recalculated prediction
  minus NDI values; retimed error fields agree as well.
- Context pressures match the final P1/P2/P3_NI_DAQ commands exactly. Logged
  simulation-clock command age is0 throughout. This checks the recorded sample
  association, not chamber pressure or physical valve latency.
- Re-evaluating101distributed/max-speed logged contexts with the frozen GRU128
  gives max force difference2.66e-15N. Original SDIRK forward propagation
  reproduces the issued three-predictor matrix to1.99e-17m. Thus the logged
  residual column is consistent with the existing frozen packet and physics.
- Logged horizon0.05s, width128. The eleven-entry shadowOptions vector itself
  was not saved; do not claim every runtime threshold independently recovered.
  Save options alongside out on subsequent trials.

## Validity and timing

3821observer rows,3814orientation valid,3821base valid. Ready3580rows;
prospective scored3560 (93.17% of all rows;93.34% of physically valid rows),
retimed scored3573. Ten valid-data segments; no observer status7 numeric guard
event. Invalid orientation samples occur at102.80,103.25,106.45,116.35,
147.30,152.35 and160.45s. Each invalid period is excluded with subsequent
warmup; missing plots are intentional, not interpolated predictions.

Host poll intervals median0.050807s/P950.050905s. The two intervals over0.10s
occur at startup (simulation0.05/0.10s),0.294/0.451s; do not describe these as
mid-run instability. Host elapsed time exceeds simulation elapsed by3.361s.
Scored interval minus nominal horizon: median0.807ms/P950.899ms/max4.810ms,
within the configured scoring policy. Prospective errors therefore retain a
small timing mismatch; retimed errors recompute the previous state/pressure/
force over actual elapsed time, known afterwards. Retimed full-row residual
RMSE[.238,.201]mm is similar, but is a different measurement of performance.

Final commands span0..3bar in each channel; physical Kill signal stays1.
Maximum absolute inferred lengths[l1,l2,l3]=[8.130,6.289,11.262]mm. Neural force
peaks[10.406,10.792]N are far inside the80N force head bound. Maximum absolute
normalized feature4.118; no ready feature exceeds5. These are input-range
diagnostics, not a model confidence bound or a stability guarantee.

## What this run does and does not validate

The operator successfully compiled/executed the forward shadow model and logged
causal predictions agreeing with the frozen implementation. This supersedes
the earlier pending installation/smoke status for this operator model version.
It does not establish a worst-case real-time compute budget: NDI poll timing
also includes acquisition/pacing, and block execution time was not logged.

No residualPressurePrediction, nominalPressurePrediction or learned inverse
diagnostic bus is logged. inverseContext is the forward observer's context;
it is not an inverse-model result. The older collection inverseModelValid and
tau/pressure signals refer to a separate nominal collection inverse, so they
cannot establish neural inverse convergence or physical inverse accuracy.
The PAC cold-history recovery limitation remains unchanged by this shadow test.

Recommended next sequence: retain these frozen weights and compare0.1/0.25/0.5s
causal rollouts offline using this recording and accurately aligned pressure
history; observe the learned inverse separately with outputs disconnected from
valves. Broader controlled experiments at different speeds/directions are
needed before residual-on PAC/PDFB conclusions. Do not infer a need for neural
retraining or model-capacity changes from this single one-step test alone.

## Retained artifacts and reproduction

- tracking_full and tracking_motion_zoom: PNG/PDF/FIG. Zoom is centred on the
  maximum measured causal-speed row, not chosen to minimize prediction error.
- diagnostics: PNG/PDF/FIG for validity, scoring, time and neural forces.
- prediction_metrics.csv, prediction_by_trial.csv, observer_status.csv.
- quality_trace.csv and quality_supplement.json; summary.json and analysis.mat.
- signal_inventory.csv distinguishes bus names/duplicate signal paths.
- source_snapshot contains analysis/helpers, raw/model/packet hashes and saved
  SLX XML/source excerpts inspected read-only. Raw recordings are not rewritten.

From the project root: addpath('residual_controller');
analyzeResidualShadowLog(inputPath,freshOutputFolder). No model is opened/run.
