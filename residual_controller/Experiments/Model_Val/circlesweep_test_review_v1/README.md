# Circle-sweep prediction causality and measurement-reuse audit

Reviewed2026-11-05. Input ../circlesweep_test.mat; runtime model
residual_FWD_INV_validation1.51,normal simulation,0..375s,ReachedStopTime,
no saved simulation error. No SLX/controller changes, hardware/model execution,
fitting, weight updates, commit or push during this audit. Offline synthetic
measurement offsets below never modify the original MAT recording.

## Conclusion

**No future-target NDI leakage was found in the logged forward predictor.**
There IS intentional measured-state conditioning: each new50ms forecast starts
from the latest measured q and causal filtered velocity. This makes an absolute
trajectory overlay a weak test of dynamics. It is not a simulation of the whole
circle starting from one initial condition. Good-looking plots do not establish
pressure sensitivity, inverse accuracy, longer free runs or NMPC readiness.

The forecast is F(q_k,vhat_k,p_k,h_before_k), targeting the next observation.
NDI q_(k+1) is used only later for scoring the stored forecast and initializing
the next forecast. Neither the forward map nor its GRU reads the future target.
Retimed predictions use the now-known elapsed interval but still only the prior
state/pressure/force; keep them labelled separately from prospective forecasts.

## Tests and evidence

1. **Every scored prediction was issued one tick earlier.** All6767scored
   predictionAtMeasurement values exactly equal a predictionLength value from
   the preceding logged observer tick. Independent target/origin host times
   agree. The target measurement is not substituted into the plotted signal.
2. **Independent input checks.** Observer/context q equals independent raw NDI
   lengths exactly. Context pressures equal final NI commands exactly. Logged
   errors equal independently recalculated prediction minus NDI exactly.
3. **Frozen model reproduction.**101logged immutable contexts reproduce GRU
   forces to5.33e-15N and the three-column forward output to6.61e-18m.
4. **Full sequential replay from initialization.** Reconstructing the observer
   from all7501raw measurement/frame/clock/pressure rows reproduces issued and
   matched predictions to1.39e-16m and neural force to1.07e-14N. There are zero
   status, readiness or finite-prediction-mask mismatches. This independently
   reconstructs velocity and GRU history from past/current measurements only.
   Source defaults plus logged width/horizon were used; the exact full options
   vector was not saved, so compatibility does not uniquely identify every option.
5. **Future-target falsification test.** Starting from the identical observed
   history after87.05s, replace only the next NDI length at87.10s with
   q+[5,-3]mm. The already-issued prediction changes by exactly0. Its error
   changes by[-5,+3]mm as required. The new forecast issued after receiving
   that altered measurement changes, as legitimate causal state conditioning
   requires. This test is artificial/offline, not an additional physical run.
6. **Weights unchanged.** The packed128-unit model SHA256 matches the frozen
   Nov3packet:
   a5efd346d2d23a0dfb83dcda9647f9ec5f2fe690364cb226c03fafb0363fcd60.
   Runtime contexts reproduce that packet. The new Nov5recording has not been
   used to update these weights; online recurrent-state updates are not training.

## Prediction accuracy and why the plot looks close

RMSE(mm),same6767prospective scored targets. Moving means causal velocity norm
greater than0.001m/s,2042targets; masks are shared across predictors.

| Predictor | All l2 | All l3 | Moving l2 | Moving l3 |
|---|---:|---:|---:|---:|
| Residual + nominal | 0.748 | 1.220 | 1.337 | 2.204 |
| Nominal | 4.296 | 4.714 | 6.321 | 7.328 |
| Constant velocity | 0.097 | 0.111 | 0.174 | 0.199 |
| Repeat latest measured length | 0.206 | 0.256 | 0.375 | 0.464 |

Even persistence (qhat_next=q_current) beats the residual predictor at50ms on
this run. Thus a large part of the visual closeness comes from starting close
to the target at every tick. The residual is better than poor nominal physics,
but predicts local increments substantially less accurately than the two simple
kinematic baselines. Its all-row95th absolute errors are[2.001,3.887]mm and
maxima[4.048,6.074]mm. These errors are visible in close-ups but can look small
against20–30mm absolute length changes on a long time axis.

See causality_audit/absolute_vs_increments.png/.pdf/.fig. The upper panels plot
absolute lengths; the lower panels subtract the SAME measured origin from each
prediction and target. This displays the actual predicted movement over one
step. The displayed window is centred on maximum recorded causal speed
(0.02669m/s at311.70s), not chosen for low model error. Gaps remain unfilled.

## No-reset conditional rollout check

320anchors, every20th eligible row from positive trial IDs, including88moving
origins. Each anchor starts from one measured q/v and observed GRU prefix.
Future states and GRU memory then propagate without NDI q/v corrections. Inputs
are the recorded future pressure commands and actual elapsed intervals. Future
NDI lengths are used only after propagation to score outputs. This is an
offline conditional rollout test, not a prospective prediction of unknown future
commands/timing. Mean horizon labels are nominal; actual medians are about
.0508/.1014/.2533/.5066s. All models use the same anchors; no.05m guard failures,
no clipping/interpolation or error-based anchor selection.

Moving-origin RMSE(mm),88matched anchors at every horizon:

| Nominal horizon | Residual l2/l3 | Nominal l2/l3 | CV l2/l3 | Persistence l2/l3 |
|---|---|---|---|---|
| 0.05s | [1.401,2.264] | [6.482,7.558] | [.173,.200] | [.381,.419] |
| 0.10s | [1.746,2.688] | [9.357,10.265] | [.340,.388] | [.772,.846] |
| 0.25s | [1.431,2.630] | [9.405,10.098] | [.767,.894] | [1.732,2.075] |
| 0.50s | [1.660,2.906] | [9.102,9.902] | [1.748,2.042] | [3.338,4.221] |

The residual remains better than nominal but does not consistently beat CV;
at.5s it narrowly beats CV for l2 and loses for l3. These results support neither
an independent highly accurate plant replica nor immediate NMPC qualification.
These anchor-specific50ms scores differ from the whole-record table because
the anchor set is sparse and requires.5s of valid future scoring data.

## Coverage and feature diagnostics

7479/7501orientation-valid,base allvalid. Ready6806; prospective6767(90.21% of
all rows),retimed6786.24segments;22invalid orientation rows; no numerical-guard
status7. Host-poll median.0508s and375s simulated duration; host elapsed is
5.444s longer. Scored horizon mismatch median.818ms,max3.827ms approximately
(exact values in summary.json). Final commands stay within0..3bar (observed
peak2.946bar),Kill1. Payload condition
was not supplied for this file and is not inferred from the filename.

Velocity features reach18.05/17.13stored training standard deviations from
their means;12.84% of ready observations have at least one feature above5 in
absolute standardized units. q/p/time features remain below5 in logged ready
contexts. This is a normalization diagnostic, not a calibrated uncertainty
bound or proof of exceeding the training minimum/maximum. Review velocity
coverage/state estimation and pressure-conditioned prediction accuracy before
designing an NMPC around this model; low visual overlay error is insufficient.

## Saved Scope wiring caveat

Read-only saved SLX inspection finds Scope input1 from Selector1, fed by
predictionAtMeasurement column1; input2 from qMeasured bus. That residual/NDI
comparison is consistent with the audited signals. Input3 comes from MATLAB
Function2, a separate length-derivative/model chain, **not column2 of the matched
prediction matrix**. Keep this third curve's identity explicit; comparing a
separate nominal free simulation with a measurement-reset residual predictor
is a different experiment from comparing the matched one-step columns. The
saved model may omit unsaved operator changes; audit conclusions above rely
on the logged outputs/full replay, not on that disk snapshot alone.

## Retained artifacts

Main folder: prediction_metrics/by_trial CSV, observer status/quality/inventory,
summary.json/analysis.mat and full/zoom/diagnostic PNG/PDF/FIG.
causality_audit: causality_summary.json/audit.mat, raw full replay and artificial
future-target test, conditional_rollout_metrics.csv and plot, feature_range.csv,
one_step_increments.csv and absolute_vs_increments plots. Source snapshots and
SHA inventories preserve analyzer/helpers, raw input, frozen packet and saved
SLX excerpts. Raw files and all controller/model code remain unchanged.

Reproduce with analyzeResidualShadowLog(inputPath,freshFolder), then
auditResidualShadowCausality(root,freshFolder) and
plotShadowPredictionIncrements(freshFolder). These functions do not run a model.
