# Frozen residual predictors on fresh November 2 recordings

Completed 2026-11-02, offline. Canonical version 2. No training, gain changes,
controller integration, Simulink execution, device access, commit or push.

## Result and decision

History-based prediction transfers usefully to these new recordings over
multi-second windows. The originally selected ridge improves the five-run mean
5 s path error from constant velocity's 3.408 mm to 1.620 mm (52.4%). However,
one-step improvement is small, the kinematic ridge ablation is competitive,
and no consistent benefit of the mechanical prior is established. These are
position-transition corrections, not the proposed generalized-force residual.
Do not inject these outputs into the force balance or proceed to PAC/PDFB
integration on this evidence alone.

The next implementation should be the explicitly versioned force-residual
model fitted through integrated length output error, with history-dependent
inputs and a declared corrected nominal baseline. Preserve these frozen
predictors as benchmarks. This cohort has now been evaluated: it is not an
untouched confirmation set for subsequent model/architecture changes.

## Frozen protocol

Five sources, exactly: StepWave, SineWave4, combinePressure, circleSweep, helix
under RL_data_collection/20261102. Earlier sine attempts are not pooled.
Raw acquisition timestamps confirm November 2. All five supplied observations,
commands and model-prediction CSVs exactly match fresh exports from their raw
MAT files; see raw_export_checks.csv and verification.txt. Export verification
does not establish unrecorded calibration, load or mounting comparability.

Models: screen_20261101_v2 forward/inverse banks and shortcut ablations.
Weights, selected model, training-only normalizer and source hashes were saved
before evaluating these recordings and verified unchanged afterward. No
normalizer or hyperparameter was fit to November 2. Every recording is fresh_test.
No new winner is selected from these results; historical forward selection
remains ridge. GRU's favorable 5 s aggregate is descriptive.

Preprocessing matches two_day_v2: inferred q=[l2;l3], physical poll-start dt,
causal velocity tau .20 s, acceleration tau .30 s, .60 s reset warmup,
10-frame history, 5-frame inverse preview; require the complete stricter shared
window k-10:k+5 mature and within one segment. Segment at invalid samples,
trial/channel boundaries, physical gap >.10 s/nonpositive dt or simulation gap
>.075 s. Command log age <=.075 s, commands 0..3 bar, base/orientation/freshness
requirements unchanged. No dropout interpolation or post-result exclusions.
Logged NDI diagnostic fields are retained in exports but do not relax masks or
enter model inputs. Sensor lengths/velocities are inferred/derived, not direct
strain/velocity measurements.

The original learned physics context explicitly uses coordinate-bounds-v1,
mode0. It is NOT silently replaced with the corrected all-actuator mode1.
Mode1 nominal is separately scored without learned corrections. The two nominal
results are nearly identical on the measured states, which remain mostly below
soft-bound activation. This does not invalidate the structural bound correction
or establish migration of every embedded controller. Recorded SLX versions
are 1.46 for step/combined and 1.47 for sine/circle/helix; do not imply identical
saved model code from a shared filename. Raw MATs remain authoritative.

Pressure inputs are NI setpoint proxies, not measured chamber pressures.
Predictions hold the NDI-row command over each physical interval, as in the
original screen. Actual NI-write host times are absent; this is an approximate
input/physical-clock alignment. Native 100 Hz commands are retained, but this
benchmark does not integrate their within-interval changes or identify lag.

## Coverage

| Recording | Sensor rows | Orientation-invalid rows | Guarded anchors |
|---|---:|---:|---:|
| StepWave | 11,001 | 12 | 10,016 |
| SineWave4 | 12,721 | 10 | 12,035 |
| combinePressure | 6,841 | 13 | 6,120 |
| circleSweep | 12,501 | 10 | 11,903 |
| helix | 10,601 | 11 | 10,004 |
| Total | 53,665 | 56 | 50,078 |

Poll-start median intervals are .05080-.05082 s. Additional gaps are segmented;
poll-start and completion gap counts can differ. The raw diagnostics exist for
sine/circle/helix, not step/combined. Earlier channel-3 MISSING diagnosis is not
automatically assigned to losses without corresponding raw diagnostics.

Training-normalized feature excursions remain in the test, not censored:
StepWave max absolute training z=14.37; 3.60% of mature rows have any feature
above |z|=5. Sine .366%, circle .0248%, helix .0784%, combined zero. This is a
scale-excursion diagnostic, not a formal training-support detector or exclusion
criterion; no test statistics are used to renormalize inputs.

## Forward prediction

Errors below are vector length RMSE in mm over the two independent coordinates
l2/l3, averaged equally across the five experiments. They are model-prediction
errors, not arm tracking or orientation errors. Pooled/per-coordinate/per-trial
and per-phase scores are retained separately. Overlapping anchors are not
independent experiments; these aggregates have no confidence/stability guarantee.

| Frozen predictor | One sensor step | 0.5 s path | 2 s path | 5 s path |
|---|---:|---:|---:|---:|
| Legacy nominal, locally frozen coefficients | 4.9183 | 8.1766 | 7.7494 | 7.0634 |
| Corrected-bound nominal, same method | 4.9183 | 8.1766 | 7.7493 | 7.0633 |
| Persistence | .10546 | .47068 | 1.5395 | 2.8965 |
| Constant velocity | .06235 | .18575 | .89300 | 3.4076 |
| Nominal + ridge residual | .05976 | .18803 | .72680 | 1.6204 |
| Nominal + MLP residual | .41796 | 2.4987 | 5.6626 | 6.6624 |
| Nominal + GRU residual | .34447 | 1.0735 | 1.1265 | 1.5112 |
| Kinematic + ridge residual | .05961 | .18679 | .67946 | 1.8064 |

One-step ridge is about 4.2% better than constant velocity on the macro score,
but slightly worse on sine (.09518 vs .09425 mm). It improves other four
one-step experiments. Phase scores show step loading .0816 vs .0941 mm and
unloading .0930 vs .1048 mm; on static/relaxation intervals constant velocity
can be better. Sine periodic phase remains slightly worse (.1000 vs .0990 mm).
Circle periodic and helix loading/unloading show modest ridge improvement.
Thus the conclusion does not rest only on plentiful stationary samples.

Multi-step protocol: 12 deterministic chronologically spread eligible initial
windows per run/horizon, 60 per horizon/model, target physical horizons .5/2/5 s.
Each ends at the first sensor sample reaching the target duration, with no gap
crossing; actual durations and selected row IDs are in rollout_windows.csv.
All predictors use matched starts within a horizon. Eligible start pools can
differ between horizons, so the curve is not an identical-start horizon study.
Only initial q/v history is observed. Later q comes from prediction, v from
the same causal difference/filter update applied to predicted q. Pressure
history and future query times are known recorded exogenous inputs; future
measured motion is used only to score, never to reset prediction. GRU uses its
original independent 10-frame sequence inference, not newly stateful training.

All 1,440 model/window propagations completed without exceeding the inherited
|q2| or |q3| > .05 m numerical domain guard. This is not the physical actuator
bound or a safety/stability guarantee. Success-only metric columns retain
explicit failure counts. This nominal baseline uses locally frozen affine
transitions with reconstructed velocity, not a new full nonlinear ode15s test.
Do not describe its error as an exact nonlinear-model score.

At 5 s, ridge beats constant velocity in each of the five run-level scores;
GRU has a lower macro score but is worse than constant velocity on combined
pressure (1.539 vs 1.397 mm). Ridge versus kinematic ridge changes with horizon:
the kinematic alternative is slightly better one-step/2 s, worse at 5 s.
Physics-prior value therefore remains mixed rather than demonstrated generally.
These are short prediction windows, not full-circle autonomous validation.

## Conditional inverse audit

Five-run mean three-pressure vector RMSE (bar): nominal .49851; previous
differential .009336; full ridge .000559; MLP .033216; GRU .030048;
pressure-history-only AR .000555; motion-only ridge .32967.

The pressure-only AR still matches the full inverse model. This confirms the
known smooth-command prediction shortcut; tiny pressure RMSE does not establish
an inverse model for novel desired motion. Observed future motion is explicitly
an offline requirement, with actual mean pressure supplying the conditional
common-pressure input. It is not a deployable measured-future policy.
No outputs clipped: full ridge max violation .00358 bar, none above .01 bar;
MLP/GRU max violations .165/.151 bar and 17.95%/24.84% above .01 bar.
No learned pressures were sent to the arm.

## Verification and retained artifacts

- Frozen protocol and before/after SHA256 checks, raw/export provenance and
  saved source snapshot; no fitting or source/model mutation during evaluation.
- Same-segment maturity assertions, forward future-target perturbation checks,
  inverse current-pressure-history exclusion, and target-free frozen inference.
- Exact raw re-export equality for all five observations/commands/predictions.
- sequences MAT files, anchor/row/window manifests, per-trial/per-phase scores,
  one-step predictions and initial-window rollout traces.
- fresh_prediction_summary and fresh_rollout_examples: editable FIG, vector PDF,
  300 dpi PNG. Source validateFreshResidualModels.m, copied in source_snapshot.
  Examples show the first eligible 5 s window per recording, selected by time,
  not error; many are near straight and deliberately retain unfavorable drift.
  The summary plot is the broader evidence. Full q2/q3 traces are retained.
- Additional assessment source assessFreshResidualValidation.m and transcript.
  Version1 is an incomplete evaluator attempt marked DO_NOT_USE; do not use it.

Reproduction requires a new output folder and the exact frozen files recorded
in frozen_inputs.csv. Tests and plots are analysis artifacts, not qualification
for unattended hardware or publication claims of improved PAC/PDFB control.
