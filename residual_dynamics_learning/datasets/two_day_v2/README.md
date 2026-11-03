# Prepared two-day residual dataset: step2

Prepared2026-11-01 offline. This is the valid preparation version. DO NOT USE
the incomplete two_day_v1 attempt, whose nominal-helper bound penalty was
transcribed incorrectly. V2 helper is checked directly against existing RHS
inside/near/outside bounds, not merely against its own inverse.

## Data separation and counts

Ten raw files remain unchanged. Per-source SHA256 is verified against step1
provenance before export. Every output directory includes cohort and run name.

| Partition | Cohort | Eligible anchors |
|---|---|---:|
| Train | October19 | 35,350 |
| Validation | October19, held-out whole trials | 13,977 |
| External test | October27 | 24,079 |

Total73,406 shared eligible anchors, with distinct forward/inverse tasks.
Overlapping windows are NOT independent observations. Development split holds
out the middle listed trial per channel/run stratum, with at least one trial
remaining for training. For two repeats this holds one of two out; sine has
roughly50%validation. Report per-run/trial performance, not just pooled rows.
No normalization or model fitting uses the external day. Test data is available
on disk but designated locked for later model evaluation, not hyperparameter fit.

## Artifacts

- configuration.mat and preparation.mat: exact options, source provenance and
  summary. summary.csv contains per-run window counts and partition counts.
- Each run/export: observations.csv, native-rate commands.csv, model predictions
  kept separately, signal dictionary/inventory and metadata from existing
  exporter. Raw MAT remains authoritative; writeNative=false avoids redundant
  per-signal CSVs. Constant input-tail semantics are declared by the exporter.
- Each run/sequences.mat: q,p,physical/simulation times, masks, split/trial IDs,
  gap segments, causal derived v/a, fixed nominal quantities, candidate residual
  labels and eligible anchor indices. No large duplicated tensors are needed.
- row_manifest.csv/window_index.csv: trace every anchor to original sensor row.
- normalizer.mat: q,dq,p means/std from mature TRAIN rows only; inverse target
  scales from TRAIN inverse anchors only. No parameters fit on validation/test.
- quality.csv and recorded_map_verification.csv: timing/domain and historical
  model-convention checks. source_snapshot contains nominal/preparation sources.

## Timing, masks and windows

Native NDI simulation timestamps anchor rows. Commands use previous-sample
alignment, max log age.075s. Original command streams (including100Hz) remain in
commands.csv; do not assume no command changes between sensor samples.
Physical dt uses recorded poll-start host time, not simulation dt. Actual NI
write wall-clock times are absent: command-to-sensor physical alignment is an
explicit approximate assumption, not a measured actuator delay.

Require fresh paired orientation, base validity, finite inferred q/commands,
pressure proxies within0..3bar, positive trial ID and bounded command age.
Theta observability is not required for valid straight q. Invalid/gap rows are
masked, never filled. Segment at every invalid row, trial/channel/split boundary,
physical dt<=0 or>.10s, or simulation interval>.075s. No window crosses a segment.

Recompute causal v (tau.20s) and a (tau.30s) from q and physical time, reset at
each segment, exclude first.60s of that segment. These are derived quantities,
not independent velocity/acceleration sensors. Filter constants are initial
preprocessing choices, not optimized arm/controller settings.

History10frames, inverse preview5frames: nominally.5/.25s in SIMULATION time.
Physical durations vary: retained dt features are necessary. Both tasks share
the stricter eligibility interval k-10..k+5 to support matched benchmarking.
All rows in that interval must be mature and in one segment/partition. These
initial history lengths require later development-only selection.

## Exact model input/target contract

Use residualBatch(S,anchors,task,cfg,normalizer). Xsequence is8x10xexamples;
Xflat isexamplesx80. Columns per sequence instant: normalized q2,q3,dq2,dq3,
p1,p2,p3, then physical dt seconds unscaled.

- Forward: q/v history ends at k; pressure history also ends at k. Y is observed
  q(k+1), with dtNext. No future response enters X; preview is NaN. Integrating
  multi-rate inputs for predictive comparisons must use retained command streams
  or explicitly declare a coarse zero-order-hold approximation.
- Inverse: q/v history ends at k but pressure history ends at k-1. Preview is
  observed q(k+1:k+5) defining an OFFLINE motion requirement; deployment must
  replace this by planned desired motion, never future measurements. p0 is a
  separate externally prescribed common-pressure condition. In data p0=mean(p_k)
  defines this conditional inverse task, not three unknown independent outputs.
  Y is delta_d in bar, where d=[p2-p1,p3-p1] and delta_d=d_applied-d_nominal.
  Inverse target mean/std are supplied separately; batch Y remains raw units.

Nominal d is computed from current observed q and causal derived v/a using the
fixed collection mechanical equations and alpha=pi*.0065^2*1e5. No pressure
clipping is applied to nominal inverse labels. A nominal pressure above3bar is
flagged, not excluded: it can reflect model mismatch the residual must learn.
Actual command labels remain bounded. Candidate acceleration residuals are
also retained but are NOT measured acceleration/force ground truth. Prefer
output-error/rollout evaluation before trusting derivative-based learning.

## Important findings and limitations

October27 step/combined runs have median physical dt~.0885/.0877s with2001/1194
gaps>.10s. Only46/26anchors remain after complete-window and warmup rules.
These subsets are small and selection-biased toward short uninterrupted periods;
they cannot support broad step/combined generalization claims. Sine retains2952,
circle11339, helix9716 external anchors. No thresholds were loosened after seeing
these counts. A later variable-time strategy or new acquisition may be needed.

All logged nominal effective-pressure force maps match alpha*A. Historical
inverse convention differs: October19 step/sine/combined match.8bar activation;
all others match0bar, within1e-6N numerical tolerance. This is an inverse-model
convention finding, NOT proof the physical actuator deadzone changed. Candidate
labels are recomputed against one explicit ZERO-activation reference model,
not copied from heterogeneous logged inverse predictions. Mechanical RHS
consistency does not establish that every historical nominal model is identical.

Calibration/mounting/load comparability across days remains unverified.
No neural model trained, no predictive test-day metrics used for selection,
no Simulink model run/edited and no hardware command sent. Source MATs preserved.

Validation: validateResidualNominal, validateResidualData(folder),
verifyResidualRecordedMaps(folder). These check original RHS agreement,
whole-trial/day/window guards, train-only scaling, pressure/future-target leakage
and the recorded force-map conventions. They do not establish trained accuracy.
