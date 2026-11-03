# Fine-tuned generalized-force GRU64: completed offline comparison

2026-11-03. User authorized fine-tuning for eventual PAC and PID feedback
model augmentation. Two weight continuations completed, with a frozen
reference. No controller or hardware changes. November2 was read by this
experiment only after all weights and checkpoint choices were frozen.

**Selected research checkpoint: Old-data continuation**, using the predeclared
equal old/Nov3 internal score. File: release/force_model.mat.
SHA256: `1c38edc36ba6ff4a41da8f974a88ca1c5ea775655426da91124a35b255a9c870`.
The selected weights were optimized on older TRAIN data; November3 internal
holdouts participated in checkpoint selection. The independently fitted
old+Nov3 replay candidate is retained in mixed_replay/force_model.mat.
Selection preceded external November2 and long-window results; it was not
changed afterward. Neither candidate replaces the earlier release files.

## Protocol and data roles

- Exact warm start from streaming/release_20261103_v1 GRU64/B seed113,
  SHA256 `8150352dd47230a2a4094b240e1e886b641b8c22ebad0cbfe300f63e9f526c6a`. No random network reinitialization.
- One recurrent layer64 units, same .5s training rollout and chunk, .5s
  burn-in, same physical-time guards, nominal mechanics/integrator,80N
  force bound and length/force regularization. Both branches of training
  retain the checked measured-state/independent-free-copy contract.
- Original October19/27 train/selection assignments and59,630-TRAIN-row
  normalizer unchanged. Old128 selection anchors retained exactly.
- November3: allfive supplied recordings audited; 28
  complete TRAIN trials and 9 whole held-out channel/trials.
  2264 new TRAIN chunks,80 new selection anchors. Deterministic
  middle listed trial per channel, the same split rule as the older data.
  Holdouts: circle3; sine1/3/5; helix2; steps3/9/15; combined2.
- Adam restarted because the original released model lacks Adam state.
  Learning rate1e-4 versus original3e-4; gradient clip1; batch32; seed113;
  exactly1200 additional updates per arm, checkpoints every150. New
  checkpoints retain Adam averages/squares, RNG, stream state and queues.
- Old-only arm:32 old streams. Mixed arm:16 old +16 Nov3 streams; equal
  stream counts, not a claim of equal frame counts with variable dt.
- No interpolation, filtering, objective/architecture changes or Nov2
  fitting. Nine Nov3 holdouts are internal previously inspected validation;
  Nov2 is previously inspected external evaluation, not blind confirmation.

Selection: .5*old+.5*Nov3 cohort score. Within each cohort, equal recording
means of .7 one-step + .3 half-second free-path vector RMSE. Older score
must remain <=110% of the initial old score; the initial checkpoint is
eligible so a worse fine-tune cannot be promoted by default. Both trained
candidates improve BOTH initial cohort scores. They trade off which cohort
improves more; only one optimization seed was run, not a significance test.

| Candidate | Extra updates selected | Old internal score (mm) | Nov3 internal score (mm) | Composite (mm) |
|---|---:|---:|---:|---:|
| Original frozen | 0 | 0.54411 | 0.87840 | 0.71126 |
| Old-data continuation | 1200 | 0.41699 | 0.75734 | 0.58716 |
| Old+Nov3 replay | 1200 | 0.51044 | 0.68348 | 0.59696 |

Selected composite improves 17.4%.
Nov3 replay gives the largest gain on NEW internal score, while old-data
continuation gives the larger old-domain gain and slightly lower composite.
This shows that additional optimization mattered; extra data was not an
automatic improvement of every metric under this fixed update budget.

## Whole-record measured-state one-step results

Equal-recording macro q2/q3 vector RMSE in mm. Same fixed initial seed113
checkpoint for before/after, not the earlier three-seed mean0.563mm.

| Cohort | Original | Old-data continuation | Old+Nov3 replay |
|---|---:|---:|---:|
| nov2_external | 0.53588 | 0.36143 | 0.44580 |
| nov3_selection | 0.93307 | 0.75485 | 0.69078 |
| nov3_train | 1.24224 | 1.04670 | 0.70759 |

Selected Nov2 one-step gain: 32.6%.
Nov3 TRAIN metrics are explicitly training diagnostics. Original CV remains
better one-step: Nov2 0.06208mm and
Nov3 heldout 0.10771mm. These are
prediction errors with measured-state refresh, not hardware tracking errors.
Per-coordinate RMSE/p95, numerical guards and force-bound occupancy are
retained in evaluation CSV/MAT files. No failed predictions were discarded.

## Short free-prediction results

Each cell is a macro mean of recording-level path vector RMSE, mm.

| Cohort/horizon | Original | Old-data continuation | Old+Nov3 replay |
|---|---:|---:|---:|
| old_internal, 0.5s | 0.43785 | 0.33577 | 0.42373 |
| old_internal, 1s | 0.56728 | 0.43555 | 0.55449 |
| nov3_internal, 0.5s | 0.77801 | 0.74922 | 0.64048 |
| nov3_internal, 1s | 0.95902 | 0.92145 | 0.85030 |
| nov2_external, 0.5s | 0.37224 | 0.27445 | 0.35224 |
| nov2_external, 1s | 0.43794 | 0.34877 | 0.48304 |

## Complete motion-window diagnostics

Free propagation has one measured initialization and observed recurrent
prefix, then integrated predicted q/v and known recorded pressures/time.
No future measured motion enters the prediction. Missing truth remains
NaN/unscored. All30 candidate/window paths complete the numerical guard.
These long diagnostic windows exceed the .5s training horizon.

**Every Nov3 window here belongs to a whole held-out trial.** Helix uses
trial2 (76.6-134.95s) instead of the prior illustration of TRAIN trial1.
Therefore compare its original4.0294mm to its updated values, not to the
previous trial1 value2.8787mm. Nov2 windows are unchanged from prior plots.
PCC paths are reconstructed from q, not independent measured sensor XYZ.

| Day/window | Original (mm) | Old-data continuation | Old+Nov3 replay |
|---|---:|---:|---:|
| Nov2 circleSweep trial3 | 0.9440 | 0.6229 | 1.7892 |
| Nov2 SineWave4 P1 trial1 | 2.3355 | 1.4374 | 2.5590 |
| Nov2 SineWave4 P2 trial3 | 1.7784 | 1.3042 | 4.2909 |
| Nov2 SineWave4 P3 trial5 | 1.7958 | 1.2315 | 1.5002 |
| Nov2 helix trial1 | 0.7627 | 0.6804 | 1.0095 |
| Nov3 circleSweep trial3 | 1.2636 | 1.0047 | 1.9239 |
| Nov3 SineWave P1 trial1 | 1.7184 | 1.2038 | 2.2712 |
| Nov3 SineWave P2 trial3 | 1.3792 | 1.4367 | 2.9648 |
| Nov3 SineWave P3 trial5 | 1.6876 | 1.7738 | 1.4328 |
| Nov3 helix trial2 | 4.0294 | 4.1306 | 3.2560 |

Selected continuation improves the Nov2 circle and allthree sine windows
and its helix example. The Nov3 circle and sineP1 improve; sineP2/P3 and
held-out helix2 worsen slightly. Mixed replay improves new-day one-step
prediction most but worsens several long rollouts, including both circles
and sineP2. These regressions are retained, not hidden by an average score.
Training/selection focused on measured-state and short free prediction;
improvement does not extend automatically to all long trajectories.

## Offline inverse and verification

Selected model solves 208/208 self-generated
feasible targets. Max force-equation residual 3.5724e-05N;
max pressure-bound violation 2.2204e-16bar.
Median/p95 observed inverse+selected-GRU time 2.044/
4.954ms; 2/208 exceed50ms.
Cold optimizer calls are retained. Identical network architecture and
different machine/JIT load mean these timings are not a claimed training-
induced speedup or qualification for the full50ms hardware cycle. Bounds
are the offline RL0..3bar domain, not new PAC runtime settings. Trial
pressures reuse fixed pre-current hidden state and commit once per sample.
This is model self-consistency, not measured physical force accuracy.

PASS: original weights/releases, old sequences/selection anchors/normalizer preserved.
PASS: training source/data and November2 source/model hashes unchanged.
PASS: architecture, loss, physics, history and normalization unchanged across arms.
PASS: 75 saved one-step rows independently recomputed, max difference 4.88498130835e-15 mm.
PASS:30 long-window rows recomputed, max difference 4.88498130835e-15 mm; all Nov3 motion windows are whole-trial holdouts.
PASS: selected model inverse bounds/state-commit checks.
No hardware/controller/stability/physical force qualification.

## Retained artifacts and downstream scope

- PREDECLARED_PROTOCOL.mat/PROTOCOL.md and data/: source/whole-trial masks,
 208 selection anchors, unchanged normalization, raw/prepared hashes.
- old_only/ and mixed_replay/: all checkpoints, optimizer/stream states,
 training curves, internal selection tables and frozen weight files.
- release/: selected research weights and model card; previous releases
 untouched. Candidate alternative remains in mixed_replay/force_model.mat.
- evaluation/: all three models, Nov2 external, separate Nov3 TRAIN and
 heldout results; one-step and .1/.25/.5/1s predictions/forces/CSV/MAT.
- motion_comparison/: fixed windows, predictions/truth/commands/time,
 scores and six before/after figure bundles. fine_tune_overview is the
 seventh bundle, each editable FIG/PDF/300dpiPNG; all regressions retained.
- inverse/: all208 cases, pressure/state/latency/Jacobian diagnostics.
- training_source/evaluation_source and hashes, verification and diaries.

The downstream goal remains matched PAC and PID FB model-off/on physical
comparisons. Controller gains/laws, SLX and hardware were not changed.
Resolve residual definition versus PAC adaptation, pressure/common-pressure
coverage and measured-velocity/runtime conventions before integration.
No explicit payload estimator, physical parameter identification, passivity
proof or closed-loop tracking improvement is established by this fit.
No migration, Git commit or push. Future confirmation requires fresh data.
