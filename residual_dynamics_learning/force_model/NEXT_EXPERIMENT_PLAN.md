# Proposed next force-GRU experiments

2026-11-03. User proposes two GRU layers and training rollout <=1s because
the intended application estimates current missing force for pressure correction.
User authorized executing this sequence on November3 ("Great. execute the
sequence"). Implementation and offline fits are now in force_model/streaming;
see its README and versioned campaign/evaluation reports for execution status.
This is not a model replacement or hardware schedule.
Current release_20261103_v1 and its original evaluation remain frozen.

Execution complete November3: 11 offline fits, old-only model selection and
frozen November2 evaluation; bounded model-inverse self-consistency checked.
See streaming/REPORT.md. Selected 1x64/.5s, not stacked32+32; CV still better
through1s. No physical-force/control qualification. The new research checkpoint
is streaming/release_20261103_v1; the parent release remains unchanged.

## 1. Align the model with the intended online use

Distinguish three time spans:

- Past observed history/burn-in initializes recurrent state. Longer past history
  is not forecasting farther ahead.
- Training loss horizon determines how integrated future length errors guide
  force learning. A short multi-step horizon is useful even if only current
  force is needed during control.
- Online output is r_tau at each actual PAC sample, intended .05s, based on
  then-available measured q/estimated v, pressure information and recurrent state.
  A .5s training horizon does not require .5s of future measurements at runtime.

The current release warms hidden state with ten observed frames, then evaluates
free-running integrated q/v. That validation mode is not identical to a live
measurement-conditioned estimator which refreshes q/v every tick. Before blaming
width/depth, define and test the state/velocity contract under measured-input
replay. Check velocity/filter conventions and units against the integrated state.

Proposed new baseline: causal measurement-conditioned state updates on valid
segments, one update per actual sample; initialize/reset explicitly at run/gap
boundaries and use the same warm-up policy in training and replay. For training,
use ordered segments with truncated backpropagation through short chunks; hidden
state can carry across chunk boundaries, with the gradient graph detached there.
Random independent windows must not silently masquerade as persistent streaming
state. No look-ahead references or future measured motion enter current-force
inputs. Both layers, if stacked, retain their own hidden state.

Retain separate primary measured-state one-step evaluation and secondary short
free-running evaluation. Neither is a closed-loop tracking or measured-force
accuracy test. CV is a kinematic sanity baseline, not an alternative residual
force estimator or pressure correction law.

## 2. Fixed objective and six-model first screen

Rebuild the 1x32 control under the new streaming/objective contract so changes
from the existing release are not falsely attributed to extra layers.
Suggested starting loss, held identical across the capacity/horizon screen:

    .5 * normalized measured-state one-step length loss
  + .5 * normalized short-rollout length loss
  + force-rate penalty + force-magnitude penalty

The two length terms are independently averaged/normalized so the longer term
cannot overwhelm the first interval merely by containing more samples. The
50/50 weights are an initial proposal, not established optimal values. Use old
development trials to adjust if necessary BEFORE the six-model comparison;
do not select them using November2. No force/acceleration truth labels invented.
Keep current small force regularization as the starting comparison, no post-hoc
curve filtering. Force-rate differences refer to chronological adjacent samples.

| Candidate | GRU architecture | Maximum physical training rollout |
|---|---|---:|
| A | 1 layer,32 units | .5s |
| B | 1 layer,64 units | .5s |
| C | 2 layers,32 units each | .5s |
| D | 1 layer,32 units | 1.0s |
| E | 1 layer,64 units | 1.0s |
| F | 2 layers,32 units each | 1.0s |

Stacking means layer1's causal hidden/output features feed layer2. It does NOT
automatically assign actuator physics to one layer and mechanics to another.
Report parameter count as well as depth/width; these candidates are practical
capacity comparisons, not exactly parameter-matched proof of a depth effect.

Use ACTUAL poll-start physical time to define .5/1s. Existing24steps are only
nominally1.2s; October27 sample intervals vary. Select the last eligible sample
within the physical duration cap and retain actual durations, not a blanket
10/20frame assumption. Keep masks/command-age limits unchanged.

Hold data/splits, input normalization, initial observed-history policy, force
bound, nominal model, integrator, loss weights, seed and optimizer policy fixed
across the six configurations. Use matched valid anchors based on the longest
needed history/horizon; otherwise shorter windows can change the sample set.
Compare equal optimizer-update budgets with early stopping on development data,
and report wall time too. Start near the selected .0003 learning rate; learning
rate adjustment is a separate development choice, not a hidden advantage to a
preferred architecture. November2 remains evaluation only.

## 3. Test history and seed sensitivity after the first screen

For the leading architecture, compare .5s versus1s observed burn-in under matched
anchors. This tests initialization separately from recurrent width/depth and
training horizon. Streaming state after burn-in can remember longer than the
burn-in itself. If poor short-horizon performance persists, inspect velocity/
state consistency and loss balance before expanding capacity further.

Repeat the leading two configurations with three fixed seeds total each;
reuse the screening seed and add two repeats. Report mean/spread and keep each
checkpoint, not only the lucky seed. No exhaustive architecture search proposed.

## 4. Data split and ranking

Existing October19/27 whole-trial training/internal selection remain the
development starting point. New speed/phi/rise-time recordings can augment a
future training pool after file/settings/quality audit, but assign their role
BEFORE fitting. Reserve complete runs/conditions and fresh repeats for validation;
keep day,calibration,pressure/mounting/load/source metadata. Do not randomly
split neighboring samples. Current November2 recordings remain evaluation-only
and are acknowledged as previously inspected, not blind confirmation.

Rank by development-set one-step and .1/.25/.5/1s length prediction, including
per-coordinate RMSE and95th-percentile error. Report angle/path diagnostics
where observable. Evaluate loading/unloading and speeds/phi/rise times separately
so quiet holds do not hide failures. Measure residual-force rate, spikes and
bound occupancy; no measured residual-force accuracy claim is possible.
Check streaming initialization/reset sensitivity and complete inference timing.
Retain optional longer free rollouts only as drift diagnostics, not the principal
criterion for a current-force estimator. Final external evaluation happens once
after model/normalizer/objective selection and freeze; no post-test tuning.

## 5. Pressure correction is a separate subsequent test

The current residual takes CURRENT pressure as an input. Therefore solving

    alpha*A*p + r_tau(q,v,p,h_previous) = nominal_desired_force

requires a bounded inverse calculation and a common-pressure choice. Keep the
committed hidden state fixed while trying pressure candidates; commit recurrent
state only once per actual applied-command/sample step. Do not advance hidden
state once per optimizer/root-solver iteration. A past-pressure-only residual
would be an alternative contract requiring separate training/validation, not
just shifting inputs of current weights.

First test inversion feasibility, pressure/rate bounds, continuity and total
latency offline under a fixed contract. Include model inference plus inversion
in the intended .05s execution budget and leave room for acquisition/control.
Self-consistent inverse success is not evidence of physical tracking improvement.
Any later arm/controller comparison remains operator-run and separately requested;
this plan does not modify SLX, run hardware or deploy a residual controller.
