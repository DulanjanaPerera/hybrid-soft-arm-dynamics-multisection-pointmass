# Residual learning for missing continuum-arm dynamics

2026-11-03 user adds two-layer GRU and short <=1s rollout proposals. Detailed
proposed comparison and causal online-estimator contract are now in
residual_dynamics_learning/force_model/NEXT_EXPERIMENT_PLAN.md. Six matched
configs:1x32,1x64,2x32 versus .5/1s physical-time training horizons; primary
measured-state one-step loss/replay, secondary short integrated rollouts,
history and seed follow-ups. No network/fit changed for this planning discussion.
Training horizon is not runtime look-ahead or past-memory length. State must
advance once per actual sample; pressure inversion holds state fixed across
candidate evaluations. Old-only/development selection and current November2
evaluation-only separation retained. Additional data split must be declared
before fitting. Keep existing release as a frozen reference.

2026-11-03 user planning update: additional circle speeds, bend-angle phi
values and step rise times are being collected for later training. Audit raw
settings/timestamps/calibration/validity when supplied, retain run/day identities
and assign whole-run/condition development and held-out sets before fitting.
The existing November2 evaluation-only contract remains; no pooling or new
training performed for this discussion. A capacity increase can
be tested, initially32 versus64 hidden units under matched data/loss/horizons;
it is not an established fix for short-horizon error. Compare history/burn-in,
rollout horizon, initialization, learning rate and objective weighting as well.
These are proposed later experiments; current weights/settings remain frozen.

Deferred workspace migration requested for later: heavy residual-model/data
artifacts to C:\Lab\Research\Residual modeling outside OneDrive, preferably
with agents/training stopped. No move or nighttime automation performed now.
Plan hash-verified copying, path/reproduction updates and preserved original
provenance before any cleanup; do not break active MATLAB/agent workflows.

2026-11-03 force-model training implemented/completed on explicit user request.
Use residual_dynamics_learning/force_model/release_20261103_v1/REPORT.md.
October19+27 training/internal whole-trial selection; November2 used only after
exact model freeze for evaluation. New GRU32 outputs two generalized forces in
N inside the corrected-mode1 mechanical ODE; integrated-length loss, force-rate
and magnitude penalties, no acceleration/force truth labels. This supersedes
the earlier force-model-unimplemented status. No existing position-model weights
reused or recast. Old normalization and raw datasets remain unchanged.
Selected lower-rate epoch38 using older holdouts only; old selection2.1103mm.
Held-out2/5s path1.8375/2.1267mm versus nominal7.0438/6.6372 and CV.8927/3.4062.
One-step2.151mm is weak. Next investigate initial-state/velocity consistency,
capacity and optimization on OLDER development data; do not tune on November2
or proceed directly to inverse/PAC/PDFB use. Input-dependent residual inversion
also remains to be implemented. No physical-force identification/stability claim.

2026-11-03 retained user preference: address GRU prediction spikes in later
model development. First distinguish nominal-transition spikes from learned
correction spikes; smoothing the correction alone can weaken cancellation.
For the intended generalized-force residual, prefer integrated-length multi-step
training plus a modest temporal smoothness penalty on the residual force, using
physical dt and adjacent ordered samples within valid segments. Do not penalize
differences between unrelated shuffled batch examples. Preserve genuine fast
loading/unloading behavior; choose penalty strength using development data.
Consider a light causal residual-output low-pass filter only if still needed,
and include its state/delay in end-to-end rollout/inverse validation, especially
at fast circle speeds. Plot-only/noncausal smoothing is not a model improvement.
Root cause and benefit remain unverified. This is a recorded proposal for later,
not authorization to retrain the old position-correction GRU or change filters
now. No code/model/training settings changed by recording this preference.

2026-11-02 fresh-data validation completed. See validations/fresh_20261102_v2
under residual_dynamics_learning. Five November2 experiments, 50,078 guarded
anchors, unchanged frozen weights/normalizer, exact raw-export equality. Ridge
one-step macro .05976 vs constant velocity .06235 mm; 5 s autonomous path
1.6204 vs 3.4076 mm, but kinematic ablation is competitive and GRU is favorable
at 5 s despite poor one-step scores. No universal architecture/prior advantage
or controller readiness established. No training/model selection on this set.
This completes the proposed fresh-cohort evaluator adapter and first frozen
evaluation; next implement the explicitly versioned generalized-force residual
with integrated-q loss, retaining baselines and train/test isolation. November2
has now been observed and must not become an allegedly untouched test for later
test-informed changes. Current weights remain position corrections, not forces.

2026-11-02 concrete nominal correction: equal nonlinear bounds now include all
three physical actuator lengths, still+-35mm. This fixes an identified omitted
l1 force rather than fitting nominal parameters. See
RL_data_collection/ACTUATOR_BOUND_CORRECTION.md. Current model source defaults
to actuator-bounds-v2 (mode1); old two-day screening/preparation and calibration
remain coordinate-bounds-v1 (mode0) explicitly for reproducibility. Existing
residual weights/context must not be switched to the new baseline implicitly.
Future corrected residual datasets/models/reports require a new version with
explicit baseline provenance. No new residual trained or hardware run.

Priority clarification2026-11-02: nominal inaccuracy is already established;
the calibration/actuator-memory trial was a bounded offline diagnostic, not a
requirement to fix the nominal model before residual learning. User challenged
that detour. Proceed with an explicit frozen nominal baseline plus a learned
history-dependent missing-dynamics residual; proposed placement is generalized
force, fitted through integrated length predictions to avoid noisy acceleration
truth assumptions. This proposed force-residual implementation differs from
the earlier screening's next-position residual and has not yet been trained.
Retain matched baselines, whole-trial splits, timing/quality controls and longer
rollout/loading-unloading validation before model/controller integration.
Do not extend nominal calibration as the main workstream without a concrete
obstacle that requires it. This supersedes the earlier next-calibration priority.

Retain potentially paper-useful figures and numerical/source artifacts in
versioned folders. See project_docs/FIGURE_INDEX.md: editable FIG, vector PDF,
300dpi PNG, source data/code/provenance, captions and limits on claims. Preserve
negative results and existing outputs; distinguish local retention from Git
commits and from publication-qualified evidence. No manuscript edits requested.

Update2026-11-02: constrained static/dynamic calibration and shared actuator-
memory diagnostic executed offline. See calibration/greybox_20261102_v2/REPORT.md
under residual_dynamics_learning. Day1 validation~2s pathRMSE original5.136mm ->
mechanics3.369 -> memory3.295, but matched constantvelocity1.055mm is better.
D reaches both upperbounds, delay/T confounded; static validation still5.583mm.
No model qualified for hardware or controller changes. Next inspect static
residual versus pressure/common mode/loading and command/state timing. Use that
diagnosis to constrain an effective actuation model or design informative
operator-run data before increasing residual capacity. Oct27 untouched by this
calibration, but already evaluated in step3 and not fresh confirmation.

Step3 update2026-11-01: preliminary offline benchmark completed, report under
residual_dynamics_learning/benchmarks/screen_20261101_v2. Ridge wins initial
validation, but cross-day forward and pressure-only inverse ablations prevent
qualification as improved missing-dynamics models. Do not proceed directly to
hardware step4 merely because weights exist. Next refine model identification,
prior/actuator-memory and informative data; collect fresh confirmation later.

## 2026-11-01 agreed execution sequence and step1 result

Step2 completion update: valid prepared dataset under
residual_dynamics_learning/datasets/two_day_v2, with its README as executable
data contract. Day1 train35,350/validation13,977 and day2 external24,079 shared
eligible anchors. No train/test pooling, invalid interpolation or statistics
fit on day2. Nominal RHS and recorded map-convention checks completed; inverse
candidate targets are derived and not force truth. Sparse day2 pressure-run
windows/timebase mismatch are limitations for step3, not reasons to loosen
quality thresholds after seeing results. No model trained or hardware run.

User proposes: (1) define input/output roles and architecture candidates from
the RL experiments, including inverse three-pressure prediction; (2) prepare
RL_data_collection for training; (3) offline nominal-plus-residual preliminary
tests; (4) operator-run forward/inverse validation; (5) targeted additional
speed/trajectory/loading data; (6) retrain and revalidate with held-out tests.

Step1 executed: see residual_dynamics_learning/MODEL_SPECIFICATION.md and
audits/interfaces_20261101. Two days explicitly preserved: root recordings
2026-10-19 and subfolder20261027 recordings2026-10-27, five each, dates verified
from saved acquisition metadata. Physical commands identified by final NI-path
provenance, not nominal inverse outputs. Forward and inverse residual contracts
are separate; inverse pressure outputs require a common-pressure policy because
two nominal generalized inputs do not uniquely specify three pressures.

Shortlist: zero residual, simple memory model, windowed MLP, small GRU; LSTM
only if empirical validation supports it. Architecture winner awaits training.
No raw recordings merged/modified, SLX edited, hardware run or new commit made
for this step. The original plan below is retained; its blanket reference-input
restriction applies to forward dynamics, while inverse preview is explicitly
permitted under the new specification's availability and label rules.

Created2026-10-28. Branch: codex/residual-dynamics-learning.
Parent checkpoint:5532a5d on codex/pac-autotuning-framework. This document starts
the research branch; no residual model has been fitted or deployed.

## Objective and scope

Learn missing dynamics as an augmentation of the existing model, then measure
prediction improvement and the effect on PDFB and PAC. The user's explicit
scope is residual MODEL learning, not a separate residual feedback controller.
Use established manually started experiments; defer continuous online gain
tuning. Do not substitute a learned error-to-command policy for a dynamics model.

A candidate representation is

    x_next = F_nominal(x, u, dt) + r_psi(x, u, input_history, dt)

where x includes the chosen mechanical state and, if justified, actuator memory.
This is a modeling template, not an implemented discretization. Define units,
states, integration and the residual's placement only after auditing the actual
nominal dynamics. A force residual and an acceleration/transition residual have
different targets and interpretation; do not silently interchange them.

## Existing implementation and data limitations

- Mechanical q=[l2;l3], l1=-l2-l3. Three pressure commands are actuator inputs.
- Lengths are inferred through NDI orientation/PCC geometry, not independent
  chamber strain measurements. Sensor/geometry error can appear as model error.
- No explicit actuator-lag/uncertainty model or payload-mass estimator exists.
  Actuator memory is a plausible first model extension, not an identified cause.
- Pressure commands are logged; delivered actuator force and chamber pressure
  are not measured. Nominal allocation residual is not a physical-force label.
  Do not train against it or tauRequested as if either were measured force.
- Adaptive stiffness/damping/history/bias estimates are compensation values,
  not uniquely identified true mechanical parameters. Separate these from the
  fixed nominal model used as the learning baseline.
- Intended PAC period .05s; audit actual measurement timestamps, physical poll
  time, repeated/missing frames and command alignment before differentiating.
  Existing filtered acceleration is noisy and is not automatically ground truth.
- Recent 5s smoke R033 completed motion but had tip-validity dropouts and rejected
  score coverage. Its control results are exploratory; don't use it as clean
  validation without explicit validity masks and purpose-specific quality checks.

## First task: data and nominal-model audit

1. Read PROJECT_STATE, TUNABLE_BLOCK, relevant VALIDATION_RESULTS sections and
   this plan. Inspect actual source and model files before using old settings.
2. Trace the implemented nominal M/C/G, stiffness, damping, history, bias and
   pressure-to-force terms. Identify the actual PDFB implementation and the
   PAC model interface; do not infer either from an acronym or paper equation.
3. Inventory available MAT runs: inputs, states, timestamps, sensor validity,
   gains/scales/bounds, payload, calibration/trim, saturation and adaptation
   switches. Record hashes and distinguish operator attestations from logs.
4. Define train/validation/test runs by whole experiment. Fix time alignment,
   units, causal filters and invalid-gap rules. Avoid random adjacent-sample
   splits, held-out preprocessing leakage and interpolation through dropouts.
5. Assess whether data can identify command-memory versus mechanical mismatch.
   If not, specify a small manual excitation experiment rather than fitting a
   high-capacity residual to an uninformative single closed-loop circle.

Deliver a reproducible audit and modeling decision before choosing architecture.

## Learning and offline validation

Start with zero-residual nominal prediction as baseline. Consider a low-order,
stable command-memory/actuator-state extension before a flexible residual if
the audit supports it. Without pressure/force sensing, call inferred lag an
effective model; do not claim unique physical identification.

If adding a learned residual, begin with the smallest justified feature set
and capacity. Dynamics inputs can include q, dq, pressure commands/history and
logged gravity/known conditions. Do not add tracking error, desired trajectory
or controller output decomposition as shortcuts to learning the controller.
Any such feature needs a physical modeling justification and explicit scope.

Prefer prediction/output-error fitting using valid position sequences and
integration when acceleration labels cannot be trusted. If derivative-based
targets are used, quantify filtering/differentiation effects and label them
as derived targets, not measured force. Treat closed-loop correlation and
mechanical/actuator confounding explicitly.

Select hyperparameters on training/validation runs only. Freeze weights,
normalization, filter parameters and feature definitions before held-out testing.
Assess one-step AND multi-step/rollout prediction, finite bounded behavior,
speed/direction coverage and uncertainty outside the training domain. Include
zero-residual and simple memory-only ablations where informative. Improved
one-step fit alone does not establish improved dynamics or stable control.

## Controller comparison after prediction validation

Use the same frozen learned dynamics augmentation for the paired model variants:

| Controller | Original model | Augmented model |
|---|---|---|
| PDFB | Baseline PDFB | PDFB with learned dynamics |
| PAC | Baseline PAC | PAC with learned dynamics |

Initially fix feedback gains within each controller's original/augmented pair,
filters, pressure bounds/rates, reference, payload, reset/calibration and scoring
windows. Match PAC initialization and adaptation settings. Do not conflate an
improved model with a simultaneous gain retune. PAC adaptation may overlap the
residual's compensation: define the decomposition and inspect double counting.

Begin with a repeatable10s circle protocol, then use5s as a harder case after
checking sensing quality. Select actual baseline settings from verified logs;
the scale1 historical tuning and current scale5 smoke settings are different
configurations and neither is universally optimal. No trial is queued here.

Report length RMSE per coordinate, velocity error, fixed-analysis roughness,
pressure limiting/faults and repeatability; XYZ/theta/phi are secondary. Report
prediction improvement separately from control improvement. Use fresh repeats
and held-out conditions; separate-run improvements do not by themselves prove
adaptation benefit. Do not inherit a PAC passivity/stability claim merely from
a learned model's prediction accuracy; document new assumptions and bounds.

## Workflow and handoff

First implementation milestone: offline dataset audit, nominal forward predictor
and held-out baseline report. Next: justified residual model fit and offline
validation. Only then prepare paste-ready controller/model interfaces for
operator-run comparisons. Residual fitting and controller integration are not
implemented by creating this branch.

Never run hardware models autonomously. Do not edit SLX unless asked. Preserve
unrelated RL edits, staged rename, paper drafts and original raw recordings.
Update PROJECT_STATE and a specialized residual-model log after substantive
findings. Commit/push only when asked; this branch was created without a push.


## 2026-11-03 - Proposed bounded missing-NDI training experiment

Current frozen preprocessing is unchanged. New data audit and explicit preview:
residual_dynamics_learning/figures/motion_predictions_20261103_v1/short_gap_audit.
53 of57 missing samples meet a one-frame<=.12s physical bracketing-span rule;
artificially hiding real observations gives q-vector RMSE.040-.098mm. This is
feasibility, not actual missing truth or evidence of training improvement.

For a subsequent separately versioned training experiment, compare strict gap
resets with bounded one-frame interpolation and/or gap-aware physics carry.
Keep raw measurements, timestamps, NI pressure ZOH, quality flags and explicit
imputed masks. Use original real-only observation losses; never invent force
or acceleration ground truth from filled lengths. Preserve trial/calibration/
command/time and longer-loss boundaries. Split whole trials/runs BEFORE fitting
and evaluate with unchanged real-only masks. November3 has now been inspected
for quality/interpolation but has not been fitted or used for model-score
selection; declare any future train/held-out roles explicitly.

Two-sided interpolation uses the next valid endpoint. Do not claim zero-delay
online availability or let future-filled inputs leak into validation recurrent
states; a causal predictor/missingness policy requires separate training/tests.
If raw rotations are later interpolated, use quaternion SLERP rather than
wrapped Euler-angle interpolation. Smoother signals do not imply accurate
velocity/acceleration or recovered missing physical dynamics. No new missing-
data runtime/training policy is authorized as an implemented default here.


## 2026-11-03 - Recommended continued-training experiment (not executed)

Prefer November3 as added development data because it extends recorded speed/
amplitude coverage, while retaining November2 evaluation-only. Nov3 now has
frozen-model prediction plots and a coverage audit, so it is inspected data;
whole-trial holdouts are internal validation, not pristine blind confirmation.
See figures/motion_predictions_20261103_v1/frozen_model_predictions_v1/REPORT.md.

Before any fitting, declare complete trial/condition TRAIN versus internal
holdouts within November3. Preserve original October19/27 split assignments.
Start from frozen GRU64/B seed113 with original normalizer, architecture,
.5s horizon and burn-in, same nominal integrator and length/force penalties.
Use a lower learning rate and replay original TRAIN examples together with
new TRAIN examples. Check old and new holdouts for both improvement and
regression. Do not silently fit a new normalizer under old weights or claim
exact optimizer continuation: release lacks Adam state; declare optimizer reset.
A new versioned continuation trainer and manifests are necessary; preserve the
old trainer's old-only whitelist/source and do not overwrite its dataset/release.

Matched-budget comparison: frozen model; extra updates on old TRAIN only;
extra updates on old+Nov3 TRAIN replay. Select checkpoint/objective/hyperparameters
using old+Nov3 internal holdouts, never Nov2. Freeze first, evaluate Nov2 once
for before/after comparison. Keep it labeled previously inspected evaluation.
Retain a subsequent independent repeat/day for final confirmation. If instead
using Nov2 for training, its subsequent plots are training diagnostics and a
different held-out cohort is required. Dates alone do not determine data value.

Keep interpolation, preprocessing/filter changes and architecture changes out
of this initial data-addition experiment. They are separate controlled ablations.
More data/updates could help, especially faster/larger circle and helix dynamics,
but neither improvement nor better closed-loop hardware tracking is guaranteed.
Neither these10s nor20s circle datasets establishes performance at2s circles.

General method references: forgetting https://doi.org/10.1073/pnas.1611835114 ;
replay https://proceedings.neurips.cc/paper_files/paper/2019/hash/fa7cdfad1a5aaf8370ebeda47a1ff1c3-Abstract.html .


## 2026-11-03 - Primary downstream application: PAC and PID FB

User clarified the purpose of residual model accuracy improvements: augment the
shared dynamics used by PAC and PID feedback controllers, then demonstrate the
physical tracking effect in matched model-off/on comparisons. Detailed plan is
RESIDUAL_CONTROLLER_COMPARISON_PLAN.md. Freeze neural weights during comparison,
match controller gains/filters/limits and PAC adaptation initialization/update
settings, and retain transient/converged/roughness/limiting measurements.

Required before implementation: trace the actual PID/integral/feedback-linearizing
path; reconcile fixed nominal residual definitions and pressure mapping; derive
PAC adaptation/residual consistency to prevent double compensation; bounded
pressure-dependent inverse with one recurrent commit per sample; assess the
actual pressure/common-pressure operating region. Better prediction is evidence
supporting this experiment, not a guarantee of better closed-loop tracking or
preserved passivity. No controller code or Simulink edits are made by this plan.


## 2026-11-03 - Fine-tuning execution completed

The authorized comparison is complete in force_model/fine_tuning/
campaign_20261103_v1/REPORT.md under residual_dynamics_learning. Exact weight
warm starts, fixed loss/physics/normalizer, fresh Adam1e-4,1200extra updates per
arm. Original split and nine Nov3 whole-trial holdouts preserved. Predeclared
guarded equal-old/new internal score selects old-only continuation; mixed Nov3
candidate retained. Mixed reduces new one-step error more but regresses some
long free trajectories; selected continuation improves circles and allfive
Nov2 long examples, with some Nov3 long regressions. Keep allresults.

No Nov2-informed refit/reselection. No interpolation, capacity, preprocessing,
controller or hardware changes. This completes the prior proposed continuation
experiment; it does not implement PAC/PID FB integration. The next model study
could investigate the short-loss/long-rollout trade-off using development data,
with a fresh confirmation cohort for resulting changes. Do not treat improved
prediction RMSE or inverse self-consistency as proven physical tracking gains.
