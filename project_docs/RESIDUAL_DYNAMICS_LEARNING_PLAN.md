# Residual learning for missing continuum-arm dynamics

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
