# Residual-model accuracy and PAC / PID feedback performance

2026-11-03. User clarified the primary intended application: augment the
dynamics model used by existing PAC and PID feedback controllers, then measure
whether improved model accuracy improves their physical tracking. This sets
the integration/evaluation objective; no controller implementation, SLX edit or
hardware run is performed by recording this plan.

## Core comparison

| Controller | Baseline model | Augmented model |
|---|---|---|
| PAC | Existing nominal dynamics | Consistent nominal + frozen force residual |
| PID feedback | Existing nominal inverse dynamics | Consistent augmented inverse dynamics |

The primary contrast is within each controller, with the model contribution
switched off/on. Comparing only PAC+residual against PID+nominal would confound
the controller and model effects. The same augmented model definition should
be available to both controller families.

Freeze neural WEIGHTS during these tests. Recurrent memory still evolves from
valid observations at each tick. PAC parameter adaptation, if enabled, remains
enabled with matched gains, initial conditions and reset policy in its pair.
Do not silently add online neural training, gain tuning or a changed filter to
one variant. A later retuned comparison can be a separate best-achievable
performance study; it cannot establish the effect of the model alone.

## Shared model and inverse contract

The current learned model is a generalized-force correction, in N:

    M(q)ddq + (C(q,dq)+D0)dq + K0q + bound(q) + G(q)
        = B*p + r0(q,dq,p,h_previous).

It was trained relative to fixed nominal K0=650*[2 1;1 2], D0=40I,
corrected all-three-muscle bounds and a zero-activation B*p pressure map.
These definitions, gravity, units, coordinates and input normalization must
remain consistent in the controller wrapper, or the residual must be explicitly
re-expressed/retrained for a different nominal decomposition. Simply attaching
these weights to a different stiffness/damping/pressure map is not a justified
model substitution.

The residual takes the candidate CURRENT pressure. Inverse use therefore solves

    B*p + r0(q,dq,p,h_previous) = required generalized force

with a specified common-pressure policy and the actual configured pressure and
slew limits. During pressure search, evaluate every candidate using the SAME
pre-current hidden state. Commit the selected state once per actual sample.
Use matching allocation policy/objective/constraints for residual-off and
residual-on comparisons, so a changed allocator is not credited to the model.
The old offline0..3bar self-consistency test is not a new hardware setting;
trace the active controller configuration before constructing an integration.

Observed q/estimated dq refresh each valid tick (intended PAC period.05s).
Existing measurement validity/recovery behavior must be explicit. The separate
offline interpolation preview is not a live measurement source or an approved
controller missing-input policy. Held-invalid NDI values cannot become truth.

## Controller-specific work before implementation

PID feedback: trace the active feedback path, integral term, derivative/filter
settings and inverse block actually used. Current repository naming includes
PDI/FB/PDFB and a common-pressure inverse accepting virtual acceleration; this
alone does not establish whether the active experiment contains PID or PD.
Use the user's intended PID FB name in the plan, and verify/report the actual
implementation before a comparison. Keep its feedback gains, integral state
reset/anti-windup and reference generation matched in each model pair.

PAC: pacAdaptiveForce currently estimates structured stiffness, damping,
hysteresis-force weights and bias. The learned residual already compensates
for errors in the fixed nominal mechanics and actuator response. Derive a
consistent residual/adaptation decomposition BEFORE adding the residual to
the adaptive force command, to avoid double compensation and a changed
parameter meaning. Check the adaptive update/regressor and residual sign
together. Adding a black-box dynamic residual does not automatically preserve
the original passivity proof. This is a required design task, not a completed
derivation or authorization to change the PAC law here.

## Evidence and trial design

1. Freeze the selected model, normalizer and controller wrapper source. Keep
   independent measured-data prediction evidence, including both one-step and
   free propagation, transient/unloading errors and relevant operating ranges.
   Nominal allocation residual is not measured physical force error.
2. Check the joint operating region: pressure/common pressure, length, velocity,
   trajectory speed/phi, payload, trim and calibration. Individual channel
   pressure bounds alone do not establish coverage of PAC co-pressurization.
   Offline inverse self-consistency is a numerical check, not physical force
   validation. Check complete execution time with room for acquisition/control.
3. Supply reviewable paste-ready model/inverse interfaces when implementation
   is requested. Retain a residual enable switch and log neural force, state
   validity, requested/commanded pressures, solver/limiting status and the
   original controller force/adaptation diagnostics.
4. Operator-run matched trials: same reference, speed, load, calibration/trim,
   start pose, gains, filters, pressure/common-pressure limits and duration.
   Fresh state reset per trial; repeat trials and interleave model-off/on order
   to reduce effects of warm-up, drift and loading history. Controller settings
   and data filenames/split roles must be recorded before each campaign.
5. Score q2/q3 RMSE and velocity RMSE, transient/cycle convergence, repeatability,
   smoothness/control effort and pressure limiting. Score theta/phi with proper
   observability/wrap rules and sensor-compensated XYZ when the reference/frame
   are available. Retain full learning curves as well as steady-cycle scores;
   PAC adaptation can change convergence and conceal early differences.

The previously used12-circle/last3 scoring protocol can be a starting point
where it matches the operator's active setup, but cycle count and settling
coverage must be checked rather than copied from historical settings. Do not
invent new acceptance limits from prediction RMSE or one selected plot.

The claim to test is: at matched controller settings and operating conditions,
the model augmentation reduces physical tracking error without an unacceptable
smoothness/pressure-limit cost. Better offline predictions motivate that test;
they do not establish its result. Keep negative results and both controller
pairs, not only the most favorable combination.

## Current status and next dependency

Frozen GRU64 has November2/3 diagnostic plots; proposed continued training is
old October TRAIN + November3 TRAIN, preserving internal whole-trial holdouts
and November2 evaluation-only, then a fresh independent confirmation cohort.
Fine-tuning is now complete: fine_tuning/campaign_20261103_v1/REPORT.md under
force_model. The equal-cohort internal score selects old-data continuation;
old+Nov3 replay is retained. Old-only gives larger Nov2/circle gains, while
replay gives lower Nov3 one-step error but worsens several long windows. No
controller integration follows automatically from this checkpoint selection.
Shared inverse/PAC consistency, operating-region and physical qualification
remain pending. All arm experiments are operator-run; no autonomous hardware.

## November3 controller-branch handoff

The user chose controller implementation and a first matched comparison before
building HIL model-update automation. Start PDFB, then PAC after deriving the
adaptation/residual decomposition. Use the128-unit short-horizon candidate
and preserve the64-unit reference; neither is controller-qualified. Exact
paths/hashes, tradeoffs and required sequence are in RESIDUAL_CONTROLLER_HANDOFF.md.
Branch requested: codex/residual-controller-integration. No hardware or SLX
changes are implied by creating the branch.
