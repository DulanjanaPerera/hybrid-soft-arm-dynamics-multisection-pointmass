# Residual dynamics MODEL with PDFB and PAC

2026-11-09: for pressure prediction for **desired lengths/velocity/acceleration**,
use [DESIRED_PRESSURE_BLOCK.md](DESIRED_PRESSURE_BLOCK.md) and paste
rlResidualDesiredPressureBlock.m. Input2 is now6x1 desiredMotion; both pressure
branches share the desired-state force target. Existing diagnostic bus is
reused; observed-state historical blocks remain available for replay.

2026-11-06: updated passive single-pass pressure conversion is documented in
[SINGLE_PASS_COLLECTION_SHADOW.md](SINGLE_PASS_COLLECTION_SHADOW.md). Paste
rlResidualSinglePassCollectionShadowBlock.m (same5inputs/3outputs/bus) for
the collection minimum-pressure convention. The older PAC-policy observation
and all actual controller wrappers remain unchanged.

For a separate passive observation block in RL_dataCollection_v2, start with
[SHADOW_VALIDATION.md](SHADOW_VALIDATION.md). Paste rlResidualShadowBlock.m;
optional rlResidualInverseShadowBlock.m observes the inverse without actuating.
Single-pass subtraction alternative: [SINGLE_PASS_SHADOW.md](SINGLE_PASS_SHADOW.md)
and rlResidualSinglePassShadowBlock.m; one GRU evaluation at actual pressure,
existing bounded allocation, separate diagnostics, no iterative pressure solve.
The user installs/wires blocks manually. No SLX modification script is provided.

2026-11-03, branch codex/residual-controller-integration. Source-level
implementation and offline component verification; no SLX edits, hardware
run, live controller qualification, training or automatic deployment.

The outer feedback, integral term, trajectory and filters stay outside this
model/inverse update. The learned force depends on candidate current pressure,
so an exact inverse must account for pressure dependence. A separate passive
single-pass block now tests the approximation of freezing that residual at
actual pressure before allocation; it is not an exact nonlinear inverse.

## Prepare once, with the model stopped

From the repository root, in the MATLAB base workspace:

```matlab
addpath(pwd,fullfile(pwd,'residual_controller'));
[residualOptions,pdfbPressureOptions,pdfbSupervisorOptions, ...
    ResidualControllerDiagnosticsBus] = setupResidualControllers(pwd);
```

This adds paths, checks the original64/128-unit model SHA256 hashes, prepares
fixed-size numeric weight packets and creates the returned diagnostic bus.
It does not load/update/start a Simulink model or connect hardware. Original
model weights are unchanged. Keep the MAT packets and all helper .m files on
the MATLAB path. Both packet filenames are constant for `coder.load`, so weights
are embedded during code generation, not trained or loaded from disk every tick.
See [MathWorks coder.load](https://www.mathworks.com/help/coder/ref/coder.load.html).

Set the MATLAB Function block explicitly to discrete0.05s. Keep the existing
outer feedback/reference/filter/physical Kill wiring. Set the new diagnostic
output type to `Bus: ResidualControllerDiagnosticsBus` in the block's Symbols
pane and use a Bus Selector/logging for its fields. The setup output must be
assigned to that exact base-workspace bus variable name.

## PDFB: FBL_commonP_circle_BRot.slx

Read-only ZIP/XML inspection found a six-input inverse and an outer virtual
acceleration law including `Ki.*ie`. Preserve that feedback/integrator.
Paste **pdfbResidualInverseBlock.m** into a duplicate inverse MATLAB Function.
The first six input ports and first eight output names/shapes match that saved
inverse. Append inputs7-13:

| Input | Type/size | Connect |
|---|---|---|
| measurementTime_s | double scalar | Existing aligned physical NDI sample timestamp; same clock origin as t |
| enablePressure | Boolean scalar | Existing pressure permit, consistent with downstream Kill |
| enableResidual | Boolean scalar | false initially; true only for a fresh model-on trial |
| enableCommonPressureAdaptation | Boolean scalar | false for the saved fixed-common PDFB baseline |
| pressureOptions | double10x1 | pdfbPressureOptions |
| supervisorOptions | double7x1 | pdfbSupervisorOptions |
| residualOptions | double8x1 | residualOptions |

First six inputs remain measured,desiredAcceleration,orientationValid,t,
gravity,baseOrientationValid. Wire the same inferred-length vector and physical
gravity/base-valid signals as before. No change to the outer Kp/Kd/Ki or errors.
Append output9 `status` (uint8) and10 `residualDiagnostic` (the bus).

**Paired nominal baseline:** model-off/on in this NEW wrapper share the same
pressure supervisor/allocator and corrected all-three-actuator bound model.
The saved embedded inverse still had the old l2/l3-only bound. This wrapper
uses the existing PAC supervisor/activation-aware allocator: it is NOT a
bit-identical replacement of every embedded legacy pressure-path behavior.
In particular, allocation/blending, timestamp freshness, explicit enable and
clock-reversal fault behavior must be qualified in a fresh model-off trial.
Do not compare the old SLX baseline directly against new-wrapper model-on and
attribute all differences to neural accuracy.

The inspected saved PDFB inverse had common1.6bar, trims
-0.0983981693363845/-0.0819816714043049bar, rails.2..3bar, s/a caps2.5bar,
slew2bar/s, pressure start10s, feedback15s, blend5s and stop220s. Its software
force map used linear pressure differences, hence the PDFB template's activation
entry is0. The template anchor-speed1e9 represents the saved absence of a speed
gate. These are audited SAVED settings, not confirmation of currently running
parameters. Review them and the stop time before installation; no gains or
reference settings are overwritten by setup. The old eight-output
`tauAvailable_N` port now reports calibrated model-equivalent force, including
the neural correction when enabled; it is not a measured chamber-force signal.

## PAC: existing pacAdaptiveTunableBlock interface

Paste **pacResidualAdaptiveTunableBlock.m** into a duplicate MATLAB Function.
Retain the current23 input connections and their actual logged tuning values,
including freezeHysteresisWeights. Append:

- Input24 `enableResidual`, Boolean.
- Input25 `residualOptions`, double8x1.
- Output29 `residualDiagnostic`, Bus: ResidualControllerDiagnosticsBus.

An optional paste-ready residualControllerTuningBlock.m supplies those two
inputs, or use Constant blocks. Its default is model-off. The existing first28
outputs are exactly identical to the original PAC wrapper in model-off mode;
this was tested with adaptation, both common-pressure modes and dropouts.
The previous PAC gain/tuning source is NOT replaced. In particular, do not
copy historic gamma/scale/pressure settings over the live configuration.

Neural weights stay frozen. PAC's existing10-parameter update/regressor,
projection/rate caps, history evolution and hysteresis-freeze switch are retained.
The estimates now provide remaining compensation around the learned baseline.
The original absolute parameter port remains for compatibility; diagnostic
parameterCorrection=theta-[650;650;650;40;40;0;0;0;0;0]. Start fresh from that
nominal theta0 when testing the learned model; do not carry learned estimates
from a previous nominal-only run. See MODEL_INTEGRATION.md for the derivation
and its limits; the old passivity proof is not automatically preserved.

## Residual options (before a fresh run)

| Row | Default | Meaning |
|---|---:|---|
|1|128|Frozen network units;64 selects the retained reference|
|2|6|Maximum damped Gauss-Newton iterations, at most12|
|3|.02|Augmented inverse residual tolerance N|
|4|1|Valid observed warmup seconds before capturing an augmented anchor|
|5|.10|Maximum gap s before resetting neural history|
|6|.001|Finite-difference step N for the existing allocator's pressure projection|
|7|20|Maximum Newton command step N, still subject to all pressure constraints|
|8|.05|Intended period s; initial feature dt only, later actual sample intervals used|

These are numerical/runtime settings, not new adaptation gains or a physical
force-accuracy acceptance threshold. Keep the default architecture/configuration
fixed within each trial. An active mode/network change after anchoring faults
the pressure supervisor; restarting a trial resets persistent/controller states.
Changing weights requires a fresh model build. Both64/128 packets have fixed128
storage; the padded64-cell output exactly matches its original model.

## State, limits and diagnostics

Each inverse candidate uses one immutable pre-current GRU state. The final
selected, post-slew pressure commits one neural update per new valid sample.
Short invalid gaps hold pressure and memory; no interpolated/imagined measurements
are inserted. A recovery interval<=.10s advances once with that actual dt and
sets historyGapHeld. This explicit bounded-gap deployment policy differs from
the strict gap-reset offline training/evaluation segmentation; hardware/model
accuracy through gaps remains unqualified. Longer gaps reset/re-warm memory;
an already active controller follows invalid grace/fault rather than using cold
history. A repeated timestamp cannot advance neural state or issue new model-on
feedback; physical measurement freshness is part of the contract.

Status uses PAC's codes0off/1trim-warmup/2blend/3control/4invalid-hold/5fault.
Faults slew toward zero and require a fresh run; no automatic rearm. Startup
and ramp-down pressure can be below the configured operational minimum by
design. Model-on inverse searches reuse the same fixed/variable-common allocator,
rails and trim-relative caps; every candidate includes the pressure slew.
It is a bounded LOCAL solve, not a global optimum guarantee. An unmet force
request sets limited and suppresses PAC adaptation via the existing gate.

Log all original signals plus ready,advanced,historyGapHeld,contractValid,
rawForce_N,effectiveResidual_N,gravityReexpression_N,inverseError_N,iterations,
converged,warm_s,units and parameterCorrection. ready denotes warmed history,
not force/model accuracy or hardware validity. modelValid retains the nominal
mechanical-input check; augmented usability also requires the diagnostic
contract/ready and supervisor status. In invalid periods the augmented
equivalent-force estimate is NaN rather than inferred from held measurements.
Interpret raw/effective residual samples with advanced,modelValid and the
contract/status flags: unevaluated diagnostic fields retain initialization
values and are not fresh force measurements.

Requested/component forces remain pre-neural control demand. In model-on mode,
tauEquivalent_N/allocationResidual_N and trim forceBias_N use the AUGMENTED
nominal force definition. None are measured physical force error. No added
neural low-pass filter, gain retuning, payload estimator or online optimizer.

The pressure sent to NI must match the selected command used to update memory.
Keep the external physical Kill path; ensure enablePressure reflects its permit.
Additional downstream saturation/offset/kill overrides can invalidate neural
history and must be resolved/logged before physical comparison. Never fake a
fresh physical timestamp by simply connecting Clock to measurementTime_s.

## Offline verification and remaining work

Canonical current results: validation_20261103_v4/verification.txt and .mat.
Cells exactly match frozen models; analytic pressure derivative agrees with
finite differences<1e-10 relative; activation/gravity/parameter definitions
checked; all28 PAC model-off outputs exactly match;48 known-target inverses
and804 complete synthetic-kinematic wrapper ticks pass bounds/slew/fault tests.
These are numerical/component tests, not physical closed-loop tracking.

MATLAB Coder emits C for both wrappers, then its standalone build stage returns
the project's previously seen toolchain error. The emitted C was independently
compiled using local GCC and replayed through2576 ticks,16 cases covering both
controllers,64/128, model-off/on and fixed/variable common. See
c_replay_20261103_v1/native_results.json. Maximum numeric difference<6e-12 and
pressure difference<3e-14bar. Strict legacy allAboveActivation can differ exactly
at the.8bar boundary due to roundoff; those explicitly enumerated logging-only
flag differences are retained in the raw traces/results. No logic was changed
to conceal them. This is standalone C evidence, not successful compilation or
execution of the user's Simulink/Desktop Real-Time model. Actual model build,
full acquisition/control timing and operator smoke validation remain pending.
