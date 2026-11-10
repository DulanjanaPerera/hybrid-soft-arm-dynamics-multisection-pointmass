# First PAC residual trial: offline fault diagnosis

Reviewed 2026-11-03. Input: ../firstTime.mat. No hardware/model execution,
controller-source changes, SLX changes, retraining or protection changes during
this review. Numeric helper replay uses the logged measurements and demands;
it does not simulate or establish the physical closed-loop response.

## Finding

Pressure shutdown starts at **20.80 s** because two missing tip-orientation
samples trigger neural-history reset on measurement recovery. There is no
tracking-error-magnitude trip in this supervisor. The inverse was converging
immediately before the fault. This recording does not establish physical model
accuracy, nor identify physical instability as the shutdown cause.

| Simulation time (s) | Orientation valid | Neural ready | Status | Action |
|---|---|---|---|---|
| 20.65 | 1 | 1 | 3 | Last valid control tick |
| 20.70 | 0 | 1 | 4 | Hold pressure |
| 20.75 | 0 | 1 | 4 | Hold pressure |
| 20.80 | 1 | 0 | 5 | Reset neural memory; latch fault; slew pressure to zero |
| 21.80 | 1 | 1 | 5 | Memory warm again, but fault remains latched |

At both invalid ticks sensorValid=[1,0], base orientation remains valid and
overrun=0. This narrows the observation to tip acquisition/validity loss, but the
specific underlying NDI rejection reason is not available in this recording.

Last accepted measurement time changes from20.65 to20.80: gap0.15s exceeds
residualOptions(5)=0.10s. residualPressureStep resets h and warm_s. An already
anchored controller requires ready; the recovered sample therefore remains
unusable during the1s warmup. The time since last usable observation exceeds
supervisorOptions(4)=0.10s, latching status5. This is an operational limitation
of the supplied recovery policy. Prior synthetic tests covered one missing tick
and three missing ticks, but missed this two-missing-then-recover boundary.

The logged simulation continues until96.7s and ends with StopCommand; its saved
ErrorDiagnostic is empty. Distinguish controller pressure shutdown at20.80s
from the later simulation stop.

## Recorded settings and numerical evidence

- Residual enabled throughout; this was not a mid-run residual switch.
- GRU128; inverse6iterations maximum, tolerance0.02N; warmup1s; history gap0.10s.
- Common-pressure adjustment enabled throughout. Controlled pre-fault mean
  pressure1.6000–1.7271bar, nominal1.6bar. Activation0.2bar, rails0.2–3bar,
  command slew6bar/s. Physical Kill signal remains1.
- Pre-fault blend/control inverse error norm maximum0.0197N. Immediately before
  dropout, error components are approximately[-1.30e-5,1.94e-5]N. Numerical
  inversion of a learned model does not measure true force error.
- Final logged NI-path pressure commands match pressureCommand_bar exactly.
  These are commands, not measured chamber pressures.
- Replaying residualPressureStep from initialization through first fault gives
  exact status/ready agreement; max pressure difference8.53e-14bar and raw neural
  force difference3.98e-13N. Replay uses filtered q/dq, logged gravity, demand,
  static+hysteresis+bias anchor input, actual options and accepted timestamps.
- An isolated original pacPressureStep replay seeded at20.65 gives statuses
  [4,4,3] at[20.70,20.75,20.80], versus residual[4,4,5]. This isolates the
  history-readiness policy difference; it is not a matched nominal arm run.
- Zero rawForce diagnostics during rejected samples are unevaluated defaults,
  not measured physical force collapses.

## Common pressure

Retain the common pressure/preload operating point. Variable-common allocation
is optional for the residual model: either fixed or variable common can be
used, but retain the same choice in matched model-off/on comparisons. Turning
it off does not resolve this history-reset fault. Neural inference still takes
all three final pressure commands as features.

## Timestamp and source provenance limitations

Recorded measurementTime_s is a validity-latched simulation clock. The saved
validMeasurementTime chart implements that policy. Host poll intervals have
median0.0508s, P950.0509s, maximum0.3089s over this recording; host elapsed time
exceeds simulation elapsed time by2.0011s. This is material to a model trained
using physical elapsed time. The source of dt used by inference needs an
explicit clock contract; do not feed unaligned host timestamps directly into
the current simulation-clock freshness check.

The saved PAC_residual.slx contains the old23-input PAC chart, whereas the
recording includes residual diagnostics and a residual switch. Thus it is not
a saved snapshot of the executed edits. Save the actual edited model before
the next trial. The numeric replay independently establishes that current
residual helper behavior explains this log, but does not prove every runtime
model connection from the stale saved SLX.

## Recommended next sequence

1. Validate forward predictions in passive shadow mode while the existing PAC
   commands the arm. First replay saved nominal logs offline. Advance GRU
   memory once per accepted physical observation using the final NI-path
   commands, actual physical dt and causal measured q/dq. Preserve missing-data
   masks and distinguish raw inferred lengths from the controller's filtered
   input representation. Compare nominal, constant velocity and augmented
   one-step predictions against subsequent valid measurements; separately
   evaluate0.1–0.5s rollouts. Never score invalid/held observations as ground
   truth. A retrospective actual-next-dt audit is not prospective timing proof.
2. Evaluate the inverse in shadow only: same measured state/history and demand,
   pressure bounds/slew/common policy, with predicted pressure disconnected from
   the valve path. Log convergence, predicted achieved force, mismatch, latency
   and pressure changes. Known-model forward/inverse round trips establish
   internal consistency, not physical inverse accuracy. Physical validation
   needs subsequent measured arm response under operator-authorized trials;
   this apparatus has no measured chamber pressure or direct force ground truth.
3. Design/test explicit dropout recovery and a separate/aligned physical-time
   input before another residual-on comparison. Regress0/1/2/3+ missing samples,
   repeated timestamps, warmup, history reset, safe transition and latched fault.
   Keep protections in place; simply raising grace/gap limits or bypassing
   readiness does not validate cold-history predictions.

No new hardware trial is scheduled or started by this analysis.

## Artifacts

fault_sequence.png/.pdf/.fig: shutdown-window plot.
fault_window.csv: recorded signals around the event.
summary.json and audit.mat: settings, replay errors, clock statistics and replay.
source_snapshot/: analysis/helper versions and SHA256 inventory.
Reproduce from the project root with addpath('residual_controller');
analyzeFirstResidualPac(pwd, fullfile(pwd,'residual_controller','Experiments', ...
 'PAC','20261103','firstTime_review_v2')); use a fresh output directory.
