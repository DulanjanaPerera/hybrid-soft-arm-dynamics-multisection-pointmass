# Residual-on circle sweep test5: tracking, invalid timing and pressure dips

2026-11-09. Offline audit of circlesweep_FWD_INV_test_5_fwdPressure.mat.
Model1.59, reached375s stop time,7501 observer samples, no hardware/model run
by the agent. Only out is saved; actual option variables are absent. Existing
source and logged diagnostics establish the reset/warm-up gates; inferred
selector behavior is checked directly against selected/applied commands.

## Main finding

The user's explanation is supported for the **physical bending dips**: abrupt
switching from corrected to lower nominal commands deepens the phi loss.
However, the initiating **sensor-validity loss or physical-time gap occurs
before the fallback pressure step**. That later step cannot explain the
initial invalid sample. The observed chain is:

```
NDI tip MISSING or physical dt > .1s
 -> observer breaks its segment / resets velocity and GRU history
 -> inverse context not ready; corrected pressure becomes NaN
 -> selector chooses original nominal command
 -> Unit Delay applies it .05s later
 -> bending falls during the long readiness warm-up
 -> corrected command returns once the observer is ready
```

The saved model still wires the validity helper's valid input to A5/base
orientation, not the new diagnostic flag. Here base validity is always1 and
the helper's finite-pressure test reproduces diagnostic readiness exactly.
Thus this wrong source should still be corrected per wiring_review, but it
does not explain the pressure dips in this recording.

## Timing and NDI evidence

-49 orientation-invalid samples; base invalid0. Every invalid sample has tip
 sensorValid0 and NDI raw tip status2=MISSING. They are genuine tracker
 invalidations, not only theta-unobservable near a straight configuration.
-9 physical sample intervals exceed0.1s. Physical dt median/P95/max:
 .050075/.051009/.635727s. Total physical elapsed379.994s versus375sim seconds.
 This localizes irregular observation timing, not its unique CPU/driver cause.
-13 overrun flags, all with orientation valid; not synonymous with the49
 missing samples. No observer numerical-guard status7. Command age always0.
-37 ready-to-not-ready transitions:32 start with missing orientation;5 start
 with a >.1s physical interval while both sensors still report valid.
 Other misses/gaps occur during existing recovery episodes.
-Median/P95/max fallback episode duration1.65/2.3825/2.60sim seconds
 (1.664/2.459/2.619host seconds). The configured design has.6s velocity
 warm-up plus1s GRU warm-up; more interruptions restart that recovery.
-Observer ready6143/7501samples(81.9%). Nonready1358samples(18.1%), including
 startup.1253samples are status3/4 warm-up, much more than the49 missing frames.

Several missing events recur at nearly identical circle phase, e.g246.3,
256.3,266.3,276.3s. This suggests a position/motion-dependent tracker issue;
the MISSING status alone cannot establish field geometry, cable, interference
or another unique physical cause. Do not loosen orientation validity or
interpolate the recurrent history to hide these events.

## Command selection and pressure drops

Reconstructing the valid/finite Switch and one-tick delay matches the logged
selected_pressure **exactly**. Post-Kill/saturation NI command vs that signal
matches4.45e-16bar. Kill is1 throughout; it did not cause these drops.
These are commanded pressures, not measured chamber pressures.

All37 first fallback commands appear exactly.05s after loss of readiness.
32 transitions lower at least one pressure channel by>.1bar; maximum drop
1.513951bar.33 transitions occur during active circles,4 during other phases.
26 episodes during a stable requested phi are followed by a>.1rad fall in
valid measured phi within2s. This association and the ordering support the
command-handoff mechanism; they are not a controlled physical causality test.

Example at182.30s (desired phi1rad):

| Simulation time | Logged event |
|---|---|
|182.30s|Tip MISSING; observer becomes invalid|
|182.35s|P1 command drops2.5164 ->1.0024bar;P3 drops1.9303 ->1.0453bar|
|Recovery interval|Fresh measured phi falls from about1rad to about.22rad|
|184.70s|Observer ready again after a further interruption during recovery|
|184.75s|Corrected command resumes through the Unit Delay|

The existing slew cap protects only newly computed corrected candidates.
The original nominal fallback is selected **after** that cap, so an abrupt
source switch can bypass it. This is why the1.51bar step can appear despite
the corrected branch's nominal .30bar-per-.05s allowance.

## Tracking progress and remaining limits

First three active circle reference profiles are identical to test4; a separate
comparison of all common desired-length samples gives max difference0m.
Both recordings use the same requested circle schedule over this interval.
These are separate runs, not randomized matched repeats with all physical
conditions controlled. Errors use valid measured samples, no NDI interpolation.

| Requested phi | Nominal test4 l2/l3 RMSE(mm) | Residual-on test5 l2/l3 RMSE(mm) |
|---|---|---|
|.2rad|1.812 /1.799|1.275 /.917|
|.8rad|6.961 /6.975|4.328 /3.220|
|1rad|8.203 /8.173|3.875 /3.923|

For phi1rad, within test5, phi MAE is.157rad when residual-ready versus.464rad
during fallback/warm-up; these subsets are not independent experiments and
are affected by preceding pressure/state history. All valid samples give
phi MAE.238rad,theta MAE.365rad. Length errors use actual orientation-derived
q=[l2;l3]; theta error wraps to[-pi,pi] and excludes theta-unobservable samples.

At phi2.5rad, lengthRMSE remains[11.988,11.285]mm andphiMAE.607rad. The full
run has1012/6143valid corrected requests reaching max-pressure scaling and
122/6143 reaching corrected-branch slew limiting. Higher amplitude performance
is therefore not solved by eliminating fallback steps alone.

## Recommended next code change (not implemented in this audit)

Replace the raw source Switch with a **stateful command-handoff supervisor**.
Every ordinary command change, including fallback and re-entry, must pass a
shared final slew limiter after source selection. Preserve physical Kill/stop
priority after that limiter; an emergency stop must not wait for smoothing.

For brief sensing losses, a bounded grace hold of the last actually applied
command can prevent the immediate step. If recovery is not promptly available,
blend/rate-limit toward the existing protected fallback, rather than holding
an old rotating-trajectory command indefinitely. Blend residual re-entry too.
Holding through the entire1.6s warm-up or retaining an old directional residual
while the reference rotates is a separate unvalidated policy, not a default.

This will reduce switching jerks; it does not recover a missing observation
or guarantee constant phi throughout a long re-warm interval. Keep the GRU's
readiness/reset policy intact while improving command handoff; investigate
the position-dependent tip misses and timing gaps separately. Simply increasing
gamma, changing NN weights, or raising timing/validity limits does not address
the demonstrated switching mechanism. No runtime source changed in this audit.

## Retained artifacts

-overview and selected event01/14/18 PNG/PDF/FIG, raw invalid gaps retained.
-fallback_events.csv: all37 entries, cause/timing/pressure/phi/recovery fields.
-trace.csv,analysis.mat,summary.json,status_counts.csv.
-tracking_by_trial.csv,separate_nominal_test4_metrics.csv.
-reference_comparison.json andsource/raw/model hashes/source snapshot.

Plot readiness flag means Switch selection **before** Unit Delay; actual
command selection changes one tick later. Case14 highlights a large phi loss;
case18 highlights the largest fallback pressure step. The event phiDrop metric
uses the median valid phi in the preceding.25s and minimum valid phi over2s,
not an instantaneous measured derivative or a guaranteed response time.

Initial v1 stopped on a Kill-signal alignment assertion and is markedincomplete.
Canonical v2 treats only the sparse event-logged Kill signal as a held control
value; sensor/observer data retain their original sampling/missing policies.

Reproduce with analyzeDesiredPressureRun(root,freshFolder). No SLX edit,
model/hardware/acquisition run,retraining,commit orpush was performed.
