# Circle/helix inverse and NDI validity audit — 2026-11-05

The similar pressure outputs have two demonstrated causes: the hypothetical
inverse is clipped around the actual command every tick, and the projected
local inverse solver can stop on flat regions created by that clipping.
The learned force is not zero. The smooth displayed XYZ traces hide genuine
tip-tracker losses by holding the preceding valid coordinates.

Only saved data and numerical helpers were executed. No SLX, NDI, controller,
network, hardware, acquisition or pressure settings changed. No training,
commit or push. The newly added MATLAB files are offline analysis helpers.

## Recordings and counts

Both operator recordings contain `out` from `residual_FWD_INV_validation`,
version 1.51, normal mode, with ReachedStopTime and no saved simulation error.
All counts below use the common 0.05 s logged observation grid; near-straight
theta counts are samples, not distinct events.

| Quantity | Circle sweep | Helix |
|---|---:|---:|
| Duration / rows | 375 s / 7501 | 290 s / 5801 |
| Base sample invalid | 0 | 0 |
| Tip sample invalid | 17 (0.23%) | 28 (0.48%) |
| Orientation invalid | 17 | 28 |
| Configuration invalid at NDI grid | 17 | 28 |
| Theta unobservable with valid orientation | 9 | 796 |
| Forward ready / scored | 6908 / 6873 | 5056 / 5030 |
| Inverse finite/valid | 6908 | 5056 |
| Nominal and augmented pressures equal | 5510 (79.76%) | 3546 (70.13%) |
| Inverse converged | 156 (2.26%) | 245 (4.85%) |
| Maximum pressure difference | 0.60 bar | 0.60 bar |

Equality tolerance is 1e-8 bar. Convergence means the model's numerical
generalized-force mismatch is at most 0.02 N; it is not measured physical
force accuracy. A finite inverse output (`valid`) does not imply convergence.

## NDI rejection stage: missing tip transforms, not angular-only rejection

All 17/28 orientation-invalid samples have raw tip sensor status 2 = MISSING,
driver-usable false and raw-position-finite false. The base remains valid.
There are no orientation failures with both sensor samples valid. Configuration
invalidity coincides with those same losses on the NDI grid, with no additional
XYZ-fit rejection while orientation is valid.

The logs contain no out-of-volume, partial-volume, repeated-frame or pair-frame
rejection at these events. None coincides with the logged overrun flag.
This identifies where data were rejected, but does not identify why the
tracker returned MISSING. Do not infer a unique speed, field, cable or sensor
cause, or conclude that real-time timing is fully qualified. The maximum poll
intervals were 0.3827/0.3103 s and there were 2/1 overruns elsewhere.

`sensorValid` contains base then tip. A scalar or first element remaining one
does not establish tip validity. Diagnostic reject bit 4096 appears on both
columns when the pair is invalid; it does not mean both sensors failed.

`configurationValid` belongs to the separate XYZ configuration fit.
`orientationValid` belongs to the quaternion-based estimator. Their agreement
here is explained by the upstream missing tip, not by their being identical
definitions.

## Why XYZ looks continuous

At every missing-tip sample, internal relative `position` and masked raw tip
coordinates (`raw5/raw6/raw7`) are nonfinite. Both `meas_position` and
`meas_positionTip_offsetCorrected` stay finite and exactly equal their previous
values: 17/17 and 28/28. The saved model contains downstream last-valid-position
functions, which explain these held values. The saved SLX was read as a ZIP/XML
file only; it was not loaded into Simulink or executed. Unsaved operator model
changes cannot be established from the saved SLX.

See each case's `raw_vs_held_xyz.png` (also PDF/FIG): raw position has a gap
where the public XYZ traces remain continuous. Corrected and uncorrected XYZ
have different reference offsets; their absolute values need not coincide.
These plots do not show that XYZ tracking survives the loss better than
orientation. Keep freshness flags alongside held display outputs. For an
honest measurement plot, mask held samples with NaN using tip validity; do not
feed held or interpolated points to the GRU as fresh observations.

## Theta observability is a different issue

The 796 helix and 9 circle valid-but-unobservable theta samples have bends
below the calibration's 0.5 degree straight-angle tolerance. Their phi ranges
are 0.000189–0.008720 rad (helix) and 0.000890–0.008286 rad (circle).
Near straight, the bending-plane direction becomes ill-defined; holding the
previous theta is intentional. This is not a missing-transform count.

The forward observer gates on orientation/base validity and timestamps, not
thetaObservable or configurationValid. It accepts every orientation-valid
sample here, including every near-straight theta-unobservable sample, and
accepts none of the orientation-invalid samples. Thus those 796 samples do
not cause the reported GRU resets. Actual missing samples do reset history;
the 0.6 s velocity plus 1.0 s GRU warmup amplifies their effect on readiness.
Predictions during invalid/unready periods are withheld, not qualified as
valid predictions. Do not loosen validity or readiness to conceal the gaps.

## Why pressure predictions overlap

The block targets the nominal EoM left side at **measured q/v** and the supplied
desired acceleration. It does not receive desired length or a tracking-error
feedback demand. It is a diagnostic inverse, not a full trajectory controller.
With the force residual on the RHS, its nonlinear solve seeks `B*p+r(q,v,p,h)
= target`. A nonzero residual can increase or decrease pressure differences;
it does not necessarily add positive pressure to every actuator.

The default hypothetical slew allowance is 6 bar/s times 0.05 s = 0.30 bar.
Both inverses are constrained around the **actual current command**, not their
own preceding hypothetical output, at every tick. The observed maximum step
is exactly 0.30 bar. The nominal output is limited at every valid tick, and
the augmented inverse is flagged limited at 6835/6908 and 4957/5056 ticks.
Distinct unconstrained demands can consequently collapse to the same pressure
edge. A shadow sequence cannot progressively walk to a remote solution under
this reanchoring rule.

Options were not saved as complete input vectors. Replaying 120 distributed
valid contexts per file with source defaults and variable common pressure
reproduces nominal/augmented pressures within 2.09e-9 bar and convergence flags
exactly. Fixed common pressure does not reproduce the outputs. This supports
the source-default explanation, without uniquely identifying every hidden
option. The learned current-force norm medians are 10.42 N and 8.94 N.

Offline copied-option probes with the slew restriction removed (rails and
correction caps retained) separate the two predictions by up to about 0.83 bar
in selected equal-output cases. At helix 39.20 s, logged outputs both equal
[2.4336,0.3000,1.5443] bar. Offline no-slew nominal is
[1.9423,0.5925,2.2652], augmented [2.6844,0.1502,1.9654] bar; the latter model
inverse error is 0.000217 N. This is a computation, not a recommendation to
remove the live slew protection and not physical pressure validation.

## Confirmed projected-solver limitation

Sixteen selected equal-output cases were probed using the identical nominal
allocator, common-pressure policy, rails, correction caps and actual-command
slew box as the logged solve. All had a zero local projection Jacobian;
15 had better candidates reachable through that same projection. The existing
local method stops after one iteration at these plateaus. At helix 39.20 s,
a candidate [2.4336,0,1.5443] bar reduces mismatch from 4.9088 to 2.2984 N
without changing the limits. This is a real limitation of the supplied solver,
in addition to genuinely constrained targets, not evidence the NN is disabled.
The coarse probes are not global optima and do not show all targets feasible.

The earlier 5^3 pressure-box probes in each `offline_constraint_counterfactuals`
table are broader feasibility checks; only `same_projector_plateau_probes.csv`
establishes improvement under the exact projector policy.

## Recommended next work

1. Keep genuine sensor validity and distinguish near-straight theta from tip
   loss. Display raw gaps/held-data flags. Use existing tracker diagnostics to
   compare the same path slower and inspect repeated affected poses/placement
   in operator-run trials; do not assume a cause from these recordings alone.
2. Improve and verify the **disconnected inverse** on these saved contexts:
   add a bounded plateau escape or suitable pressure-space search, log
   constrained versus diagnostic unconstrained solutions separately, retain
   immutable GRU history and pressure protections. No production fix has been
   applied in this review. The shared controller allocator is affected too.
3. Then evaluate pressure sensitivity/forward consequences and actual tracking
   with a defined outer feedback demand and explicit dropout recovery design.
   Current inverse curves and short-horizon overlays do not establish that a
   residual-on controller is ready.

## Reproduction and retained evidence

`analyzeInverseNdiValidation(root, freshFolder)` performs extraction, default
replay, validity diagnostics and copied-option counterfactuals.
`probeInverseProjectionPlateaus(root, folder)` adds identical-projector probes
and raw-versus-held figures to an existing audit. No model execution needed.

Retained: summary CSV/JSON, per-event CSVs, decoded schema/reject-bit tables,
immutable-context audit MATs, counterfactual tables, full-run inverse/validity
and missing-sample close-up PNG/PDF/FIG, source snapshot and SHA256 inventory.
