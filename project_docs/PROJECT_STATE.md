Last Updated: 2026-10-07
Current Milestone: 1 rad circle undertracking with lower-rail and symmetric-correction limits

# Current project state

Read this document first, then the specialized document needed for the task.
Check it against newer relevant recordings, saved model contents, and project
activity before relying on it. Historical statements do not override this state.

## Project objective

Develop and experimentally validate dynamics and feedback control for a single
hanging, backbone-constrained pneumatic continuum-arm module using paired NDI
orientation/position measurements and three pressure-regulator commands.

## Current architecture / model

- Two independent length-change coordinates `q=[l2;l3]` in metres;
  `l1=-l2-l3`. NDI lengths are inferred from calibrated orientation and
  constant-curvature geometry, not independent muscle strain.
- Active length 0.17411 m, PMA radial offset 0.013 m, mass 0.1 kg;
  pressure-area radius 0.0065 m. The ideal generalized force is
  `13.2732*[P2-P1;P3-P1]` N for pressures in bar.
- The bounded common-pressure inverse uses `K=650*[2 1;1 2]` N/m,
  `D=40*eye(2)` N s/m, downward arm-frame gravity, and the smooth
  +/-35 mm bound penalty. These are model assumptions, not independently
  identified physical parameters. The offline animation runner currently uses
  `D=20*eye(2)` and rotated world gravity; do not assume all models match.
- NDI: COM12, 921600 baud, paired sensors 0A base / 0B tip, nominal 0.05 s
  polls. The recorded calibration is
  `ndi_orientation_recordings/calibration_20260930_190505_898.mat`;
  verify mounting and preflight before physical use.

See [MODEL_ARCHITECTURE.md](MODEL_ARCHITECTURE.md) for interfaces and data flow.

## Current implementation status

The earlier small-step experiment is `FBL_commonP.slx` with paste-ready source
files in `common_pressure_testing/`. Its 240 s test ramps to trial-4 trim
`[1.6656,1.6492,1.4852]` bar at mean 1.6 bar, captures the static model-force
anchor at 50 s, blends incremental feedback over five seconds, and requests
separate 0.01 rad U/V steps. Corrections in both differential settings are
bounded to +/-0.04 bar. The normal command slew limit is 0.8 bar/s; V2 holds
a brief invalid input and ramps to zero after a sustained fault. No additional
0.8 bar deadzone offset is added to this all-active allocation.

Read-only SLX archive inspection after the V2 recording found
`commonPressureFeedbackInverseBlockV2` in the main `MATLAB Function1` chart;
`MATLAB Function8` retains a V1 copy. The model was saved before the V2 MAT
recording. The chart historically named thetaCircle contains
`commonPressureFeedbackReferenceBlock`. Full wiring, Kill state, and
post-Kill pressure behavior have not been verified from this recording.

The current circle experiment is the user-maintained `FBL_commonP_circle.slx`.
Its saved reference now has the requested one-input `thetaCircle(t)` with
internal `des_phi=1`, `circle_time_s=30`, `num_circles=5`, and circle start
10 s. The saved inverse uses `maxDeltaS_bar=maxDeltaA_bar=1`, 1.2–3 bar
active pressure bounds, a 2 bar/s command slew, active trim/feedback from
t=0, a five-second blend, and venting from 220 s. The latest circle MAT
logs Kp=[3000;3000], Kd=[120;120]. The saved SLX was written after the
trial-4 MAT; use the MAT logs for actual run settings and do not assume every
earlier trial used the latest saved reference.

Earlier models and workflows remain available: nominal feedback linearization,
trackingWorks/circle trials, stiffness-function variants, inverse analysis,
P1 stiffness fits, and offline stiffness sweeps. Their settings differ by model
and recording. No controller/model/code was changed by this migration.

## Validated components

- Paired NDI calibration/preflight and explicit measurement-validity handling
  have recorded offline/measurement checks.
- Causal velocity-estimator comparison supports the 0.20 s filtered-difference
  baseline, without establishing true velocity accuracy.
- Early feedback-linearization functions passed synthetic checks and an early
  hardware-free model copy passed a 0.5 s smoke simulation. This does not
  validate subsequent saved model revisions.
- User-recorded hardware experiments now exist; the October 1 statement that
  no closed-loop hardware test had occurred is historical.
- Common-pressure source functions passed documented offline numerical checks;
  an earlier physical recording exposed a fault-path exception absent from the
  nominal slew checks, and V2 subsequently held one invalid frame. The latest
  circle runs expose bend undertracking. Results and limitations are in
  [VALIDATION_RESULTS.md](VALIDATION_RESULTS.md).

## Latest circle result and current task

The latest inspected recording,
`FB_circle_phi_1rad_trial6_good_5s_0.5bar.mat`, logs five 5 s circles at
desired phi=1 rad with Kp=[3000,3000], Kd=[22,22] and a user-reported
0.5 bar active command floor. Mean measured phi=0.9883 rad, phi MAE=0.0485
rad, and measured theta advances 4.989 turns. Pressure limiting falls to
65/500 circle samples; 93/500 have at least one command below 0.8 bar,
but no sample reaches 0.5 bar. Theta MAE remains 0.3173 rad and l2/l3
RMSE=2.7649/3.0615 mm. Circle time and Kd also changed, so the floor's
independent effect is unmeasured. Eleven circle samples are NDI/model
invalid; NI-path and controller software commands match exactly.
An additional retrospective comparison estimates a repeatable 0.23–0.25 s
theta phase lag across the five laps. Shifting desired theta and lengths
back by 0.24 s in offline comparison reduces theta MAE to 0.069 rad and
l2/l3 RMSE to 0.82/1.11 mm; this is not a controller improvement or a
uniquely identified actuator delay. The 2 bar/s command slew is reached
in 65/499 adjacent circle intervals. Keep the current gains/floor as the
baseline while checking repeatability and speed dependence before tuning.
The newer 10 s/circle MAT at the same logged phi target=1 rad and gains
has phi MAE=0.0485 rad, theta MAE=0.2713 rad and only 5/1000 limited
samples. Its retrospective pure time shift is 0.40 s, showing that a
single fixed 0.24 s physical delay is not established. Across the 5 s
and 10 s runs, theta error fits an offset about 0.15–0.20 rad plus a
speed-related coefficient about 0.12–0.13 s. A 0.12 s theta-reference
advance is a proposed isolated test, not a verified compensation; a
remaining angular offset must be evaluated separately. The saved SLX
currently sets `des_phi=1.5`, whereas both inspected MATs logged 1.0.
The user subsequently clarified that the desired implementation is a
compensator added to the existing FBL/PD controller, not a manual change
to the `thetaCircle` timing. The reference-advance code shown in discussion
was not adopted in the model. A bounded desired-velocity feedforward or a
model-based predictor before pressure allocation is under consideration;
its gains and the separate speed-independent theta offset are unvalidated.

The earlier `FB_circle_phi_1rad_trial4.mat` logs five 30 s circles at desired phi=1 rad.
Measured theta advances 4.977 turns, but mean phi is 0.7425 rad, falling
from 0.7540 in lap 1 to 0.7360 in lap 5. The 1.2 bar lower rail is reached
on 2881/3000 circle samples; `pressureLimited` is true on 2894/3000, and
an offline reconstruction of the saved allocation gives a +/-1 bar symmetric
correction clamp on 1835/2990 valid samples. No channel reaches the 3 bar
upper rail. Mean actual phi varies by desired-theta sector from about 0.57
to 0.93 rad. Ten invalid NDI frames are isolated and the logged pressure
command holds through them. Trial 3 (one 25 s circle) likewise averages
0.7693 rad phi and is limited for 468/500 samples. The reported 5–10 s
"pre-circle hold" is actually a phi ramp, so neither run includes a steady
1 rad fixed-theta hold. Controller and NI-path pressure-equivalent logs match;
actual chamber pressure remains unmeasured. Earlier 0.2/0.1 rad trials had
undertracking without active limiting; the file named `FB_circle_phi_1rad.mat`
historically logged only 0.1 rad. See
[VALIDATION_RESULTS.md](VALIDATION_RESULTS.md) and
[OPEN_ISSUES.md](OPEN_ISSUES.md) for evidence and next diagnostics.
For subsequent circle MAT files, run
`common_pressure_testing/analyzeCommonPressureCircle.m` and share its compact
printed report; it reads saved logs only and detects the actual reference and
lap count. It includes logged Kp/Kd, validity, force and pressure summaries.
Usage is in the [experiment README](../common_pressure_testing/README.md).

The preceding small-step baseline remains useful:

`common_pressure_testing/feedback_testing_results_V2_newGains.mat` spans
240 s with Kp=[1400;1400], Kd=[120;120] and 0.01 rad U/V steps. At its sole
invalid NDI frame (22 s), V2 held the controller and NI-path commands and
resumed without a step; sustained-fault behavior remains an offline check.
Zero-target length-error norms were about 0.020–0.028 mm, while U/V holds
were 0.090/0.134 mm. The bending-plane angle was unobservable during
3298/3300 feedback samples, limiting directional conclusions from inferred
lengths near straight. Full figures and the earlier V1 fault history are in
[VALIDATION_RESULTS.md](VALIDATION_RESULTS.md).

## Open issues

1. V2 single-invalid-frame command hold is recorded; sustained-invalid fault
   latch and actual chamber pressure remain unverified physically.
2. Circle bend undertracking, especially l3, plus residual small-step error,
   trim drift and history dependence.
3. Visually observed local buckle is outside the tip-only/two-coordinate model;
   distributed shape stability is unverified.
4. Saved-model reference reproducibility: the nominal 1 rad file contains a
   0.1 rad target, and the saved reference still hardcodes that value. Chamber
   pressure remains unknown; the circle MAT files include controller and
   post-Kill NI command logs.
5. Parameter consistency across models and the pending historical Inv_Dyn_val
   +/-30 versus +/-35 mm bounds discrepancy.
6. Earlier tapped-delay predictor lacks simultaneous-actuation recursive rollout
   validation for NMPC.
7. Bending-plane angle is unobservable near straight; orientation-derived
   lengths use the held angle, limiting directional tracking assessment.

Details, attempted approaches, and next actions: [OPEN_ISSUES.md](OPEN_ISSUES.md).

## Important current assumptions and collaboration rules

- The user chooses to treat command pressure as chamber pressure for analysis
  because chamber sensing is unavailable. Prefer post-limiter/post-Kill logs
  where available; retain this as an assumption, not a measured fact.
- A stiff backbone constrains the length sum; common pressure is an ideal
  bending-force nullspace, but may alter prestress and local response.
  Tip Z change alone does not prove axial shortening.
- Trial-4 trim is a local candidate, not a globally straight or buckle-free pose.
  Geometric U=V=0 and a measured zero-pressure reference are different targets.
- Discussion is the default. Questions and result inspections do not authorize
  edits to code/models or hardware operation. Maintain the relevant specialized
  documentation when recording material findings in the project log; this does
  not authorize controller/Simulink or hardware changes.

## Next steps

For the 1 rad test, separate lower-pressure clipping, +/-1 bar symmetric
correction clamping, and direction-dependent response before changing gains.
An offline allocation review can test whether a higher common command would
preserve requested differential pressure within 1.2–3 bar; this is not a
validated physical operating point because common pressure changes prestress
and the known local buckle. A measured fixed-theta 1 rad hold, if later
chosen, would distinguish static bend capability from circle-motion lag.
Retain the earlier V2 one-frame fault result, and verify its sustained-fault
path separately. Assess tip-only buckle limits and near-straight theta
observability when comparing trials. Consider bounded integral/disturbance/trim
correction only after evidence and explicit implementation direction. NMPC
remains a later candidate pending held-out recursive predictor validation.

## Documentation maintenance

After significant work, append a concise dated entry to
[PROJECT_LOG.md](PROJECT_LOG.md), update this state and its date/milestone, and
update relevant issues, decisions, results, architecture or technical notes.
Update for material changes, not file touches. This state is edited in place;
the log is append-only after migration. Preserve reversed decisions and their
reasons; remove resolved issues from the active list and record their resolution.
For routine work read this state plus the relevant specialized file. Search
history by date, function, experiment or term instead of routinely loading the
entire log. Local experiment READMEs retain operational instructions and source
result artifacts retain detailed provenance.
