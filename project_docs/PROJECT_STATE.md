Last Updated: 2026-10-12
Current Milestone: Adaptive K/D/history-force/bias path implemented and checked offline; new Simulink integration and physical validation pending

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

On branch `codex/passivity-adaptive-controller`, isolated stage-1 force code is
in `passivity_adaptive_control/`. `pacFixedBlock` accepts q/dq, desired q/dq/ddq,
arm-frame gravity, and validity flags. It implements fixed-parameter passivity
control with no integral, adaptation or payload. Numerical model-identity and
matched-model ODE checks passed. `pacCommonPressureBlock` now adds fixed-mean
allocation, a constant calibrated trim-force bias, slew/rails and a freshness/
fault supervisor. Command-path tests and MATLAB Coder C source generation
passed; duplicated-SLX wiring, closed-loop sampled pressure behavior and
hardware behavior remain unvalidated. Existing Simulink models are unchanged.
See that folder's README before use; the force output is not a valve command.

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

The latest inspected Ki comparison is
`StaticTest/changeKi_dynamically.mat`: Ki=[700,700] through 43.31 s,
then [0,0], while Kp=3000 and Kd=42 stay fixed. The reference is five
10 s circles at phi=1.5 rad, so the switch occurs during lap 3 rather
than during a static hold. Comparing complete on laps 20-40 s with off
laps 50-70 s, phi signed error worsens from +0.010 to +0.089 rad and
phi MAE from 0.0917 to 0.1080 rad after Ki is disabled. Wrapped theta
MAE stays about 0.22 rad and l2/l3 MAE about 2.4/3.0 mm. With Ki on,
phi error varies more and `pressureLimited` is true for 81% of command
samples versus 6% after Ki is disabled. The integral state is not logged.
This supports integral correction of systematic phi underbend but leaves
tracking tradeoffs and pressure-aware anti-windup unresolved. The saved
BRot model inspected by Codex still lacks the user's running Ki/gait edits;
use MAT logs to identify tested settings. Detailed results are in
[VALIDATION_RESULTS.md](VALIDATION_RESULTS.md).

In a later user-recorded `integral_term/plot_data_110s.mat` run, pressure
ceased around 123 s because tip sensor 0B stayed invalid for about 0.45 s,
tripping the inverse controller's >0.10 s latched-fault path. Commands
ramped to zero and did not resume when tracking recovered. The reference
continued; the saved model StopTime is 250 s and the active inverse/integral
time bounds are 220 s, with no inspected 120 s cutoff. See
[VALIDATION_RESULTS.md](VALIDATION_RESULTS.md) and [OPEN_ISSUES.md](OPEN_ISSUES.md).

A 10 cm X-Y square reference was proposed with Z following fixed-backbone
geometry, a 15-25 s corner approach, 25-30 s hold, and laps beginning at
30 s. The user has since saved `square/initial_test.mat`, which contains
logged `des_position` and `meas_position` XYZ. The new
`square/plotSquareTracking.m` opens the six tracking figures and uses those
logged signals directly in its 3D figure. Its offline plotting check passed;
square tracking performance has not yet been evaluated. Codex did not edit
or run the Simulink model.

The earlier circle baseline remains relevant:

An earlier inspected recording,
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

1. V2 brief-invalid hold and sustained-invalid fault latch are recorded in
   command-path logs; tip 0B intermittently dropped out, and actual chamber
   pressure remains unmeasured.
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

For PAC, wire `pacCommonPressureBlock` in a user-created duplicate only after
checking the last-successful-measurement timestamp, discrete 0.05 s rate,
gravity convention and downstream Kill/NI path. Keep the reference held until
status=3; calibration can be delayed by motion. Preserve the 0.2-3 bar saved
rails for now, while logging nominal activation crossing and force residual.
The current command floor does not guarantee all-active force mapping.
There is no payload, adaptation or integral in this new path. See
`passivity_adaptive_control/README.md` for the full interface and timing.

For the current integral controller, compare repeated matched circle or
fixed-pose trials at lower and higher Ki, log the integral state, and separate
the cause of `pressureLimited` before increasing gain. Check reversal and
zero-reference recovery; implement pressure-aware anti-windup if saturation
persists. For the earlier 1 rad test, separate lower-pressure clipping,
+/-1 bar symmetric correction clamping, and direction-dependent response.
An offline allocation review can test whether a higher common command would
preserve requested differential pressure within 1.2–3 bar; this is not a
validated physical operating point because common pressure changes prestress
and the known local buckle. A measured fixed-theta 1 rad hold, if later
chosen, would distinguish static bend capability from circle-motion lag.
Retain the earlier V2 one-frame fault result, and verify its sustained-fault
path separately. Assess tip-only buckle limits and near-straight theta
observability when comparing trials. NMPC
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

First PAC circle recording reviewed: engagement and NI-path commands worked,
with no circle pressure limiting; bending remained near zero. Nominal offline
feedback gains provide only 20 N/m restoring authority and require experimental
tuning. See VALIDATION_RESULTS.md (first PAC circle recording).

Static PAC test confirms weak feedback with gains Ks=10, lambda=2; reusable
report/plots are in PAC_circle/analyzePacStatic.m (see local README).

PAC experimental baseline now lambda=200 /s, Ks=20 N s/m (reported trials).
Static offsets depend on direction: phi MAE .04786/.00935/.03287 rad at
theta=0/pi/2/pi; pi/2 theta MAE .16734 rad despite good phi. No hold limiting.
Adaptation remains pending; offset versus stiffness/coupling is not identified.

Ten-hold sweep shows repeat theta0/.2 actual phi drops .14517 to .10945
with same gains/trim bias. History/time effects confound static stiffness fit;
repeatability and zero-reference baselines are next before parameter selection.

Adaptive implementation in passivity_adaptive_control/ADAPTIVE_INTEGRATION.md.
Structured stiffness, damping, history force weights and bias adapt; mass/gravity
and Bouc-Wen internal coefficients stay fixed. Existing nominal files and SLX
unchanged. New Dhat*v compensation changes baseline; compare learning on/off
using same wrapper. No sampled physical stability/performance guarantee.

First adaptive-wrapper physical baseline oscillates even learning disabled:
45-60s phi ~.43rad for desired.2, ~90% limiting. lambda200/Ks20 carryover
is unsuitable for new Dhat*v path. Retune fixed baseline before learning;
proposed lambda50/Ks40 starting comparison is not yet experimentally validated.

Reduced-gain lambda50/Ks40 adaptive wrapper now shows no static-hold limiting.
Learning reduces phi MAE .04794 to .03342 rad in separate45-60s runs;
not fully converged and repeatability/other directions remain unvalidated.

Adaptive sweep .3rad accurate across three directions (phiMAE .0033-.0050),
but final repeatedtheta0/.2 degrades .1802 to .1364 actual. Transfer/readaptation
remains open. .5rad application test proposed with unchanged gains/rates.

High-phi sweep (.5/.8/1.5) oscillation begins on .8 unloading110-115 s.
Slew limiting starts113.95 before invalid orientation115.70 and activation
crossing128.50. Zero-reference115-125 phi .0186-.462, slew86.5%; k3 drops52%
during unloading but K stays SPD. Delayed force/pressure loop plus rapid transient
learning is leading hypothesis, not proven cause. Mean-hold reports conceal
oscillation. See VALIDATION_RESULTS.md high-phi onset analysis. No code changed.
Repeat high-phi trial2 confirms .8 unloading failure: first limiting113.80s vs113.95s; zero-hold phi .016-.459,88% limited. See validation repeat entry; causal mechanism still unresolved.
Freeze comparison: trial2 parameters stop changing102.91s yet unload oscillation persists; ongoing updates not necessary. Trial1 also adapted, not off-baseline. New .01 logs with .05 measurement timestamps warrant multirate execution check before tuning.
Latest MATLAB restart clears NDI slowdown. Frozen trial now quiet after .8 unloading, residualphi~.14-.22; major limiting/oscillation starts larger loading ramp~128s near activation crossing128.15. See latest validation; mechanism not uniformly unloading.
Latest1.0 target trial oscillates already on .8 return115-125, all pressures above activation; ongoing adaptation/deadzone crossing not necessary. Next controlled timing/gain comparison; see validation entry.
Trial3 logs .05s,lambda50/Ks40: unloading oscillation persists85% limited115-125; timing alone insufficient. Ks20 comparison still pending.
Ks20/lambda50 trial4 suppresses unloading oscillation:115-125SD.0119,no limiting;1radMAE.04165. Provisional baseline, residual zero bend~.2 and pi2 theta error.511 remain; see validation.
Full adaptive Ks20 sweep tracks all holds well, recovered oscillation265-280 remains. Next target starts275, so recovery not evidence of self-decay at zero. k3 hits projection0; see validation. Retainlambda50/Ks20.
User repeat/20s zero hold remains oscillatory. k3=0 does not imply singular K; negative k3 requires coupled SPD projection and altered physical interpretation. No bounds changed.
Lambda25/Ks20 removes observed oscillation but slow directional readaptation, final repeatphiMAE.1708. New1s ramps confound lambda-only comparison; retain conservative baseline and matched-reference tuning.
Lambda40/Ks20 shows improved compromise vs25 and50: brief return transient, allholdsnolimits,finalrepeatMAE.0183 vs25.1708. Provisional baseline pending repeat; see validation.
Adaptive1rad30s circles:phiMAE~.03,thetaMAE.235->.195 over7laps,lengthRMSE~2mm; k2hits50,activationresidual gateslearning. See circle validation; retain speed for phase diagnostics.
Repeat10circles1rad30s plateaus thetaMAE~.28rad;activationresidual gates~2/3lap3 despite nolimits. Allocation consistency priority before gains; see validation.
Activation-consistent fixedmean allocator implemented; exact feasibleforce mapping, constrainedfallbackflags residual. Pressurebounds/gainsunchanged. Repeat30scircle; see architecture/validation.
Newfixedmeanactivationallocator didnotimprovecircleerrors; exactcandidatefrequently violateslow/highrails. Variablecommon bounded allocation nextproposal, not implemented.
