Last Updated: 2026-10-06
Current Milestone: Bounded common-pressure feedback validation around the 1.6 bar trim

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

The active experiment is `FBL_commonP.slx` with paste-ready source files in
`common_pressure_testing/`. Its 240 s test ramps to trial-4 trim
`[1.6656,1.6492,1.4852]` bar at mean 1.6 bar, captures the static model-force
anchor at 50 s, blends incremental feedback over five seconds, and requests
separate 0.01 rad U/V steps. Corrections in both differential settings are
bounded to +/-0.04 bar. The normal command slew limit is 0.8 bar/s, with
an invalid-input exception described below. No additional 0.8 bar deadzone
offset is added to this all-active allocation.

During this documentation migration, read-only SLX inspection confirmed that
the saved inverse chart contains `commonPressureFeedbackInverseBlock` and
the chart historically named thetaCircle contains
`commonPressureFeedbackReferenceBlock`. This supersedes the earlier October 6
log observation that the saved file still held the old scripts. Full wiring,
current gain values, Kill state, and compilation were not verified here;
equivalence to the recorded run remains to be established.

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
  the latest physical recording exposed a fault-path exception absent from the
  nominal slew checks. Results and limitations are in
  [VALIDATION_RESULTS.md](VALIDATION_RESULTS.md).

## Latest feedback result and current task

`common_pressure_testing/feedback_testing_results.mat` spans 240 s with
4800/4801 valid NDI samples. Recorded gains were `Kp=[1100;1100]` and
`Kd=[100;100]`, constant during that run. U and V steps moved in the intended
directions but left persistent error: desired/measured U was 0.0100/0.0054 rad
in the U hold, and desired/measured V was 0.0100/-0.0032 rad in the V hold.
The local correction clamp was inactive throughout 3300 feedback samples.
The run does not establish a suitable gain increase or integral/bias design.

At 14.40 s a single invalid NDI pair drove pressure immediately to zero
(about 33.3 bar/s), then pressure recovered at 0.8 bar/s. The source's early
return resets pressure state and bypasses the normal slew limiter.
The next technical priority is to discuss and explicitly authorize a fault-path
fix, preserving an appropriate persistent-invalid fail-to-zero behavior, then
verify fault injection and the exact saved model before another physical run.

## Open issues

1. Invalid-input slew bypass remains in source and saved inverse script.
2. Residual small-step tracking error, trim drift, and history dependence.
3. Visually observed local buckle is outside the tip-only/two-coordinate model;
   distributed shape stability is unverified.
4. Saved-model wiring/reproducibility and missing post-Kill logs in the latest MAT.
5. Parameter consistency across models and the pending historical Inv_Dyn_val
   +/-30 versus +/-35 mm bounds discrepancy.
6. Earlier tapped-delay predictor lacks simultaneous-actuation recursive rollout
   validation for NMPC.

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
  edits to code/models or hardware operation. The October 5 standing automatic
  update authorization named the project log only; this request authorizes this
  documentation restructuring and supplies maintenance rules, not code changes.

## Next steps

Resolve the invalid-sample command behavior; establish exact saved-model wiring,
gain/settings, reference, validity logging and Kill path; validate offline fault
cases before a user-directed experiment. Then compare repeatable long holds and
small reference steps, assessing drift, pressure authority and local buckle.
Consider bounded integral/disturbance/trim correction only after evidence and
explicit implementation direction. NMPC remains a later candidate pending
held-out recursive predictor validation.

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
