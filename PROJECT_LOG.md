# Project research log

## Collaboration preference

The user prefers discussion by default. Do not modify code, Simulink models,
or other files unless the user explicitly requests the change. Questions,
discussion, and requests to inspect results do not authorize file edits.
This preference was explicitly recorded at the user's request on 2026-10-01.

This log records reproducible decisions, datasets, methods, results, and open
questions for the single-module NDI/dynamics/control experiments. Add dated
entries as work proceeds. Numerical results should remain linked to the raw
recordings and the code/settings that produced them.

## 2026-10-01 — Causal NDI length-velocity estimator comparison

**Question.** Before feedback-linearization control, compare a filtered
finite difference with a simple Kalman estimate of `dl2` and `dl3`. The older
six-state UKF in `Controlling/Dynamic_jointspace_feedbacklinearization` used
different dynamics, a 0.01 s step, and an angle measurement model; it was not
transplanted into the current controller.

**Repository state.** Work was performed on `codex/feedback-linearization-clean`
at HEAD `42fc8a0`, with existing user changes and recordings preserved. The
analysis code is isolated in `estimator_comparison/`; no Simulink or hardware
control files were edited, and no hardware was operated.

**Input data.** Six 30 s recordings, each with 601 NDI samples, from
`forwardModelValidation/`: `P1_sineWave.mat`, `P2_sineWave.mat`,
`P3_sineWave.mat`, `P1P2_sineWave.mat`, `P1P3_sineWave.mat`, and
`P2P3_sineWave.mat`. Logged commanded pressure spans 0–3 bar. The NDI
`lengthChnage_m` signal contains orientation-derived constant-curvature
length changes, not independent muscle strain. Ten of 3606 samples had
`orientationValid=false`; both estimators skip measurement updates for those
samples. Median `pollStart_s` spacing was 0.0502–0.0509 s by trial, with
occasional intervals up to 0.266 s. The median reply-completion (`hostTime`)
interval was about 0.057 s and includes BX latency and software overhead.
`pollStart_s` is used
as the practical sample timestamp; its relation to the internal NDI frame
time has not been independently calibrated.

**Data SHA-256.** These hashes identify the recordings used for this entry.

| Recording | SHA-256 |
| --- | --- |
| `P1_sineWave.mat` | `3FB571EA84C2A8A8E330A64E467CBF3EBFDA8914D9C54DEBB10827B80A368136` |
| `P2_sineWave.mat` | `55FA4FB00B18A3423780801DACA5A90EDDC41C94E11FB0D8D8727E783DF07D1B` |
| `P3_sineWave.mat` | `4BE75DDB50A862416ED04A073801E06E7DD2A959E4A99798EDC1F1BB26B04581` |
| `P1P2_sineWave.mat` | `8542BB769CF93F896E91A56E50717F314EFA388538DF54D3AFB3AB3C43B30B2B` |
| `P1P3_sineWave.mat` | `AE85F38004B9B52070D3DBBA439D99EB6490D2C8069684D083665D65F98B53EF` |
| `P2P3_sineWave.mat` | `D9F4A2A5F79E38C8D88400EE0D44D689855FA94FE9EBE026EED0680EECCFEFCE` |

**Method.** `estimator_comparison/compareNdiEstimators.m` performs two
causal estimates from `[l2,l3]`: a backward difference plus 0.12 s low-pass
filter, and a four-state constant-velocity Kalman filter. The selected Kalman
settings are 0.05 mm measurement standard deviation and 0.003 m/s² per-step
acceleration standard deviation. Those values were chosen after an exploratory
4×4 parameter sweep on the single-actuator runs; paired runs were also
inspected, so they are not a blinded holdout. A gap over 0.5 s resets the
velocity estimate. Results, settings, and six figures are saved under
`estimator_comparison/results/`.

**Initial 0.12 s observations.** The table reports velocity RMS (mm/s) in the first three
simulated seconds with commanded pressure at most 0.3 bar, followed by
one-step length-prediction RMS (mm) against the next NDI sample.

| Trial | Valid | Low-command filtered / KF (mm/s) | Prediction filtered / KF (mm) |
| --- | ---: | ---: | ---: |
| P1 | 601/601 | 0.043 / 0.042 | 0.087 / 0.082 |
| P2 | 600/601 | 0.071 / 0.068 | 0.106 / 0.097 |
| P3 | 598/601 | 0.064 / 0.060 | 0.111 / 0.121 |
| P1P2 | 600/601 | 0.075 / 0.073 | 0.119 / 0.123 |
| P1P3 | 596/601 | 0.077 / 0.073 | 0.088 / 0.115 |
| P2P3 | 601/601 | 0.049 / 0.048 | 0.126 / 0.120 |

The Kalman estimate is slightly quieter in the low-command segments, while
the filtered difference has slightly lower mean one-step prediction error
across all six runs (about 0.106 versus 0.110 mm). Neither metric provides
ground-truth velocity, and low commanded pressure does not prove the arm was
stationary. The figures show similar velocity trajectories with some sharp
NDI-derived fluctuations during motion. These results do not justify claiming
that one estimator is more accurate. The simpler filtered difference remains
the first-controller baseline; the Kalman filter is retained as a comparison
and a possible replacement if later tests show a useful noise/delay benefit.

**Next measurements.** Capture an intentionally motionless hold to estimate
measurement noise, then compare controller-relevant delay and stability at
the 0.05 s update rate. Keep NDI validity and pressure-command signals in all
closed-loop logs. Do not interpret filter smoothing as a correction for
calibration, deadzone, or model stiffness errors.

## 2026-10-01 — Filter time-constant sensitivity

The user observed visible oscillation in the filtered velocity traces and
asked whether stronger smoothing would be useful. Without changing the
analysis defaults or operating hardware, the same six recordings were
reprocessed at four causal low-pass time constants. The table averages each
metric across the six runs. Roughness is the RMS sample-to-sample change in
the estimated two-coordinate velocity vector; it is a smoothness measure,
not an accuracy measure.

| Filter time constant (s) | Low-command velocity RMS (mm/s) | One-step prediction RMS (mm) | Successive-sample velocity change RMS (mm/s) |
| ---: | ---: | ---: | ---: |
| 0.12 | 0.0633 | 0.1061 | 0.6135 |
| 0.20 | 0.0534 | 0.1027 | 0.4064 |
| 0.30 | 0.0464 | 0.1011 | 0.2868 |
| 0.50 | 0.0381 | 0.1009 | 0.1827 |

Increasing the time constant to 0.20–0.30 s substantially reduces visible
jitter and also improves this particular one-step prediction proxy. That
proxy does not measure controller phase lag or true velocity accuracy. A
0.20 s filter is a cautious starting point for the first controller test;
0.30 s is worth comparing offline and in low-gain closed-loop tests. The
comparison script default was changed from 0.12 to 0.20 s and its saved
figures, CSV, and MAT results were regenerated at 0.20 s; the prior table
above remains the historical 0.12 s baseline. A 0.50 s
filter should not be chosen merely for its smooth plot because its delay may
reduce derivative-feedback effectiveness. A sinusoidal pressure command need
not produce sinusoidal velocity when deadzone, stiffness variation, and
transients are present. No Simulink files or hardware were changed in this
sensitivity check.

## 2026-10-01 — First connected feedback-linearization model

Before editing, committed the user's new `feedback_linearization.slx`, updated
validation models, six pressure trials, and estimator-analysis results as
`ee6cd0e` so the starting point is recoverable. The new model had placeholder
derivative, error, controller, and inverse-dynamics functions; its controller
and inverse model were not connected to the valves.

`feedback_linearization.slx` now uses reduced lengths `[l2;l3]` throughout.
The desired length path is sampled at 0.05 s and uses causal first-order
filters with 0.10 s velocity and 0.15 s acceleration time constants. The
measured path takes NDI's inferred `[l1;l2;l3]`, selects `[l2;l3]`, and uses
the physical BX poll timestamp with 0.20 s velocity and 0.30 s acceleration
filters. Invalid readings hold the derivative state; the NDI orientation
validity flag is also connected to the inverse model so a missing reading
commands zero pressure. Errors are desired minus measured position, velocity,
and acceleration, computed directly from these filtered signals without
differentiating an error signal. The virtual acceleration is
`ddl_des + Kd.*(dl_des-dl_meas) + Kp.*(l_des-l_meas)`. Editable initial gains
are `Kp=[4;4] 1/s^2` and `Kd=[4;4] 1/s`.

The inverse uses the same one-section mass, Coriolis, downward gravity,
damping, nonlinear bound penalty, and coupled `1350*[2 1;1 2]` N/m stiffness
as the forward model. It maps generalized force through
`B = pi*(0.013/2)^2*1e5*[-1 1 0;-1 0 1]` N/bar, chooses the smallest
nonnegative three-pressure solution, applies the provisional 0.8 bar valve
deadzone to active channels, and proportionally scales effective pressure
when needed to keep all commands at or below 3 bar. Named logged signals
include the command, desired and available generalized forces, saturation,
and inverse-model validity. The valve Kill constant was saved as `0` (zero
pressure) because the existing theta/phi sine sources span 0–3 rad, much
larger than a reviewed first hardware target. Change the reference and
review a low-gain run before enabling pressure.

Offline MATLAB function checks passed for missing frames, sine-wave
derivatives, error sign, bounded pressure, force reconstruction, saturation,
and invalid-frame zero pressure. A temporary model copy with both NDI and NI
hardware subsystems removed compiled successfully and simulated for 0.5 s;
the logged synthetic pressure range was 0–0.8002 bar. No hardware controller
run has been performed. The two derivative time constants and 0.8 bar
deadzone are starting values, not identified plant parameters. NDI-derived
lengths are inferred from orientation and are not independent muscle-strain
measurements.

## 2026-10-01 — Renamed checkout and handoff to a new Codex task

**Repository location and Git state.** The experiment-computer checkout was
renamed to `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl`.
It remains the same Git repository on `codex/feedback-linearization-clean`,
with origin at `https://github.com/DulanjanaPerera/hybrid-soft-arm-dynamics-multisection-pointmass`.
The rename did not change commits or branch history. Commit `03781f0`
updated the absolute calibration-file paths in `NDI_readingOnly.slx`,
`Dynamic_model_Val.slx`, and `feedback_linearization.slx`; all three paths
were checked against the existing calibration MAT-file. Commit `c34351c`
is the user's feedback-model signal-routing and logging revision. Local
commits are ahead of origin; no push is implied by this handoff. The old
Codex task `01a0e929-9d85-74d1-ba85-cfd8002b4d07` retains the deleted
directory as its default working directory even though the saved local
project `module_dynamic` points to the new directory. Start a new task in
that existing project, using the current local checkout with Worktree
unchecked. The old task can be read for relevant history, but this log and
the repository are the durable technical record.

**Physical and measurement contract.** The arm is a single hanging module
with active length 174.11 mm, actuator radius 13 mm, a 55 mm base-sensor
offset, and a 50 mm tip-sensor offset. NDI sensor `0A` is base and `0B`
is tip. The NDI arm-frame mapping is `x_arm=x_NDI`,
`y_arm=z_NDI`, `z_arm=-y_NDI` (the 90-degree X rotation already in the
reader). The NDI port is COM12 at 921600 baud and the nominal BX polling
period is 0.05 s; a valid pair is from one BX poll. The current straight
calibration is `ndi_orientation_recordings/calibration_20260930_190505_898.mat`.
It was recorded from 100/100 accepted pairs with 2.30 mm straight XYZ
reference error; a separate preflight passed 98/100 paired samples with
0.0502 s median poll interval. Calibration is a sensor-mount reference,
not a correction for a bent pose. Verify the mounting and run preflight
again before a new hardware-control session. A lost or invalid pose must
not be treated as a measurement; the reader holds display values but
retains explicit validity flags. The three Simulink files above now refer
to the renamed absolute calibration path.

**Model and controller contract.** The reduced state is
`q=[l2;l3]` in metres with `l1=-l2-l3`. NDI `lengthChange_m` is
inferred from calibrated orientation and constant-curvature geometry;
it is not independent muscle-strain sensing. The nominal one-section
model uses downward gravity `[0;0;-9.81]` m/s^2, mass 0.1 kg,
`K=1350*[2 1;1 2]` N/m, `D=40*eye(2)` N s/m, and the existing
nonlinear length-bound penalty. `feedback_linearization.slx` has
desired and measured derivative blocks, direct desired-minus-measured
length/velocity/acceleration errors, a virtual acceleration
`v=ddl_des+Kd.*edl+Kp.*el`, and inverse dynamics. Initial editable
gains are `Kp=[4;4]` 1/s^2 and `Kd=[4;4]` 1/s. Desired lengths are
sampled at 0.05 s and filtered with 0.10 s velocity and 0.15 s
acceleration time constants. Measured lengths use valid NDI BX poll
timestamps and 0.20 s velocity and 0.30 s acceleration filters. Invalid
readings hold filter state and gate the inverse command to zero.

**Pressure contract.** The inverse maps generalized force through
`B=pi*(0.013/2)^2*1e5*[-1 1 0;-1 0 1]` N/bar, chooses the
minimum-common-pressure nonnegative P1/P2/P3 representation, adds a
provisional 0.8 bar command deadzone on active channels, and caps the
command at 3 bar by scaling effective pressure when needed. Positive
commands are valve commands, not measured chamber pressure. The
deadzone, stiffness, damping, and pressure response have not been
identified from closed-loop data; tubing and valve dynamics are not in
this baseline. The existing sine-wave reference sources remain
`amplitude=1.5`, `bias=1.5`, `phase=3*pi/2`, hence reach 0–3 rad
for both theta and phi. These are too large as an unreviewed first
hardware target. The user's later dashboard edits set the theta slider
range to about -3.14–3.14 rad and phi slider maximum to 3.14 rad;
the separate, currently commented phi Constant is 0.5. Check which
source is actually wired before a run. The NI valve subsystem can
energize physical hardware. A subsequent dashboard save had left the
valve Kill Constant at 1; it was reset and saved as 0 for this handoff.
This preserves the user's slider edits while leaving valve output gated.

**Validation performed and limits.** Synthetic MATLAB checks covered
missing-frame hold, sine-wave filtered derivatives, error sign,
nonnegative pressure, force reconstruction, saturation, and invalid-frame
zero pressure. A throwaway model copy with both NDI and NI subsystems
removed compiled and simulated for 0.5 s; its pressure signal ranged
from 0 to 0.8002 bar. That smoke test preceded the user's later
signal-routing/logging commit `c34351c` and dashboard edits; the
updated model was loaded to verify Kill=0, but has not had a second
offline compile or simulation. The actual hardware-bearing model has
not been run as a closed-loop controller or validated against hardware
timing.
The six 0–3 bar forward-model trials and estimator comparison remain in
`forwardModelValidation/` and `estimator_comparison/`; the earlier
entries above record their analysis and limitations. Do not infer
controller performance from the offline smoke test or from commanded
pressure alone.

**Instructions for the next task.** First verify the new working
directory, branch, HEAD, remote, and `git status`; preserve any
uncommitted Simulink edits. Read this log, `README.md`, the three
Simulink models, `ndiSupportFnc.m`, `AuroraReader_s2_Position_only.m`,
`AuroraDriver_2.m`, and `armS_single_entry.m` /
`armS_single_dynamics.m` / `armS_single_core.m`. Confirm current
block wiring, units, calibration paths, signal logging, valve Kill
state, and desired-reference source before changing code. For offline
checks, use a temporary model copy with the NDI and NI subsystems
removed; do not treat an ordinary simulation of the real model as
hardware-free. Before a user-directed physical test, establish the
straight reference/preflight, verify COM12 ownership and valid paired
samples, choose a modest theta/phi target, inspect the 0–3 bar command
path and kill gate, and log validity, poll timing, desired/measured
lengths and angles, pressure commands, saturation, and model-valid
diagnostics. Compare model predictions with those measurements before
changing physical parameters or introducing a learned residual model.

## 2026-10-02 — Length-bound expansion

At the user's request, increased the modeled reduced-coordinate length bounds
from +/-30 mm to +/-35 mm in Dynamic_model_Val.slx and feedback_linearization.slx,
the offline simulation runner, and the forward/inverse stiffness sweep.
Windows denied writes to Dynamic_model_Val_FWDnINV.slx and Inv_Dyn_val.slx;
those two saved models remain unchanged pending release of their file locks. The
nonlinear penalty magnitude and transition sharpness are unchanged. Existing
model wiring and user edits are preserved. This is a model assumption change,
not verification of a physical length limit. Historical inverse-analysis
settings remain at +/-30 mm. No hardware was operated.

There is no separate 15 mm cap on l1/l2. For the symmetric pressure case the
model yields l1=l2=-l3/2; a roughly -30 mm l3 plateau consequently produces
roughly +15 mm l1/l2. Plot units are mm: 0.015 m is 15 mm, not 0.015 mm.

**File-lock follow-up.** After the user reported closing both models and
requested completion, updated Dynamic_model_Val_FWDnINV.slx forward and
inverse bounds to +/-0.035 m. Verified only those two Stateflow XML members
changed and all XML parses. Windows still denied writing Inv_Dyn_val.slx;
that model remains pending at +/-0.030 m. No hardware or simulation was run.

## 2026-10-05 - Follow-up experiments, plotting, and hysteresis/MPC direction

**Standing collaboration preference.** The user prefers discussion by default and
does not authorize edits to code, Simulink models, or other files unless explicitly
requested. On 2026-10-05 the user separately authorized ongoing updates to
PROJECT_LOG.md whenever new project code is written or an important technical
discussion occurs. This standing authorization applies to this log only.

**Prior identification as a possible PAC baseline (discussion, 2026-10-04).**
The earlier continuum-arm system-identification handoff report is at
C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\Controlling\Dynamic_jointspace_feedbacklinearization\ChatGPT\continuum_arm_hysteresis_sysid_handoff_report.md.
Its P2/P3-chirp simulation fit gave a coupled stiffness matrix approximately
[2576.3, 1271.0; 1271.0, 1886.0] N/m and damping diag([16.15, 18.50])
N s/m. Its fitted Bouc-Wen *mechanical-force* state used
hdot = alpha_h*dq - beta_h*abs(dq).*h - gamma_h*dq.*abs(h), with
alpha_h=65.59, beta_h=49.04, gamma_h=-19.93. This is separate from
pressure-command deadzone and pressure-side hysteresis. These values are
candidate initial estimates for a future passivity-based adaptive controller
(PAC), not validated drop-in parameters. The older simulation passed
[0,l2,l3] to its dynamics and used direct [F2,F3] input on P2/P3-only
trials; the current model contract uses l1=-l2-l3 and [F2-F1,F3-F1].
Those force mappings agree when P1=0 but differ when P1 participates.
A coupled adaptive K must remain positive definite, and K/D adaptation
does not by itself learn valve memory. The user's more recent forward/inverse
comparison suggested scalar stiffness near 650 N/m for its trajectories;
that estimate and the older full matrix need comparison on common held-out
data, geometry, and pressure conventions.

**Closed-loop test testing_20261004_1311.mat (analysis, 2026-10-04).**
The user ran the hardware-bearing feedback_linearization_trackingWorks.slx
and saved feedbacklinearization_tests/testing_20261004_1311.mat. Earlier
statements above that no closed-loop hardware run had occurred describe the
October 1 handoff only. Inspection of the saved trackingWorks model, whose
timestamp predates this recording, found inverse-block stiffness
K=250*[2 1;1 2] N/m, despite the later 650 N/m forward-model estimate.
Its separate forward block has a different K; always check the block used
by the recorded run. No model was edited in this analysis.

Over the 95-125 s hold, desired/measured phi averaged 1.768/1.049 rad.
Reduced length errors averaged [+9.37,-4.01] mm. Logged
acceleration-correction force averaged [10.04,0.79] N, combined gravity
and stiffness [5.97,-0.41] N, Coriolis-plus-damping near zero, and requested
generalized force [16.01,0.39] N. Mean pressure command [P1,P2,P3] was
[0,2.006,0.829] bar. Inverse validity remained true, pressure limiting
remained false, and reconstructed available force matched requested force
only by the controller's assumed pressure/deadzone mapping; actual chamber
pressure was not measured. P3's modeled effective pressure was about
0.029 bar. At the measured pose, the first generalized gravity component
was about 1.03 N and K=250 stiffness contributed about 4.94 N; K=650 would
contribute about 12.84 N. Gravity alone does not explain the plateau.
With PD-only correction, a force-model or valve mismatch can leave steady
tracking error without the pressure command continuing to rise.

**Continuous circle and plotting code (2026-10-05).** A continuous circular
reference uses unwrapped theta_des(t)=2*pi*t/T and a nonzero desired phi;
the configuration-to-length sine/cosine terms make desired lengths repeat
without numerically resetting theta. Saved tests are
feedbacklinearization_tests/circle_1_round.mat (62.65 s),
circle_2_round.mat (124.20 s), and circle_multi_rounds_20s.mat (64.05 s).
Desired theta is logged unwrapped; measured orientation theta wraps at
+/-pi, and desired phi rises as high as 1.5 rad.

At the user's request, feedbacklinearization_tests/plotCircleTracking.m was
added. Passing out directly or a circle MAT path opens seven figure groups:
desired/measured theta and phi, all three actuator length changes, three
pressure commands, 3D desired/measured tip-sensor XYZ, angle errors, length
errors, and XYZ component/norm errors. It returns data and figure handles
without saving plots or running Simulink/hardware. Desired XYZ uses the
constant-curvature sensor-offset geometry [0.055,0.17411,0.05] m; measured
XYZ comes from NDI position. Measured theta/phi and lengths use orientation,
a distinct measurement from XYZ, so their reconstructed task positions need
not coincide exactly with measured XYZ. The script masks invalid samples,
unwraps measured theta for trajectory comparison, and uses the shortest
signed angular difference for theta error. It was run offline against all
three circle MAT-files and with an out object passed directly; each
produced seven figures.

**Deadzone and forward-model MPC idea (discussion, 2026-10-05).**
The inverse pressure allocator sets inactive channels to zero, but when a
channel requests any positive effective pressure it commands
p_effective+0.8 bar. This makes a mathematical command jump from zero to
about 0.8 bar at activation, independent of any mechanical hysteresis.
Valid channel switch-ons in circle_multi_rounds_20s.mat started around
0.80-0.83 bar. At least two channels were commanded in 98.8% of its
0.05 s samples. Actual chamber pressure and true valve switching
thresholds were not logged. Adding the old mechanical Bouc-Wen force term
to the left side of the EoM would not remove this allocator discontinuity.
A separately validated pressure-side hysteresis/valve state would belong
on the actuator-input side.

The old project's tapped-delay forward MLP uses pressure, pressure-rate,
loading-state and previous-length histories to predict active-actuator
length. Its direct length-to-pressure inverse performed poorly. A proposed
nonlinear MPC would instead optimize future pressure sequences using a
validated forward predictor, penalizing future length error, pressure
changes and valve switching under pressure and rate limits; apply only
the first command, then replan from new measurements. This could resolve
the multivalued inverse caused by loading/unloading history and pressure
redundancy. The old forward models were trained as active-only outputs,
whereas current circles usually command multiple actuators together, so
their cross-actuator and multi-step predictive accuracy is unproven.
Before implementation, evaluate one-step and recursive multi-step
predictions on held-out current circle data, feeding *predicted* future
lengths into the tapped history during rollout. If performance fails,
identify a full three-pressure/two-length forward model using
simultaneous-actuation trials. Keep pressure/valve dynamics distinct from
passive mechanical K, D and hysteresis to avoid double counting.
No controller or Simulink files were changed during this discussion.

## 2026-10-05 - Backfill of inverse and planar stiffness validation

The offline inverse-dynamics analysis lives in
inverDynamic_test/analyzeInverseDynamics.m and reads
inverDynamic_test/simulationresulta_20261001_1936.mat. Its saved report
found 2,659 samples over 167.836 s of physical BX poll time, 99.85%
valid measurements, 7.48% pressure-limited samples, three qualifying
reference holds, two stationary tails, and one target reached and
maintained under its stated 1 mm/two-second rule. The analysis reports
tracking, force decomposition, validity, and timing but cannot identify
dynamic terms from zero-rate holds or establish actual chamber force
from software-reconstructed available force. Outputs and assumptions are
documented in inverDynamic_test/README.md and results/.

The forward/inverse comparison was replayed offline in
FWD_INV_COmparision/sweepStiffness.m using
FWD_INV_COmparision/Comparision_results.mat. Only the scalar in
K=k*[2 1;1 2] was swept from 1350 down to 200 N/m while retaining
D=40*eye(2), the 0.8 bar assumed deadzone, 6.5 mm effective-pressure
radius, 13 mm actuator offset, gravity, and +/-35 mm bound penalty.
The saved sweep report selected 650 N/m from this grid, with 3.55231 mm
RMS error in the independent reduced coordinates over 600 valid
samples. The 1350 baseline replay agreed with the logged forward
trajectory to within 4.61e-7 m. This is a same-recording conditional fit;
it cannot independently identify material stiffness, pressure calibration,
hysteresis, timing, or actuator symmetry. The score table and plots are
under FWD_INV_COmparision/stiffness_sweep_results/.

The P1-only stepped-ramp experiment and analysis are under
planar_stiffness_test/. Each pressure level is held for ten simulation
seconds over three loading/unloading cycles. The offline
estimatePlanarStiffness.m fits a static effective scalar stiffness per
hold from commanded P1 and orientation-derived lengths. It deliberately
assumes zero deadzone, and its post-Kill command is not measured chamber
pressure. Loading and unloading are pooled by pressure level.
fitStiffnessCurve.m builds an empirical PCHIP from those means;
stiffnessFromPhiBlock.m is a MATLAB Function block version with endpoint
extension. The saved mean estimates range from 4132 N/m near phi=0.094
rad to about 953 N/m near phi=2.056 rad, with very large low-angle
standard deviations (1843 and 2597 N/m at the first two levels).
The user reported poor agreement, particularly on unloading, when using
the phi-only stiffness function on stepped-triangle and sine tests.
These values therefore reflect pressure/deadzone/hysteresis and model
assumptions as well as any mechanical stiffness. The apparent low-angle
jump must not be treated as proven material stiffening. The fit is
plane-specific and not validated as a controller or general forward law.
The test files, curves, mean tables, and caveats are in
planar_stiffness_test/README.md and results/.

## 2026-10-05 - Common-pressure allocation and possible co-contraction

The user proposed a roughly 2 bar command on all three PMAs at neutral,
then bending by lowering the opposing channel while increasing the other
two without crossing the valve deadzone. This is a discussion and
experiment proposal; no controller or Simulink code was changed.

The current ideal input map is tau = f*[-1 1 0;-1 0 1]*p_effective,
where f=pi*(0.0065)^2*1e5 = approximately 13.273 N/bar. Its nullspace
is span([1;1;1]): equal *effective* pressures cancel from the two
modeled bending forces. When all command pressures exceed the assumed
0.8 bar deadzone, an equal common command also cancels algebraically.
At [2,2,2] bar the nominal effective pressures are [1.2,1.2,1.2]
bar and modeled bending force is zero. For desired theta=0, the
configuration-to-length convention requires P2 and P3 above P1.
For example, [1,3,3] bar gives differential pressure [2,2] bar and
nominal generalized forces about [26.55,26.55] N while all three
channels remain above the nominal deadzone.

Starting at 2 bar does not inherently reduce the maximum nominal
differential pressure if the common level may change: approaching
[0.8,3,3] reaches the same 2.2 bar difference as the existing
minimum-common-pressure allocator with P1 inactive. A guaranteed
above-deadzone floor, a required margin around the real switching
threshold, a fixed mean pressure, and finite pressure rates do reduce
the reachable differential range. The allocator should choose common
pressure within per-channel lower/upper bounds and report infeasible
requested forces rather than silently assume 2 bar is always possible.

This nullspace exists in the current reduced *force map*, not
necessarily in the physical arm. Co-pressurizing PMAs may change
incremental bending stiffness, axial shortening, stored energy, valve
dynamics and hysteresis; the current EoM holds K constant and omits
a common-length coordinate, so it predicts none of those effects.
Pressure-dependent stiffness is plausible for antagonistic PMAs but
has not been identified for this module. NDI orientation-derived
lengths mainly reflect bending; log tip XYZ/Z as well to detect
common-mode axial motion. A useful validation is to vary common
pressure at fixed pressure differences, then apply small differential
perturbations at several common levels while logging loading and
unloading, measured pose, and actual pressure if available.

A read-only inspection of feedback_linearization_trackingWorks.slx
on 2026-10-05 found its inverse block had since been changed to
K=650*[2 1;1 2] N/m. It still uses the hard zero-versus-
p_effective+0.8 bar channel activation rule. The earlier K=250
finding applies to the timestamped 2026-10-04 13:11 recording, not
automatically to the current saved model.

## 2026-10-05 - Backbone constraint and pressure nullspace clarification

The user clarified that this module has a backbone restricting axial length. In the ideal fixed-backbone, symmetric three-PMA geometry, physical actuator length changes satisfy dl1+dl2+dl3=0; the third length is dependent on the other two. This is the reason for the two-coordinate bending model. Equal actuator forces/pressures have zero generalized work in those two bending coordinates: with dl1=-dl2-dl3, F1*dl1+F2*dl2+F3*dl3=(F2-F1)*dl2+(F3-F1)*dl3. Thus the common-pressure direction is a force-allocation nullspace under the ideal pressure-to-force model, not an omitted free axial degree of freedom.

Correction to the preceding discussion: for a sufficiently stiff backbone, common pressurization should not be described as causing appreciable common-mode axial shortening. It may instead change backbone load/prestress, PMA bulging and effective incremental bending stiffness, so the physical pressure-to-bending map may still depend on the common level. Tip Z may change from bending at fixed backbone arc length; Z change by itself is not evidence of axial shortening. Any proposed common-pressure baseline should therefore be evaluated by holding pressure differences fixed and measuring orientation/length response, then comparing small differential bends at several common levels. Actual chamber pressure sensing would separate valve behavior from mechanical effects. This is a discussion only; no controller or Simulink files were changed.

## 2026-10-05 - Two-computer Git handoff discussion

The project is used from two computers. Local branch codex/feedback-linearization-clean was at 028c3ef and tracked origin/codex/feedback-linearization-clean, but was 10 commits ahead of the last fetched remote-tracking ref (0664270). New local branch codex/common-pressure-testing was created from the same 028c3ef and had no upstream. Modified and untracked experiment files remained in the shared checkout and were not committed on either branch. Discussion: pushing branches is appropriate for the two-computer workflow, but only committed files transfer; review and commit intended changes before relying on the second computer, fetch remote changes before pushing, and avoid force-pushing shared branches. No commit, fetch, or push was performed during this discussion.
