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

## 2026-10-05 - Common-pressure test design (discussion)

Read-only inspection of the current feedback_linearization_trackingWorks.slx inverse block confirms forcePerBar=pi*(0.013/2)^2*1e5=13.2732 N/bar and tau=forcePerBar*[p2-p1;p3-p1] for effective pressures. The allocator chooses the minimum nonnegative effective pressure triplet and then commands each active chamber at effective+0.8 bar. Therefore at least one chamber is nominally inactive for every ordinary request; existing controller runs cannot directly identify a three-chamber common-pressure bias effect.

For desired generalized force tau=[tau_l2;tau_l3], define d2=tau_l2/forcePerBar and d3=tau_l3/forcePerBar. If all commands stay above the individual deadzones, every ideal solution is [p1;p1+d2;p1+d3]. A fixed command mean pbar=2 gives p1=2-(d2+d3)/3, p2=p1+d2, p3=p1+d3. With a chosen command floor pmin (strictly above measured valve activation) and 3 bar ceiling, p1 is feasible only in [max(pmin,pmin-d2,pmin-d3), min(3,3-d2,3-d3)]. The interval must be checked before allocation. A flexible common level can retain more bending range than a hard fixed 2-bar mean. The logged 2026-10-04 hold force [16.01,0.39] N illustrates the arithmetic: a 2-bar mean would command approximately [1.588,2.794,1.618] bar, subject to the assumed equal deadzone; this is not a proposed replay of that recorded hardware run.

Recommended first experiment: with the controller disabled, ramp equal P1/P2/P3 commands through a few validated common levels, hold each and record orientation-derived lengths, theta/phi, NDI XYZ, command after Kill/saturation, validity and timestamps; actual chamber pressures would be especially valuable if instrumentation becomes available. At each common level apply small positive and negative differential pressure perturbations while preserving the same mean, then compare equilibrium bend and incremental force-versus-length response. Repeat loading and unloading. The fixed backbone constrains the sum of physical length changes; any tip Z change alone does not prove axial shortening. No Simulink/controller code was edited and no hardware was run during this design discussion.
## 2026-10-05 - Common-pressure test folder and staged experiment plan

At the user's request, moved the newly created common_pressure_testing.slx duplicate from the repository root to common_pressure_testing/common_pressure_testing.slx, byte-for-byte (SHA-256 67FB1A2FDE01FBD30337AF8A5DDE4188007EFFA6FD7C4149CC84119D942AE60B). The original Dynamic_model_Val_FWDnINV.slx was already modified independently in the checkout and was left untouched. No Simulink blocks were edited or model run.

Read-only inspection of the duplicate found P1/P2/P3 root Constant blocks saved at zero; their Mux drives both the forward model and the pneumatic-valve path. The inverse block's pressure output is used for comparison/display and is not wired as the valve command. The forward model uses the assumed 0.8 bar deadzone, 0-3 bar command range, K=650*[2 1;1 2] N/m, and tau=forcePerBar*[-1 1 0;-1 0 1]*p_effective. The model contains a valve Kill dashboard control; inspect its live/saved state before hardware use. No chamber-pressure sensor is available, so command-side tests cannot isolate regulator pressure dynamics from mechanical effects.

Proposed sequence, pending explicit request to modify blocks/scripts or run hardware: (1) begin straight, unpressurized and pass the normal NDI preflight; verify Kill and zero commands; (2) with open-loop pressure commands, slowly apply equal P1=P2=P3 plateaus from a low validated level toward the proposed 2 bar baseline, returning to zero, and log command and post-Kill pressure, NDI validity/timestamps, orientation theta/phi, orientation-derived three lengths, and XYZ; quantify neutral-pose drift, settling and repeatability; (3) at several common levels, apply small positive/negative fixed-mean theta=0 planar perturbations p=[c-2*delta/3;c+delta/3;c+delta/3], with every channel within a measured above-deadzone margin and 3 bar; repeat loading/unloading; (4) compare equilibrium bend and local differential compliance across common levels, separating model-predicted response from measured response. Only after these tests should a common-pressure allocator or pressure-dependent stiffness term be selected. Save future recordings and analysis scripts in common_pressure_testing/ rather than the repository root.
## 2026-10-05 - Single-block common-pressure protocol code

At the user's explicit request, created common_pressure_testing/commonPressureSequenceBlock.m as a copy-paste MATLAB Function block implementation; no Simulink model or NI/NDI hardware was run or modified. It takes simulation Clock time t (seconds), returns scalar P1_bar/P2_bar/P3_bar for the model's existing three-input Mux, and phaseID for logging. The fixed 300 s schedule is 0 bar for 10 s, then successive 90 s groups at common command levels 1.2, 1.6, and 2.0 bar, then a 5 s return ramp and 15 s zero hold. Each common-level group ramps to the level and holds equal pressure, then tests +0.2 bar and -0.2 bar differential pressure with intervening equal-pressure holds. Cubic smoothstep ramps last 5 s; all active-level commands stay between 1.067 and 2.133 bar, above the nominal 0.8 bar deadzone. Initial/final activation ramps necessarily cross the deadzone. Positive differential uses P1=c-2*delta/3 and P2=P3=c+delta/3, maintaining constant command mean and P2-P1=P3-P1=delta.

MATLAB batch evaluation of the pure function at 6001 samples (0:0.05:300 s) passed bounds, above-deadzone active interval, representative positive/negative phases, invalid/negative-time zero output, and final zero output. This was not a Simulink compilation or hardware validation. The user should connect the three scalar pressure outputs to the current P1/P2/P3 Mux inputs, connect t to a Clock, log phaseID and pressure/NDI signals, verify the valve Kill state, and set Stop Time to 300 s when ready for a physical experiment.
## 2026-10-05 - Common-pressure 300 s data analysis

At the user's explicit request, added common_pressure_testing/analyzeCommonPressureTest.m and README.md. The offline function reads test_300s.mat/out.logsout, selects duplicate-named model signals by block path, aligns all channels to the 0.05 s NDI time grid, checks orientation/sensor/configuration validity, and summarizes each phase from its final five simulation seconds. It computes signed U=(l2+l3)/0.013 and cross-plane V=(l3-l2)/(sqrt(3)*0.013), measured and forward-model lengths/angles, commanded pressure, regulator-path command, NDI XYZ, physical poll timing, phase-tail slope/range, and positive/negative differential response at each common level. Outputs are CSV tables, a text report, five PNG figures, and an aligned MAT under common_pressure_testing/results/test_300s/. No Simulink model or hardware was opened or run for this analysis. MATLAB batch execution and export checks passed for 17 phases, three common levels, and 6001/6001 valid NDI samples.

In this recording, equal command pressure produced nonzero NDI bending although the symmetric forward model predicts zero: mean neutral U was 0.1065, 0.1702, 0.1958 rad and cross-plane V was 0.0717, 0.1444, 0.1523 rad at common commands 1.2, 1.6, 2.0 bar (averaging the three neutral phase tails at each level). The measured signed-U gain between +0.2 and -0.2 bar differential phases was 0.3096, 0.4649, 0.3897 rad/bar, versus 0.9641 rad/bar from the saved forward model. Returning to equal pressure after the positive step shifted U by +0.0490, +0.0574, +0.0374 rad; after the negative step it shifted by -0.0597, -0.0743, -0.0493 rad. The largest within-tail five-second U range was 0.00477 rad, so a tail average is not proof of equilibrium. Logged post-Kill/regulator-path command matched requested command to the printed precision (max difference 0.00000 bar), but there was no chamber-pressure measurement. Median physical NDI poll interval was 0.0508 s and physical poll span 304.58 s over 300 simulation seconds.

The conditional apparent scalar K values in the report (about 2199, 1464, 1747 N/m) assume the nominal equal pressure-force gain, symmetric K=k*[2 1;1 2], and equilibrium; they are not measured material stiffness. Unequal actuator/valve gains, preload, pressure-side dynamics, cross-plane response, and hysteresis/drift could explain the observed mismatch. The existing controller EoM cannot reproduce equal-command bending with its ideal pressure nullspace and no bias term. The user-provided common_pressure_testing.slx and independently modified Dynamic_model_Val_FWDnINV.slx were not edited.
## 2026-10-05 - Common pressure and PD feedback discussion

The 300 s open-loop recording supports testing common pressure as a way to keep all three valves above the assumed 0.8 bar activation threshold during small differential commands. It does not establish a deadzone-free physical pressure response: the initial/final ramps cross the threshold, actual chamber pressures were not measured, and valve hysteresis/lag or unequal channel gains may remain. Equal-pressure commands generated substantial nonzero bend and cross-plane response, increasing with common level; the symmetric model treats equal pressure as an exact force nullspace and misses this behavior. The measured differential U gain was also less than half the modeled gain at all three common levels.

With a constant unmodeled generalized force, ideal PD feedback linearization has an approximate steady-state relation Kp*e_ss=disturbance (up to error-sign convention); the derivative term vanishes at rest. PD therefore can reduce a bias-dependent tracking error if the pressure allocation has sufficient authority and gains are stable, but cannot guarantee zero steady-state error. Higher gains alone may aggravate pressure-rate limits, delays, oscillation, or saturation. Before closed-loop use, characterize a neutral trim pressure triplet and local differential response at a practical common level, account for actuator bounds in the allocator, and then test conservative PD tracking. A feasible integral, disturbance-estimate, or learned-equilibrium trim term may be needed for zero constant-target error. This entry records discussion only; no controller, Simulink, or hardware changes were made.
## 2026-10-05 - Proposed implementation sequence for common-pressure control

Discussion-only plan following the 300 s common-pressure recording. First choose a feasible command-pressure floor from measured valve behavior, with margin above the switching region, and compare candidate common levels. The current 2.0 bar equal-command hold bent the arm more than 1.2 bar, so 2.0 bar should not be assumed to be the best operating point. At each candidate level, identify an unequal P1/P2/P3 neutral trim that holds the chosen reference orientation, if such a trim is feasible within the command floor and 3 bar ceiling. Repeat approach from increasing and decreasing pressure and allow for settling; commands are logged but chamber pressures are not measured.

Next, characterize small positive and negative differential perturbations around the neutral trim in both bending directions. Estimate the local command-pressure-to-bending/length response and check whether it changes with common level, direction, and pose. Use this calibration to form a feedforward trim and a pressure allocator. An initial ideal allocator can map desired generalized force tau to pressure differences d2=tau_l2/forcePerBar and d3=tau_l3/forcePerBar, then choose the common offset so all three commands remain within the verified floor and 3 bar; the calibrated neutral triplet may be used as the offset reference. The allocator must detect and report infeasible demands rather than silently clip one valve or change the intended pressure differences. The current minimum-common-pressure allocator deliberately leaves at least one channel inactive and would need modification for this approach.

Only after open-loop calibration should feedback-linearization tracking be tested with conservative PD gains and the new allocator. Log desired/measured bend and lengths, controller force components, raw allocated pressure, post-Kill pressure command, saturation/feasibility flags, NDI validity and timestamps, and chamber pressure if instrumentation later becomes available. Compare identical reference trajectories using the old allocator and the common-pressure allocator, focusing on deadzone crossings, steady-state error, overshoot, settling, and pressure headroom. Pure PD cannot guarantee zero steady-state error from a persistent unmodeled bias; if a stable residual remains and authority is available, add a slow bounded integral/disturbance-estimate or calibrated equilibrium trim. Do not modify the dynamic stiffness term merely to force a fit to these command-side tests. No code, model, or hardware was changed for this discussion.
## 2026-10-05 - Reusing the common-pressure Simulink model for neutral-trim calibration

Discussion only. The existing common_pressure_testing/common_pressure_testing.slx can be reused for step 1 by replacing the prior time-sequence pressure generator with a new MATLAB Function block driven by Clock time. Its P1/P2/P3 outputs can feed the existing root pressure Mux, which drives the forward model and NI valve command path; a phase ID should be logged alongside the existing orientation, lengths, XYZ, and command signals. Keep the original 300 s generator as a separate file for reproducibility rather than overwriting it.

The previous 300 s recording already probes equal common levels and a symmetric P2=P3 differential, but has no independent P2-versus-P3 perturbation. To identify a two-direction neutral trim, propose short held + and - perturbations in both symmetric and antisymmetric pressure directions at candidate common levels, with ramps and returns to the baseline. A fixed-mean parametrization is p=[c-2*s/3, c+s/3-a, c+s/3+a], where s changes P2/P3 together relative to P1 and a separates P2 from P3. Verify every generated command stays above a selected active-channel floor and below 3 bar. At each hold, use valid NDI orientation-derived independent lengths or signed U,V rather than theta near a straight pose, and estimate a local two-input response. Solve offline for an unequal trim candidate that minimizes both bending components, then run a separate hold to verify it and check loading-direction dependence. The sequence is open loop; it should not automatically search toward an unverified pressure setting during the hardware run. No Simulink, controller, or script files were changed in this discussion.
## 2026-10-05 - Neutral-trim calibration block code

At the user's explicit request, added common_pressure_testing/commonPressureTrimSequenceBlock.m for copy-paste into a new MATLAB Function block in the existing common_pressure_testing.slx. It accepts Clock time t and outputs P1_bar/P2_bar/P3_bar plus phaseID. The deterministic 387 s open-loop schedule has 10 s initial zero, then three 119 s groups at common command levels 1.2, 1.6, and 2.0 bar, followed by a 5 s ramp to zero and 15 s final zero hold. Each group contains a 5 s baseline ramp, 10 s baseline hold, and four 26 s pulses: symmetric +0.2/-0.2 bar and antisymmetric +0.1/-0.1 bar. Each pulse ramps out for 4 s, holds for 10 s, ramps back for 4 s, and holds the baseline for 8 s. Pressure commands use p=[c-2*s/3, c+s/3-a, c+s/3+a], preserving the mean c. All active-level commands are 1.067-2.133 bar and under 3 bar; initial/final ramps cross the nominal deadzone. Phase IDs identify each ramp/hold; 900 is final ramp and 901 final zero hold.

A 0.05 s-grid independent arithmetic sweep verified 7741 samples, active minimum 1.066667 bar, maximum 2.133333 bar, expected antisymmetric command at t=82 s, 57 unique phase IDs, and zero output after 372 s. MATLAB batch startup failed even for disp(1) with 'System Error: File system inconsistency', so MATLAB execution and Simulink compilation were not verified. No Simulink file or NI/NDI hardware was opened or run.
## 2026-10-05 - 387 s common-pressure neutral-trim analysis

Analyzed the user-recorded common_pressure_testing/trim_test_387s.mat offline with new common_pressure_testing/analyzeCommonPressureTrimTest.m. The script reads saved out.logsout only, aligns signals on the 0.05 s NDI grid, uses orientation/sensor/configuration validity, summarizes the final 5 s of 29 held phases, and writes phase_summary.csv, trim_estimates.csv, report.txt, two PNG plots, and analysis_results.mat under common_pressure_testing/results/trim_test_387s/. MATLAB batch execution passed after sandbox escalation; no Simulink model or NI/NDI hardware was run. The new user-created common_pressure_trim_testing.slx was not edited.

All 7741/7741 NDI samples met the logged validity checks. Initial/final zero-pressure signed bending was U=0.03016/0.03085 rad and V=0.01488/0.01515 rad, so the analysis uses their mean (U=0.03051, V=0.01501 rad) as the reference pose instead of assuming exactly zero angle. Logged command and regulator-path signals agreed to printed precision, but no chamber pressure was measured. The raw logged Kill dashboard state was 1; this alone does not establish physical valve state. Median physical NDI poll interval was 0.0508 s and physical poll span 393.23 s over 387 simulation seconds.

At common commands 1.2, 1.6, 2.0 bar, mean equal-command U,V were respectively (0.1001,0.0738), (0.1655,0.1442), (0.1908,0.1479) rad. A local two-input command-to-bend Jacobian was estimated from +/- symmetric s=0.2 bar and antisymmetric a=0.1 bar holds relative to adjacent equal-pressure holds. The [dU/ds,dV/ds] columns were (0.2802,-0.0422), (0.4305,-0.0777), (0.3683,-0.0667) rad/bar; [dU/da,dV/da] were (0.1024,0.4298), (0.1125,0.5973), (0.1035,0.5392) rad/bar. Matrix condition numbers 1.57,1.39,1.47 indicate locally independent responses in these two command directions; the Jacobian is not a pressure-to-force map.

Linearized pressure-triplet estimates to return to the zero-pressure reference were [1.328,1.292,0.981] bar at mean 1.2; [1.766,1.766,1.268] at 1.6; and [2.236,2.172,1.592] at 2.0. All require extrapolation beyond at least one tested +/- perturbation, so none is validated for direct use as a fixed operating trim. The 1.2 candidate also puts P3 below the provisional 1.0 bar active-channel floor, although still above the nominal 0.8 bar deadzone. The 1.6 candidate stays within provisional 1-3 bar bounds and has more headroom than 2.0, making it a reasonable next validation target, not a verified trim.

Equal-pressure U varied across the five baseline holds by 0.0585, 0.0705, 0.0502 rad at the three levels; corresponding V spans were 0.0498, 0.0430, 0.0319 rad. These history-dependent returns are much larger than the initial-to-final zero-pressure drift and cannot be ignored when selecting a fixed trim. Largest within-hold 5 s U/V ranges were 0.00861/0.00649 rad. Next validation should approach a candidate in bounded increments at one chosen common level, hold and measure both signed bending coordinates, then repeat from a different loading history; do not infer true pressure or stiffness from this recording alone.
## 2026-10-05 - Next test after neutral-trim analysis

Discussion-only recommendation: verify the 1.6 bar mean-pressure candidate before changing the feedback-linearization allocator. The current local linear estimate is P1=1.766, P2=1.766, P3=1.268 bar, targeting the recorded zero-pressure reference U=0.03051,V=0.01501 rad. This candidate is inside provisional 1-3 bar command bounds but extrapolates beyond the tested antisymmetric +/-0.1 bar perturbations and is not yet calibrated. Use the existing open-loop test model with a controlled sequence: begin at zero pressure to record the reference, ramp and hold [1.6,1.6,1.6], approach the candidate in bounded intermediate steps, hold long enough to judge settling, return to baseline or zero, then approach the same candidate from a different pressure history. Log phase ID, command and post-Kill regulator-path signals, signed U,V, orientation-derived lengths, phi/theta, XYZ, validity, and physical poll timestamps. Compare the two candidate holds with the zero-pressure reference and with each other; large history dependence would argue against treating one fixed triplet as a complete equilibrium calibration. Use the measured response to update the local trim estimate if necessary. Only after repeatable verification should a common-pressure allocator be integrated with conservative feedback-linearization PD tracking. No code, Simulink model, or hardware was changed for this discussion.
## 2026-10-05 - Open-loop validation block for the 1.6 bar trim estimate

At the user's explicit request, added common_pressure_testing/commonPressureTrimValidationBlock.m as paste-ready MATLAB Function block code for a new open-loop Simulink test; no Simulink model was edited and no NI/NDI hardware was run. The block takes Clock time t and outputs P1/P2/P3 command pressures and phaseID. Stop Time is 199 s. Route A holds zero, equal 1.2 bar, equal 1.6 bar, intermediate [1.65,1.65,1.50] and [1.70,1.70,1.40] bar, then estimated trim [1.766,1.766,1.268] bar for 20 s. After returning through equal 1.6 and zero, Route B reaches equal 1.6 directly from zero and then ramps directly to the same estimated trim for a second 20 s hold. The run ends with a ramp and zero hold. Smoothstep ramps and unique ramp/hold phase IDs allow offline comparison; candidate hold IDs are 62 and 102, zero-reference holds are 12, 82, and 112. All target triplets preserve a 1.6 bar mean in the trim approach except the separate 1.2/zero phases. The candidate remains an unverified extrapolated estimate; the block is time-driven and does not close the loop on NDI.

MATLAB batch evaluation of the pure function on 3981 points from 0 to 199 s at 0.05 s spacing passed output bounds (0 to 1.766 bar), both 20 s candidate holds, expected phase IDs, and final zero output from 184 s. This was not a Simulink compilation or physical pressure validation.
## 2026-10-05 - 199 s neutral-trim validation result

Analyzed the user-recorded common_pressure_testing/trim_test_199s.mat offline with new common_pressure_testing/analyzeCommonPressureTrimValidation.m. It aligns saved logs to the NDI grid, summarizes the final 5 s of 11 hold phases, compares two candidate holds with nearby zero-pressure references, and writes hold_summary.csv, report.txt, two PNG plots, and analysis_results.mat under common_pressure_testing/results/trim_test_199s/. MATLAB batch execution passed; no Simulink model or NI/NDI hardware was run or edited. 3980/3981 samples passed logged validity checks; median physical NDI poll interval was 0.0508 s. Command and regulator-path software signals agreed to printed precision, but chamber pressure was not measured. Raw logged Kill dashboard state was 1 and is not independent proof of valve state.

The zero-pressure reference varied through the run: initial phase 12 U,V=(0.03209,0.02219) rad, middle phase 82 (0.02662,0.01480), and final phase 112 (0.02043,0.00879). Route A target was the mean of initial/middle zero holds, (0.02935,0.01849); Route B target was the mean of middle/final zero holds, (0.02353,0.01179). At the estimated candidate [1.766,1.766,1.268] bar, phase 62 reached U,V=(-0.07201,-0.13355) rad and phase 102 reached (-0.07344,-0.13616). Their absolute difference was (-0.00143,-0.00261) rad and length differences were [+0.019,+0.020,-0.039] mm, so the two routes landed near each other; however phase 102 still had final-five-second slopes about -0.00132 and -0.00143 rad/s, and the zero reference drifted more than their difference. Both candidate holds clearly overshot the desired pose: residuals were (-0.10136,-0.15204) and (-0.09697,-0.14795) rad relative to local zero references. Tip XYZ residual norms were 26.68 and 25.72 mm. Nonnegative phi alone would obscure that the signed U,V bend reversed direction.

The intermediate [1.700,1.700,1.400] bar hold (phase 52) was much closer: U,V=(0.04730,0.00706) rad, residual (+0.01795,-0.01143) rad relative to Route A zero reference, combined signed-bend residual 0.02128 rad, and tip XYZ residual norm 6.57 mm. Its V final-five-second slope was -0.00117 rad/s, so this point is not a validated equilibrium. The earlier small-perturbation Jacobian extrapolated poorly to the predicted candidate; this demonstrates a nonlinear/history-dependent command-to-pose response but does not isolate stiffness, valve behavior, or actual pressure. The next calibration should use smaller independent symmetric and antisymmetric perturbations centered around the observed near-target [1.7,1.7,1.4] bar point, with longer holds and repeated approaches. A one-dimensional adjustment along P1=P2 alone cannot be assumed to satisfy both U and V targets.
## 2026-10-05 - Manual versus automated neutral-trim search

Discussion only. The user proposed manually changing pressure commands and phase IDs during a variable-length test. Recommended a short manual local search around the observed near-target [1.700,1.700,1.400] bar point, followed by an automated replay of the selected setting from multiple pressure histories. This separates efficient trim discovery from repeatability validation. Use a fixed common command mean c=1.6 bar and two tunable differential inputs s,a with p=[c-2*s/3,c+s/3-a,c+s/3+a]; the observed near-target point has s=a=-0.15 bar. Vary one differential input at a time by small increments (roughly 0.02-0.03 bar initially), hold long enough to assess the ongoing settling seen in the 199 s data (about 20-30 s or more if U/V still drift), and compare signed U,V with a nearby zero-pressure reference. Keep commands within the chosen active-channel floor and 3 bar ceiling and apply smooth ramps or rate limits rather than sudden manual jumps.

For variable-duration testing, log the actual three command values, post-Kill regulator-path commands, U,V/lengths, XYZ, validity, and timestamps. A manually changed trial ID is useful as an annotation but should be a separate tunable input, not a source-code edit during a hardware run; offline analysis can segment stable pressure plateaus from the logged commands even if no fixed phase schedule exists. Use explicit trial IDs only after each target is set, and evaluate settled windows rather than a fixed offset from the ID change. Revisit zero pressure during the run because the measured reference drifted over the prior 199 s experiment. Once a promising triplet is found, repeat it from at least two loading histories before treating it as a controller equilibrium. No code, Simulink file, or hardware was changed for this discussion.
## 2026-10-05 - Manual common-pressure trim block

At the user's explicit request, added common_pressure_testing/commonPressureManualTrimBlock.m as paste-ready MATLAB Function block code for manual local trim search. It takes tunable absolute differential settings s_bar and a_bar, trialID, and enablePressure, and returns P1/P2/P3 target commands, phaseID, and limitFlag. The common mean is fixed at 1.6 bar. The observed near-target [1.7,1.7,1.4] bar corresponds to s=a=-0.15 bar. When disabled it outputs zero while preserving the trial label. Active target commands are clipped to a provisional 1.2-3.0 bar interval and limitFlag identifies clipped/invalid input; clipped trials should not be used for trim estimation because clipping changes the intended pressure differences. The code has no internal ramp or NDI feedback. For a physical manual experiment, route the combined command through a Simulink Rate Limiter before both the valve path and forward model and log the post-limiter command, along with the target, phaseID, limitFlag, orientation-derived lengths, signed U,V, XYZ, validity, and timestamps. A distinct trialID labels each manual hold; offline analysis should select settled plateaus using the actual logged pressure and time.

Pure MATLAB batch evaluation passed the near-target, equal 1.6 bar, disabled-zero, and clipped-input cases. This was not a Simulink compile or hardware validation. No Simulink model or NI/NDI hardware was edited or run.
## 2026-10-05 - Manual s/a search range and trial count

Discussion-only first-pass manual design for commonPressureManualTrimBlock.m: keep mean c=1.6 bar and start at s=a=-0.15 bar, corresponding to the observed near-target [1.70,1.70,1.40] bar. Probe one absolute parameter at a time with s in {-0.18,-0.15,-0.12} bar while a=-0.15, and a in {-0.18,-0.15,-0.12} while s=-0.15. These are absolute settings, not accumulated increments. The four off-center probe commands are: s=-0.12 -> [1.68,1.71,1.41]; s=-0.18 -> [1.72,1.69,1.39]; a=-0.12 -> [1.70,1.67,1.43]; a=-0.18 -> [1.70,1.73,1.37] bar. All remain above the block's provisional 1.2 bar floor and below 3 bar.

Recommended one variable-duration run with 11 labeled holds: zero reference; center; s=-0.12; center; s=-0.18; center; a=-0.12; center; a=-0.18; center; final zero reference. Give each center return a new trialID because pressure history may shift the pose. Wait for the logged post-rate-limiter pressure to settle, then hold each level about 30 s or longer if U/V still drift. Analyze final stable windows against adjacent center holds and zero references; the previous 20 s candidate hold retained about 0.0013-0.0014 rad/s bending slopes near its end. If the +/-0.03 bar probes prove too large or too small relative to observed motion/noise, adjust the subsequent search range based on the logged response. Do not treat this first pass as a final neutral-trim calibration; use it to estimate a local two-direction response and then verify a refined candidate from multiple approach histories. No code, model, or hardware was changed for this discussion.
## 2026-10-05 - Clarification of manual trim labels

In the proposed 11-hold manual search, "center" means the center of the local s/a parameter search: s=a=-0.15 bar, giving commanded [P1,P2,P3]=[1.7,1.7,1.4] bar at common mean 1.6. It does not mean the arm is straight or that the triplet is a verified neutral trim; the prior test left signed-bend error. "Zero reference" means set enablePressure=0, making the manual MATLAB Function block output [0,0,0] regardless of s and a, and hold after the downstream Rate Limiter has reached zero to measure the unpressurized arm pose. The arm's measured U,V at zero pressure need not be numerically zero. With enablePressure=1, the block computes active pressure targets from s and a; trialID is a separate manual label passed through as phaseID in either state. This is a terminology clarification only; no controller, block code, Simulink model, or hardware changed.
## 2026-10-05 - Clarified goal of the manual common-pressure search

The user clarified the intended goal: adjust the two differential settings s and a at a 1.6 bar common mean to find a pressure triplet that holds the arm straight while all three chambers remain active, then repeat that same setting from different pressure histories to assess reproducibility. The preceding offline analyses used nearby zero-pressure NDI pose as a practical reference; this pose had small nonzero signed U,V, so matching it is not exactly the same as a calibrated geometric U=V=0 straight pose. The manual search should judge straightness using both signed bending components (or independent length errors) against whichever straight reference is physically calibrated, not phi magnitude alone. Repeatability means returning to one selected pressure triplet after zero/equal-pressure or alternate loading paths and comparing settled pose, rather than merely collecting multiple unrelated adjustments. This is a discussion clarification only; no code, model, or hardware changed.
## 2026-10-05 - Read-only diagnosis of manual-trim Simulink wiring and scope

At the user's request, inspected the saved common_pressure_testing/common_pressure_trim_testing.slx read-only as an SLX archive and queried block parameters in MATLAB without starting a simulation or NI/NDI hardware. In the saved model, Dashboard Slider3/Slider4 are bound to the Value parameters of root Constant/Constant1, which feed s_bar/a_bar inputs of MATLAB Function4 containing the supplied commonPressureManualTrimBlock code. The saved slider values are s=-0.1677590521 and a=-0.1111801047; trialID and enablePressure Constants omit Value and inherit the Simulink Constant default of 1. The function's three pressure outputs feed Mux -> Rate Limiter (rising +0.8, falling -0.8 bar/s) -> both forward-model input and the pneumatic-valve input path. The model is saved in normal simulation mode with StopTime=inf. The old root P1/P2/P3 Constant blocks are commented out. No obvious disconnected pressure path was found. The function's phaseID and limitFlag outputs go to Terminator blocks; their signal-logging status was not confirmed by this archive inspection.

At the saved s/a settings, the block formula predicts target command [1.71184,1.65526,1.43290] bar. The user's 'pressures comparision' screenshot shows the first three des_pressure traces changing slightly as the sliders move and ending near those values. The scope's first input is the post-Rate-Limiter software pressure command, upstream of the pneumatic-valve Kill/saturation and NI analog outputs. The other three pressureCommand_bar traces are the inverse-model comparison output, not the manual valve command. Thus the screenshot confirms s/a tuning reaches the command path but does not establish actual chamber pressure or physical motion. Increasing the Rate Limiter to 0.8 bar/s speeds transitions but does not increase steady pressure differences: a +0.03 bar s change alters [P1,P2,P3] by [-0.02,+0.01,+0.01] bar, and a +0.03 bar a change by [0,-0.03,+0.03] bar. The next diagnostic is to compare logged post-Kill regulator-path commands and measured NDI U,V/lengths during a held slider step; if command changes but NDI does not, actual valve/chamber behavior remains unresolved without pressure sensing. No model or code was edited for this diagnosis.
## 2026-10-05 - Working pressure assumption and coarser manual trim search

The user explicitly chose to treat commanded pressure as chamber pressure for this project because chamber-pressure sensing is unavailable. Future analyses should state this as the working assumption and use logged post-limiter/post-Kill command where available, without repeatedly blocking progress on unavailable pressure sensing. The user observed that the earlier +/-0.03 bar s/a adjustments are too small for practical manual tuning. Revised discussion-only search: first reset to s=a=-0.15 at c=1.6 bar, giving [1.700,1.700,1.400] bar. Then vary s and a together to -0.10 and -0.20, returning to the starting point between holds. This preserves mean1.6 and P1=P2 while producing [1.667,1.667,1.467] and [1.733,1.733,1.333] bar. Earlier measured holds [1.65,1.65,1.50] and [1.70,1.70,1.40] showed substantial signed-bend change, so these larger joint changes should bracket the dominant response more clearly than isolated 0.03 bar changes. After identifying the better region, vary a independently by about +/-0.05 bar at fixed s to correct the remaining cross-plane bend; a shift of 0.05 changes P2 and P3 by opposite 0.05 bar while P1 stays fixed. Use logged U,V/lengths rather than visual motion alone, hold until response settles, and repeat a selected triplet from different histories. The model's current slider limits [-0.2,-0.1] permit this first coarse search. No code, Simulink model, or hardware was changed for this discussion.

## 2026-10-05 - Manual common-pressure cycling result

Inspected common_pressure_testing/manual_test_results.mat offline. It contains 635.25 s of data, four manually selected pressure triplets, 17 pressure-on episodes interleaved with 17 zero-pressure intervals, and 12706/12706 valid NDI samples under the logged orientation/sensor/configuration validity flags. The logged post-rate-limiter des_pressure and post-Kill regulator-path software commands matched exactly in this recording. Per the user's experimental working assumption, commanded pressure is treated as chamber pressure in the following interpretation. No Simulink model or hardware was run or edited, and temporary analysis artifacts were removed after inspection.

For each episode, the final 5 s of pressure-on and final 3 s of zero pressure were summarized (shorter if the interval was brief). Signed bend was reconstructed from measured lengths as U=(l2+l3)/0.013 and V=(l3-l2)/(sqrt(3)*0.013), in rad. The approximate command triplets and outcomes were: trial 1 [1.653,1.598,1.549] bar, late on episodes |(U,V)|=0.106-0.115 rad; trial 2 [1.653,1.616,1.531], late |(U,V)|=0.086-0.088; trial 3 [1.693,1.618,1.489], late |(U,V)|=0.037-0.039; trial 4 [1.666,1.649,1.485], six on episodes |(U,V)|=0.023-0.028. The first pressure-on episode at each newly selected setting differed materially from subsequent pressure cycles, so the first result alone should not be used as repeatability evidence.

Trial 4 is the best pressurized setting of the four for approximate geometric straightness and repeatability. Its selected s=-0.09840 and a=-0.08198 bar held the same command triplet through six on episodes. The five repeats after the first had mean U=0.02293 rad (sample SD 0.00103), V=0.00690 rad (SD 0.00247), and bend magnitude 0.02404 rad (SD 0.00107). The six zero-pressure intervals returned to U=0.01670 rad (SD 0.00022) and V=0.00857 rad (SD 0.00034). Relative to the immediately preceding zero-pressure interval, the five later on holds had signed-bend change magnitude 0.0047-0.0080 rad and NDI tip-position shift 3.8-4.3 mm. Zero pressure itself had a smaller absolute bend magnitude (mean 0.01877 rad) than trial 4 on pressure, so trial 4 is not a verified U=V=0 neutral trim and is not exactly the unpressurized pose. The last several trial-4 holds lasted only about 10-13 s, and some retained final-five-second U/V variation up to roughly 0.009/0.007 rad; longer holds would be needed to establish a settled equilibrium. This recording establishes promising repeatability under repeated zero-to-same-pressure cycling, not repeatability across other loading histories or an exact straight pose.

## 2026-10-06 - Local buckling observation and stability discussion

The user observed that the tip appeared straight at the identified common-pressure trim but the arm showed a small local buckle; they did not observe a global shortening of the backbone or sudden jumps while manually varying s_bar and a_bar. The location and material of the buckle have not been measured, so its mechanical origin is not yet identified. A constrained backbone can retain its axial length while carrying internal axial force; an actuator, sleeve, or other soft component can also wrinkle or kink locally without global backbone buckling. Equal-pressure nullspace in the two-coordinate ideal model therefore does not imply a stress-free configuration or unchanged distributed shape.

Smooth response during the tested manual variations is evidence against a gross snap-through on those paths, not proof that no branch switch, limit point, or local instability exists at other common-pressure levels, bend angles, load directions, or pressure histories. High bending stiffness may increase an instability threshold but does not by itself rule one out. The current constant-curvature, tip/length-based two-coordinate model cannot represent local curvature concentration or multiple centerline shapes sharing a similar tip pose. Before treating this trim as a closed-loop operating point, compare side-view whole-arm shape through slow increasing/decreasing pressure sweeps, mark buckle location/magnitude, and look for abrupt shape/pose changes, distinct loading/unloading branches, or loss of incremental stiffness. Maintain the already validated pressure envelope and stop escalation if the local deformation grows or jumps. This entry records discussion and the user's observation only; no controller, Simulink model, or hardware was changed.

## 2026-10-06 - Next common-pressure step without shape tracking

The user does not have a camera tracking system and wants to proceed with existing NDI tip/orientation-derived length logging. Treat trial 4's approximately [1.666,1.649,1.485] bar triplet (s=-0.09840, a=-0.08198 at mean 1.6 bar) as a candidate local operating point, not as a proven global straight configuration or guarantee against local buckling. First repeat zero-to-trim approaches with sufficiently long trim holds (roughly 30-40 s or until U and V settle), zero intervals between them, and at least one different pressure approach history. Compare final-window U,V/lengths, tip XYZ, commanded/post-Kill pressures, validity, and tail slopes; visually note whether the small buckle remains unchanged. The previous trial-4 on holds were often only 10-13 s and some were still moving.

If the trim remains smooth and repeatable, identify a local two-input pressure-to-signed-bend map around that triplet by perturbing s and a separately in both directions, returning to the same trim between holds. Start with small, pressure-limited perturbations and adapt their magnitude to the NDI response while retaining the validated command range. Estimate the 2x2 local Jacobian from settled differences in U,V, including loading/unloading comparison. Only then use the triplet and local response to design the common-pressure allocation/feedback test. NDI tip and length signals cannot certify the distributed arm shape, so any visible increase or sudden shift of the buckle remains a stop criterion. This is a discussion plan only; no code, model, or hardware was changed.

## 2026-10-06 - Paste-ready local common-pressure calibration block

At the user's explicit request, added common_pressure_testing/commonPressureLocalCalibrationBlock.m and its usage notes in common_pressure_testing/README.md. The single MATLAB Function block accepts Clock time t and outputs P1/P2/P3 targets plus phaseID (positive IDs for holds, negative during ramps). It schedules 15 s zero; trial-4 trim at s=-0.0983981693, a=-0.0819816714 and mean 1.6 bar; separate +0.04/-0.04 bar s perturbations and +0.04/-0.04 bar a perturbations, each followed by a return to trim; and a final zero hold. Probe holds last 35 s, intermediate trim holds 30 s, initial/final trim holds 40 s, transitions 5 s. Set Simulink Stop Time to 395 s. Pressure remains zero after 395 s. Active plateau commands range from 1.4452 to 1.6923 bar; generated smoothstep slew is below 0.5 bar/s. The existing downstream Rate Limiter and Kill path should remain connected, and the user should log phaseID and the existing pressure/NDI validity, lengths, U/V, and XYZ signals. This is an open-loop experiment, not a controller or a guarantee against the visually observed local buckle.

Offline pure-MATLAB evaluation from 0 to 400 s in 0.05 s steps passed finite and 0-3 bar bounds, active-plateau floor above 1.2 bar, zero output after 395 s, and maximum generated pressure rate 0.4996 bar/s. It printed the expected phase pressures: trim [1.6656,1.6492,1.4852]; s+ [1.6389,1.6625,1.4986]; s- [1.6923,1.6358,1.4719]; a+ [1.6656,1.6092,1.5252]; a- [1.6656,1.6892,1.4452] bar. No Simulink model was edited, compiled, or run; no NI/NDI hardware was accessed.

## 2026-10-06 - Local calibration 395 s result

Analyzed user-recorded common_pressure_testing/local_calibration_395s.mat offline without editing the model or creating analysis scripts. The saved run covers 395 s; 7900 of 7901 NDI samples passed the logged orientation, paired-sensor, configuration, and finite-length checks. Post-Rate-Limiter commanded pressure exactly matched the logged post-Kill regulator-path software commands. Pressure-to-motion interpretation retains the user's working assumption that command pressure represents chamber pressure. Each result below uses the last 10 s of a positive phaseID hold; U=(l2+l3)/0.013 and V=(l3-l2)/(sqrt(3)*0.013).

At the repeated trial-4 command [1.6656,1.6492,1.4852] bar, the five trim holds had (U,V) approximately (0.0023,-0.0216), (0.0102,-0.0240), (0.0006,-0.0232), (0.0048,-0.0131), and (0.0040,-0.0269) rad. Across them, U varied by 0.0096 rad and V by 0.0138 rad despite identical logged commands. The initial zero-pressure pose was (0.0186,0.0108) rad and final zero pose (0.0175,0.0081), close to the preceding manual run's zero reference. In contrast, the preceding manual run's later trial-4 on holds averaged (0.0229,0.0069) rad at the same pressure command. This run-to-run pressurized-pose difference is materially larger than the prior within-run repeatability and indicates an unmodeled state/history dependence under the working pressure assumption; the logs alone cannot identify its physical cause.

The s+ hold [1.6389,1.6625,1.4986] bar reached (U,V)=(0.0227,-0.0257) rad; s- [1.6923,1.6358,1.4719] reached (-0.0099,-0.0202). The a+ hold [1.6656,1.6092,1.5252] reached (0.0067,-0.0007) rad, the closest of this run to geometric U=V=0 (bend magnitude 0.0067 rad), whereas a- [1.6656,1.6892,1.4452] reached (0.0028,-0.0417) rad. Subtracting neighboring trim holds to reduce drift gives approximate signed-bend slopes J=d[U,V]/d[s,a]=[[0.397,0.069],[-0.079,0.490]] rad/bar. These directions are distinct and locally well-conditioned, but the same-command trim spread and residual final-ten-second slopes up to roughly 0.00035 rad/s limit the accuracy of J as an equilibrium map.

At 237.30 s during the constant-pressure a+ hold, the NDI-derived U,V vector changed by about 0.0063 rad in one 0.05 s sample and stayed on a slightly different orientation trajectory. Logged pressure, NDI validity, paired frame progression, and tip XYZ did not show a corresponding large change; phi remained near 0.009 rad while theta changed sharply, where bending-plane angle is sensitive near straight. This observation is not sufficient to diagnose a mechanical snap or rule one out. The a+ setting deserves a repeated, longer hold from more than one approach history, with visual attention to the known local buckle, before it is selected as the trim or used to tune feedback. No controller, model, or hardware was changed by this offline analysis.

## 2026-10-06 - Decision to try bounded common-pressure feedback

The user proposed trying closed-loop common-pressure control now, with the controller correcting imperfect static trim. This is reasonable as a limited local test because the saved 395 s calibration showed separate, well-conditioned responses to s and a perturbations around mean 1.6 bar. It does not establish that the existing PD/feedback-linearization controller will remove all equilibrium error: same-command trim pose varied within and between runs, several holds retained drift, and PD without integral or learned bias can leave steady-state error. The visually observed local buckle also remains outside the two-coordinate model.

Read-only inspection of feedback_linearization_trackingWorks.slx confirmed the current inverse allocator computes an absolute generalized force tauDesired_N, converts it using forcePerBar=pi*(0.013/2)^2*1e5, selects the minimum nonnegative effective-pressure solution, scales it to a 3 bar command ceiling, and adds the 0.8 bar deadzone only to active channels. That allocation normally leaves at least one chamber inactive. A common-pressure experiment requires replacing that allocation with an all-active, mean/trim-aware allocation; simply adding an unequal measured trim triplet to its existing output would change the delivered generalized force, while adding the 0.8 bar offset a second time would be wrong. The first closed-loop test should start from the better repeated trial-4 trim [1.6656,1.6492,1.4852] bar, wait for pressure/NDI settling, enable feedback without a command jump by initially targeting the measured pose, and then request only small static U/V changes. Bound pressure corrections to the locally tested neighborhood (s and a perturbations about +/-0.04 bar initially), preserve the existing rate limiter, active-pressure floor, 3 bar ceiling, and Kill path, and log requested/limited/post-Kill pressures, tracking error, NDI validity, and any visible buckle change. If corrections saturate or the buckle changes abruptly, stop rather than assuming the controller will compensate. The a+ [1.6656,1.6092,1.5252] point was closer to geometric straightness in one run but had only one hold with a small unexplained orientation step; keep it as a later candidate to repeat. This records a discussion/implementation direction, not a code or hardware change.

## 2026-10-06 - Paste-ready bounded common-pressure feedback test

The user authorized implementing the first common-pressure feedback test, then explicitly directed that Simulink files not be updated. Accordingly, no .slx file was copied, edited, compiled, or run. Added two paste-ready code files in common_pressure_testing/: commonPressureFeedbackInverseBlock.m replaces the root inverseModel MATLAB Function block code, and commonPressureFeedbackReferenceBlock.m replaces the continuous-circle reference block code in a user-made copy of feedback_linearization_trackingWorks.slx. The inverse block preserves the existing dynamic model, adds Clock time as a fourth input, and keeps all eight original outputs in order. The reference block takes Clock time and emits theta/phi; the user must connect the new phi output in place of the saved 1.5 rad ramp/switch reference. Set Stop Time to 240 s.

The inverse block schedules 0 bar until 10 s, ramps via an internal 0.8 bar/s per-channel limiter to the repeated manual trim [1.6656,1.6492,1.4852] bar, and holds until 50 s. At 50 s it captures only the nominal static G+Keff*q force at the measured pose, gradually engages feedback over 5 s, and allocates tauDesired minus that anchor into increments of s and a around the empirical trim, each clamped to +/-0.04 bar. The all-active mean remains 1.6 bar and no additional deadzone offset is added. Zero target begins at 220 s, reaching zero through the internal limiter before the 240 s stop. The reference is geometric zero initially, a smooth 0.01 rad U step from 90-130 s, and a smooth 0.01 rad V step from 160-200 s. This deliberately tests only a small neighborhood; it does not guarantee that imperfect trim or the visually observed local buckle will be corrected.

The implemented first test uses geometric zero as the initial desired pose, rather than capturing the measured pose as proposed in the preceding discussion. The five-second blend on incremental feedback and pressure slew limiter prevent a command jump at engagement while still allowing the controller to act on the initial trim error. This is a deliberate implementation choice, not a measured guarantee of stability.

Offline pure-MATLAB evaluation of both blocks at 0.05 s spacing passed finite 0-3 bar bounds, initial and final zero pressure, <=0.8 bar/s pressure slew, and intended reference timing. With a representative measured bend U=0.004,V=-0.027 rad, the feedback correction reduced P2 and raised P3, the direction expected to increase V; no saturation occurred in that synthetic case. These are unit-level numerical checks only. The user must paste, connect, and inspect the blocks in their Simulink model, preserve the physical Kill path, and observe logged pressure/NDI signals and visible buckle during the experiment. No NI/NDI hardware was accessed by this work.

An additional offline extreme-acceleration input saturated the local correction as intended: after settling, pressure stayed [1.6656,1.6892,1.4452] bar, mean 1.6000 bar, limited flag true, and maximum observed slew 0.8 bar/s.

## 2026-10-06 - Read-only audit of earlier tapped-delay pressure-to-length project

Read the two handoff documents in Controlling/Dynamic_jointspace_feedbacklinearization/ChatGPT and the actual dataset builder, main_NARX training/validation, dataset-generation, and offline forward-predictor code. The forward architecture includes all three pressures, causal backward pressure differences, and tanh(dP/0.05) loading indicators at delays 0 through 8, followed by three previous lengths of the particular output: 84 inputs per scalar network. main_NARX.m sets useOutputHistory=true and useOnlyActiveTrial=true and trains separate l1/l2/l3 networks with fitnet([64 32],'trainbr'); the handoff's [32 16] is an example, not the current script configuration. Training data comes from separate P1/P2/P3 20 s chirps. All three pressure channels are features, but each network is trained only on its corresponding active-actuator trial, so simultaneous-pressure coupling is not established.

main_NARX_validation.m reconstructs features using recorded output history, and predict_fwdMLP_lengths_offline_codegenMatch.m updates its output history from measured Lhat after every prediction. These are measurement-assisted predictions, not autonomous multi-step rollouts. Despite the latter function's Leq/equilibrium naming, the implemented prediction is current length L(k) from P(k) and past lengths; a physical equilibrium interpretation is unproven. For NMPC, initialize histories from real past data, compute pressure-rate/loading features from candidate future commands, and feed predicted lengths into future output-history taps without mutating the real history during optimization. Verify pressure-to-output sample alignment and the real sampling period; the predictor's default Ts fallback is 1 s and should not be relied on. Future desired references can come from known trajectories or a constant current target; future measured lengths are not required.

Before reuse in common-pressure control, the priority is held-out recursive prediction testing with simultaneous-actuation data, using two independent length outputs or another construction enforcing the backbone constraint. The current validation plots alone do not demonstrate that capability. The dataset builder smooths length labels with centered Savitzky-Golay filtering and training normalizes before block division, so validation should also audit preprocessing and avoid leakage from held-out data. No earlier-project files, controller code, models, or hardware were changed; this is a factual audit and discussion entry.

## 2026-10-06 - First bounded common-pressure feedback recording

Analyzed user-recorded common_pressure_testing/feedback_testing_results.mat offline without running hardware or editing controller/model code. The run spans 240 s at a 0.05 s NDI sample grid; 4800/4801 NDI samples passed orientation, paired-sensor, configuration, and finite-length checks. The logged Kp=[1100;1100] and Kd=[100;100] stayed constant; their sparse event timestamps were not gain changes. Logged pressureCommand_bar is the controller output; this MAT lacks post-Kill regulator-path logs. The user's working analysis assumption remains that commanded pressure represents chamber pressure.

The 40-50 s trim tail at [1.6656,1.6492,1.4852] bar had signed bend U=0.0089,V=-0.0066 rad (bend magnitude about 0.0111 rad), but U was still drifting at -0.00049 rad/s. With feedback active and zero reference at 80-90 s, command averaged [1.670,1.643,1.487] bar and measured U=0.0038,V=-0.0053 (bend magnitude about 0.0065 rad). Improvement cannot be attributed entirely to feedback because the trim was still settling. The repeated zero target after the U step (150-160 s) measured U=0.0033,V=-0.0068, and after the V step (210-220 s) U=0.0022,V=-0.0046.

During the 115-125 s U hold, desired (U,V)=(0.0100,0) rad while measured (0.0054,-0.0065), with mean length-error norm 0.112 mm. Relative to the preceding zero-target hold, U moved about +0.0016 rad toward the +0.010 step. During the 185-195 s V hold, desired (0,+0.0100) while measured (0.0016,-0.0032), with mean length-error norm 0.211 mm. Relative to the preceding zero-target return, V moved about +0.0036 rad toward the +0.010 step. Thus both responses had the intended direction but retained substantial tracking error, especially cross-plane V. Final-window U/V slopes on these two holds were small (roughly <=0.00006 rad/s), supporting a persistent rather than purely transient error over their last 10 s.

The pressure command retained mean 1.6000 bar during feedback. In the 55-220 s feedback window, relative s corrections ranged -0.00749 to +0.00288 bar and relative a corrections +0.00390 to +0.00973 bar, far inside the tested +/-0.04 bar bounds; pressureLimited was false for all 3300 feedback samples. Kp/Kd therefore did not demand enough pressure change to remove the observed error, rather than the new local pressure clamp blocking a larger command. Any gain increase or integral/static bias change should be judged against rate limits, drift, and the known local buckle; this recording alone does not establish the right gain.

At 14.40 s one NDI sensor pair became invalid for a single 0.05 s sample, orientationValid/configurationValid and inverseModelValid went false, and the logged command jumped from about [1.666,1.649,1.485] bar to zero immediately (33.3 bar/s over the sample), then ramped back at 0.8 bar/s. The paste-ready inverse block's invalid-input early return bypasses its normal internal slew limiter. This is a concrete code issue to fix before another hardware run, while retaining a fail-to-zero path for persistent invalid input. All other observed pressure ramps were at or below 0.8 bar/s. No controller or model file was changed during this result analysis.

The on-disk FBL_commonP.slx currently contains the old inverseModel allocator and thetaCircle script, although the MAT log clearly shows the 240 s common-pressure test behavior. The tested block edits therefore appear to exist only in the running/unsaved model state (or a different model state). Preserve the working Simulink edits in a saved model before closing if this run is to be reproducible; no .slx save was performed by Codex.

## 2026-10-06 - Documentation system migration and saved-model clarification

Reorganized project documentation into project_docs/PROJECT_STATE.md,
PROJECT_LOG.md, DECISIONS.md, VALIDATION_RESULTS.md, OPEN_ISSUES.md,
MODEL_ARCHITECTURE.md, and TECHNICAL_NOTES.md at the user's request.
The full prior root log was copied byte-for-byte before this entry was appended;
its SHA-256 was 31A09742DA0F9EDEF769C9733C263897B4CA066C44B0778180429C8B5A8C1EA7.
The root log is now a navigation pointer; local experiment READMEs and result
reports retain detailed operating instructions, artifacts and provenance.
Read PROJECT_STATE first and the relevant specialized document; search this
append-only history when needed. Maintain the files according to the user's
material-update rules recorded in PROJECT_STATE, preserving the discussion-first
preference and requiring explicit direction for code/model edits and hardware.

Read-only inspection of the current FBL_commonP.slx archive found the supplied
commonPressureFeedbackInverseBlock script in chart_23 and the supplied
commonPressureFeedbackReferenceBlock script in chart_146 (the chart formerly
called thetaCircle). This supersedes the prior entry's on-disk old-script finding:
the model file was modified later on October 6. Wiring, current gains, Kill state,
compilation and equivalence to feedback_testing_results.mat were not verified.
The invalid-input early return still bypasses pressure slew in the saved inverse
script and standalone source. Current offline runner settings also differ from
the common-pressure inverse: K=650 and D=20 with world-to-arm gravity rotation,
versus K=650 and D=40 with directly specified arm gravity in the inverse.
No technical code, Simulink model, recording or analysis result was modified,
no MATLAB tests rerun, and no hardware operated during migration.
