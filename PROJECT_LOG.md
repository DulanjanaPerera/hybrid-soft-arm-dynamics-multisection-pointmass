# Project research log

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
