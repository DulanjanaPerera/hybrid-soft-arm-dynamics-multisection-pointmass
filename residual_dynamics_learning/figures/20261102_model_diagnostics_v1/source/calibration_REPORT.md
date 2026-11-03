# Offline mechanical calibration and effective actuator memory

Completed 2026-11-02 on codex/residual-dynamics-learning. This is an offline
model-identification experiment. No PAC/PDFB settings, SLX files or hardware
were changed. No fitted model is qualified for deployment.

## Main finding

Constrained mechanical calibration improves predictions relative to the original
mechanical model. A shared first-order command-memory extension adds only a
small improvement. Both are still worse than matched persistence and constant-
velocity baselines. High fitted damping and the lag/delay tradeoff prevent a
physical-parameter interpretation. Do not copy these values into PAC tuning.

## Protocol and data

- Prepared input: datasets/two_day_v2, backed by the audited raw-source hashes.
  Only the October 19 cohort is loaded here. October 27 is not used for fitting,
  model selection or evaluation in this calibration; it has already been
  evaluated in the earlier neural screening and is not fresh confirmation.
- Static fit uses 15 StepWave **training** holds, one median per trial. Eligible
  holds are mature valid phase-2 samples, final two physical seconds, at least
  15 samples, command range <=0.01 bar and median speed <=0.0005 m/s per axis.
  Damping is absent from this equilibrium fit. Additional holds are diagnostics.
- Dynamic fit uses ten active training windows, two per waveform, 20 steps
  (about one physical second). Validation uses 30 windows, six per waveform,
  40 steps (about two seconds). Whole training/validation trials are disjoint.
  Windows stay within valid mature segments. Selection requires a pressure
  vector change >=0.05 bar within 20 frames and spaces eligible anchors evenly.
  These small selected subsets do not represent every trajectory segment.
- Initial q and causally filtered dq come from the first frame. Subsequent q/dq
  are integrated without measured-state resets. Known recorded pressure commands
  are inputs, not predicted controls. Scoring uses observed future q only.
- Native command timestamps are mapped from simulation time to sensor physical
  time by interpolating the two recorded clocks. No q samples are interpolated
  across gaps. Actual NI-write wall times and chamber pressures were not logged;
  this clock mapping is approximate and limits physical-delay identification.
- The existing armS_single_dynamics and armS_single_core are integrated with
  ode15s, including nonlinear M/C/G, structured K, diagonal D and bound forces.
  Validation uses tighter tolerances and a 0.01 s maximum step. A 45 mm per-axis
  domain event rejects incomplete rollouts; failed candidates cannot win.

## What was fitted

Fixed: q=[l2;l3], nominal geometry L=0.17411 m and r=0.013 m, mass=0.1 kg,
gravity, bound force, and alpha=pi*(0.0065 m)^2*1e5 N/bar with
A=[-1 1 0;-1 0 1]. No activation offset or learned force map is added here.

First fit positive structured stiffness and constant generalized-force bias
using bounded linear least squares. Then freeze them and fit diagonal damping
using integrated training output error. Finally freeze mechanics and fit one
shared first-order effective command state, T*zdot=p(t-delay)-z, at each fixed
delay in [0,0.1,0.2] s. Memory is initialized from preceding recorded commands;
it is not measured chamber pressure. Candidate parameters are saved before
validation. This is a staged, limited-budget fit, not joint global optimization.

| Quantity | Original | Calibrated | Fit bounds |
|---|---:|---:|---|
| k1 (N/m) | 650 | 1349.866 | 1 to 6000 |
| k2 (N/m) | 650 | 630.474 | 1 to 6000 |
| k3 (N/m) | 650 | 1008.796 | 1 to 6000 |
| d2 (N s/m) | 40 | 250 | 0.5 to 250 |
| d3 (N s/m) | 40 | 250 | 0.5 to 250 |
| b2 (N) | 0 | 0.810708 | -8 to 8 |
| b3 (N) | 0 | 1.821780 | -8 to 8 |

Both damping estimates reach the upper bound. A separate **training-only**
profile with shared D=40/100/250/500/1000 gives costs
5.3833/4.9337/3.9428/2.7434/1.5689. This profile does not refit parameters or
evaluate expanded bounds on validation. Increasing damping suppresses predicted
motion; with persistence outperforming mechanics, these results are consistent
with compensation for other mismatch, not unique identification of damping.

Memory fit: delay 0/0.1/0.2 s gives T=0.11242/0.04330/0.01066 s, with nearly
flat training costs 3.87885/3.85730/3.86537. Mechanical-only cost is 3.94280.
Validation selects delay=0.2 s among these candidates, at the grid boundary,
with T near its 0.01 s lower bound and below the roughly 0.05 s sensor interval.
Delay and T trade off; do not report a physically identified 0.2 s actuator
delay or a resolved 10.7 ms time constant.

## Held-out results

Metric below is the mean of each window's vector path RMSE:
sqrt(mean((q2_pred-q2_obs)^2+(q3_pred-q3_obs)^2)), in mm, excluding the initial
frame. It is neither per-coordinate tracking RMSE nor whole-arm XYZ error.
Each waveform contributes six windows. Overlapping windows are not independent
experimental repeats. All model rollouts completed without a domain failure.

| Predictor | Path RMSE over ~2 s (mm) | Mean endpoint error (mm) |
|---|---:|---:|
| Original nonlinear model | 5.136 | 5.339 |
| Calibrated mechanics | 3.369 | 3.905 |
| Calibrated + memory, delay 0 | 3.324 | 3.853 |
| Calibrated + memory, delay 0.1 s | 3.314 | 3.840 |
| Calibrated + memory, delay 0.2 s | 3.295 | 3.813 |
| Persistence | 1.338 | 2.346 |
| Constant velocity | 1.055 | 1.944 |

Mechanics reduces path error about 34%; memory adds about 2% relative to
calibrated mechanics. Constant velocity remains better at every checked
horizon: ~0.5/1/2 s path errors are 0.154/0.372/1.055 mm, versus
2.295/2.791/3.295 mm for the selected mechanical-memory candidate.

Per-waveform ~2 s original -> calibrated -> calibrated+selected-memory mm:
circleSweep 7.81 -> 5.50 -> 5.35; combinePressure 1.85 -> 1.39 -> 1.37;
helix 6.98 -> 4.04 -> 4.03; SineWave 4.09 -> 1.81 -> 1.77;
StepWave 4.95 -> 4.11 -> 3.96.

Independent static-equilibrium solves start from zero, not observed validation
q. Mean vector error across four validation holds improves 7.604 -> 5.583 mm;
19 diagnostic training holds improve 9.343 -> 5.638 mm. These residual errors
are substantial even at near-zero speed. The 15 fitting holds are loading
holds; this does not identify quasi-static unloading hysteresis separately.

## Interpretation and next work

Priority update after discussion2026-11-02: the nominal model's inaccuracy was
already known; this was a limited parameter/memory diagnostic. The calibration
workstream is now paused in favor of missing-dynamics residual learning around
a frozen nominal baseline. The diagnostic suggestions below are retained as
possible targeted checks, not prerequisites to make the nominal model accurate.
See the latest RESIDUAL_DYNAMICS_LEARNING_PLAN.md and FIGURE_INDEX.md. Figures
and unfavorable comparisons are retained for transparent research reporting.

There is an equilibrium mismatch as well as dynamic mismatch. A pure lag model
cannot remove steady-state error under a constant input. Possible contributors
include the effective command-to-force map, pressure-dependent stiffness,
loading/unloading memory, reference/calibration geometry and the reduced
coordinate assumption. This experiment does not distinguish those causes.

Next, use the existing data to inspect equilibrium residuals against individual
and common pressure and loading direction, and audit command/state clock
alignment. Keep mechanical and actuation changes separate so they cannot all
compensate arbitrarily. If coverage is insufficient, design operator-run
stepped differential-pressure loading/unloading holds at several common
pressures, with synchronized command/sensor timestamps and chamber-pressure
measurement if available. No experiment is queued or run here.

Only resume flexible residual learning after the effective low-order model
has been evaluated on informative motion with matched baselines. A strong
constant-velocity score is a sanity check, not a replacement forward/inverse
physical model. Require longer autonomous rollouts, static and dynamic checks,
and fresh independent confirmation before PDFB/PAC integration. Do not simply
raise the damping cap or neural capacity based on this validation result.

## Reproduction and artifacts

From repository root in MATLAB, with a NEW output folder:

```matlab
addpath('residual_dynamics_learning');
validateGreyboxModel();
dataset=fullfile(pwd,'residual_dynamics_learning','datasets','two_day_v2');
folder=fullfile(pwd,'residual_dynamics_learning','calibration','new_version');
calibrateResidualGreybox(dataset,folder);
assessGreyboxCalibration(folder);
```

Requires Optimization Toolbox. No Simulink model execution is involved.
validateGreyboxModel passed independent original-RHS equivalence, future-q
invariance and delayed-input causality. assessGreyboxCalibration passed
whole-trial train/validation isolation. verification.txt records MATLAB version
and dependency paths. source_snapshot contains the exact fit/simulation/check
sources and mechanical dependencies; source_hashes.csv records their hashes.

protocol.mat and identification_data.mat preserve holds, selected windows,
initial states and commands; mechanical_fit.mat, frozen_candidates.mat and
calibration.mat preserve fit stages and selection. CSVs retain individual
static errors, dynamic horizons, baselines, training memory costs and damping
profile. example_validation_rollouts.png shows one fixed example per waveform,
not all validation windows. horizon_summary.csv is a direct aggregation of the
saved per-window CSVs. The incomplete greybox_20261102_v1 is explicitly marked
INCOMPLETE and must not be used. No raw files or prior benchmark outputs changed.
