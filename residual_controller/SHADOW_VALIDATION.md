# Residual-model observation blocks for RL_dataCollection_v2

Paste-ready source only: **you add and wire these blocks manually**. No SLX was
edited or opened/run; no hardware was accessed during implementation or tests.
The outputs are diagnostics. Keep all existing pressure generation, physical
Kill, limits, filters, references and acquisition connections.

## Setup

From the project root in the MATLAB base workspace:

```matlab
addpath(fullfile(pwd,'residual_controller'));
[shadowOptions,ResidualShadowDiagnosticsBus,ResidualShadowContextBus, ...
 shadowPressureOptions,shadowInverseOptions,ResidualShadowInverseDiagnosticsBus] = ...
    setupResidualShadow(pwd);
```

Assign those exact variable names. Setup adds paths, verifies existing frozen
packet hashes and creates typed buses; it does not load/edit/run a model.
The forward block uses **shadowOptions,11x1**, not PAC's eight-entry
residualOptions. Default network is GRU128. To compare64, set
shadowOptions(1)=64 and shadowInverseOptions(1)=64 before a fresh run.

## Forward block: install first

Create a separate MATLAB Function block, e.g. `Residual model shadow`, and paste
**rlResidualShadowBlock.m**. Set its Property Inspector update method explicitly
to **Discrete, sample time0.05s**. Helper functions remain on the MATLAB path.

| Input port | Name | Type/size | Branch from |
|---|---|---|---|
| 1 | lengths3_m | double3x1 | NDI Sensor reading public output3, meas_length=[l1;l2;l3], orientation-inferred metres |
| 2 | orientationValid | Boolean scalar | NDI Sensor reading output4 |
| 3 | baseValid | Boolean scalar | NDI Sensor reading output7, baseOreintationValid |
| 4 | pollStart_s | double scalar | NDI Sensor reading output5, physical pollStart clock |
| 5 | tipFrame | double scalar | Element2 of frame, NDI Sensor reading/Sensor/MATLAB System output10; element1 is base,2 is tip |
| 6 | pressure_bar | double3x1 | Mux final P1_NI_DAQ/P2_NI_DAQ/P3_NI_DAQ after Kill/limiting, before voltage conversion inside Pneumatic Valves |
| 7 | commandTime_sim_s | double scalar | Simulation timestamp sampled with that final pressure vector at its source rate |
| 8 | t_sim_s | double scalar | Existing Reference time50ms/ZOH Clock |
| 9 | reset | Boolean scalar | Constant false, or a true pulse to reset this observer only |
| 10 | shadowOptions | double11x1 | Constant value shadowOptions |

You can bring internal frame/NI signals out through extra branch Outports or
new Goto/From tags. Preserve existing connections. Do not use desired lengths,
XYZ-fitted lengths or PAC filtered states for input1. The observer computes its
own causal velocity matching training: tau0.20s, zero initial velocity and
segmented warmup. The existing collection derivative block holds states through
invalid samples; it is therefore deliberately not reused here.

For pressure freshness, sample a simulation Clock alongside the final pressure
vector at its100Hz source. Transfer the vector and its source timestamp together
to the0.05s observer using compatible Rate Transition/ZOH handling. Do not
relabel a held pressure sample with the current observer time. Input7 and8 use
the same simulation-clock epoch. This checks software command freshness, not
physical valve latency; these are pressure commands, not chamber measurements.

Input4 must be physical pollStart_s, not validMeasurementTime or a Clock.
Host elapsed time is compared only with host time; command age is compared only
with simulation time. No comparison between these different epochs is made.

| Output port | Name | Type/size | Meaning |
|---|---|---|---|
| 1 | predictionLength_m | double2x3 | Issued now, for host time pollStart_s+0.05s |
| 2 | predictionAtMeasurement_m | double2x3 | Previously issued prediction matched to the current valid observation |
| 3 | retimedPrediction_m | double2x3 | Prior-state prediction recomputed using actual elapsed physical time |
| 4 | shadowDiagnostic | Bus: ResidualShadowDiagnosticsBus | Flags, errors, forces and timing |
| 5 | inverseContext | Bus: ResidualShadowContextBus | Immutable pre-update history for the optional inverse |

Matrices: rows[l2;l3], columns[augmented,nominal,constant velocity], metres.
Set matrix outputs to fixed2x3. In Symbols, set each bus output to its exact
`Bus: ...` type, size1 fixed. Use Bus Selectors for diagnostic fields.

For live length plots, compare shadowDiagnostic.qMeasured_m with
**predictionAtMeasurement_m(:,1)**, then columns2/3 as separate predictors.
Port1 predicts a future instant, so do not compare it directly with current q.
Use matrix Selectors for two-element columns. Reshaping2x3 to6x1 gives
[l2_aug;l3_aug;l2_nom;l3_nom;l2_cv;l3_cv]. Multiply by1000 for millimetres.

Log all three prediction outputs and these bus fields:

- qMeasured_m, velocity_mPerS, rawForce_N;
- ready, advanced, accepted, predictionValid, scoreValid, retimedScoreValid;
- errorAtMeasurement_m, retimedError_m (prediction minus measurement);
- dtPhysical_s, commandAge_sim_s, horizonMismatch_s;
- historyReset, segment, status, maximumTrainingZ;
- pressureChangeSinceOrigin_bar, predictionTarget_host_s.

Retain existing raw NDI/frame/poll/validity/diagnostics, final commands and
trajectory metadata. Full128-state inverseContext logging is optional.

## State, timing and interpretation

The frozen GRU commits once per accepted mature observation, using measured
q/dq, final commands and previous physical dt. Original nonlinear SDIRK physics
propagates with held pressure and residual force. The nominal column uses r=0;
constant velocity uses q+dq*horizon. Physics is the as-trained fixed vertical
gravity[0;0;-9.81], K/D and corrected all-three-actuator soft bounds. This is not
a newly learned rotated-base/payload model. Keep base/calibration consistent and
record base gravity to assess that assumption.

Velocity warms for0.60s, then GRU history warms for1.0s: readiness starts after
about1.6s of contiguous physical observations. Missing/invalid data, repeated
frames at a new tick, excessive clock gaps or backward clocks split the segment
and re-warm automatically. Re-executing the identical tick cannot advance or
score twice. Invalid/warm predictions are NaN, not zero. There is no latched
pressure fault, control output or automatic actuator re-arm here.

Prospective scoring accepts a previously issued0.05s prediction only if the next
physical interval matches within0.005s; inspect horizonMismatch_s. Retimed scoring
uses the same prior q/dq/p/r but the actual elapsed interval, known afterwards;
it does not use new q/dq/p as propagation inputs. Retimed accuracy is **not**
proof of prospective timing accuracy. Both assume origin pressure was held;
pressureChangeSinceOrigin_bar helps review that assumption. Rapid100Hz command
changes need an explicitly sampled command-history integrator later. Do not
interpolate missing NDI values or score held outputs as truth.

maximumTrainingZ is the largest absolute normalized feature: above5 indicates
extrapolation for review, not calibrated model uncertainty. ready denotes warm
history/finite propagation, not physical accuracy. The0.05m prediction guard
does not change the physical model's0.035m soft bounds.

Status:0reset;1invalid observation/command;2repeated tick/frame;
3velocity warmup;4GRU warmup;5ready;6new segment/gap/clock reset;
7nonfinite/out-of-guard prediction. No pressure-fault state exists here.

shadowOptions entries:

| Entry | Default | Meaning |
|---|---|---|
| 1 | 128 | GRU width64/128 |
| 2 | 0.05s | Prospective horizon |
| 3 | 0.20s | Causal velocity-filter tau |
| 4 | 0.60s | Velocity warmup |
| 5 | 1.0s | GRU warmup after velocity maturity |
| 6 | 0.10s | Maximum adjacent physical interval |
| 7 | 0.075s | Maximum adjacent simulation interval |
| 8 | 0.075s | Maximum final-command age, simulation clock |
| 9 | 0.005s | Prospective scoring timing tolerance |
| 10 | 0.05m | Numerical coordinate prediction guard |
| 11 | 3bar | Maximum accepted pressure command |

Retain the preprocessing defaults for the first comparison. Option changes
reset observer history and change the validation contract. This version is a
one-step observer; a0.1–0.5s rollout bank remains future work.

## Optional inverse observation block

For the separate **single-pass subtraction** alternative (one neural evaluation,
no nonlinear inverse solve), see [SINGLE_PASS_SHADOW.md](SINGLE_PASS_SHADOW.md).
It uses the same forward context and force target, exposes pressures before/
after slew, and has a separate diagnostic bus. Existing blocks are unchanged.

After the forward observer works, create another discrete0.05s MATLAB Function
and paste **rlResidualInverseShadowBlock.m**. It is stateless. Set inverseContext
input to `Bus: ResidualShadowContextBus`, connected from the forward block.

Inputs: inverseContext bus; desiredAcceleration_mPerS2 double2x1 (desired reduced
[l2;l3] acceleration from the desired-length derivative block);
shadowPressureOptions double10x1; variableCommon Boolean; shadowInverseOptions
double8x1. The two vectors are Constants from setup. variableCommon can be true
or false to examine the hypothetical allocator.

Outputs: residualPressurePrediction_bar double3x1;
nominalPressurePrediction_bar double3x1; inverseDiagnostic scalar fixed
`Bus: ResidualShadowInverseDiagnosticsBus`. **Only connect pressure predictions
to Scopes/logging. Never connect them to Pneumatic Valves.**

Measured q/dq plus supplied ddq define the frozen nominal LHS force. The block
solves Bp+r(q,dq,p,history)=target using the immutable pre-update history and
bounded allocation/slew, without committing any neural state. Nominal comparison
solves Bp=target. This is a model-inverse diagnostic, not PAC's adaptive/trim-
anchored demand, outer feedback or a replacement controller.

Hypothetical defaults: common1.6bar, zero trim,0..3bar, s/a caps2.5bar,
slew6bar/s, activation0 for the original linear training pressure reference.
No extra deadzone is inserted around the learned residual. These settings do
not change actual RL commands. Set them before a run if another hypothetical
allocator is desired. Network width must match the forward block.

Log valid, converged, limited, nominalLimited, iterations, targetForce_N,
predictedForce_N, inverseError_N, nominalInverseError_N, pressureDifference_bar.
Convergence is internal numerical consistency, not measured force accuracy.
Three pressures are not uniquely determined by two coordinates: exact agreement
with all three recorded commands is not an inverse-accuracy criterion. This
block adds no gravity re-expression, adaptive weights or trim calibration.

## Operator run and save

Use an existing known collection experiment and start it manually when ready.
No gain change/new excitation is needed to observe predictions. First verify
q/frame/clock/final-pressure branches and bus logging. Let the observer warm,
then view errors with validity/ready flags. Save your edited model yourself.

After a run, with out and options in the base workspace, for example:

```matlab
save('shadow_SineWave_r01.mat','out','shadowOptions', ...
     'shadowPressureOptions','shadowInverseOptions','-v7.3');
```

Use a new filename per trial. Configure signal logging before running: save
cannot recover signals that were never logged.

## Verification and remaining checks

Canonical numeric evidence: shadow_validation_20261103_v4. Exact trained-flow
and independent sequential-cell agreement, zero-residual nominal reduction,
prefix causality, wrapper identity, widths64/128,1/2/3-invalid recovery,
repeated tick/frame, clock/gap, stale command, base invalidity and fixed/variable
common inverse consistency all pass. The inverse context remains unchanged.

Saved firstTime replay:1935rows,1825warm predictions,1814prospective scores,
1822retimed scores; observer ready again at23s after recorded shutdown20.80s.
Pre-fault prospective RMSE: augmented[0.694,0.622]mm; nominal[1.537,1.747]mm;
constant velocity[0.111,0.098]mm,380samples. Most are startup/trim/blend; only
0.65s of circle control precedes the fault. This is not proof of improved
physical tracking. Post-fault mostly stationary data are scored separately.
Saved replay assumes current final-command samples, previously verified equal
to the NI path; actual source timestamp wiring still needs operator verification.
Numeric observer replay median/P95/max0.705/0.798/1.799ms, excluding actual
acquisition, Simulink scheduling and optional inverse overhead.

MATLAB Coder emits both C sources, then the standalone build stage hits this
environment's existing `Index exceeds array bounds` compiler error. This is not
a successful MATLAB/Simulink build. Independently compiled GCC components pass
1326native forward/inverse reference ticks over widths64/128, dropouts, warmup,
clock gaps, reset, stale commands and both allocation modes; see
shadow_c_replay_20261103_v2/native_results.json. Actual model build/timing and
physical predictive accuracy remain operator checks. No physical run occurred.

Reproduce in fresh directories using validateResidualShadow,
checkResidualShadowCodegen, writeResidualShadowCReplay and
checkNativeResidualShadow.py. Source/fixture versions, generated C and SHA
inventories are preserved with this review.
