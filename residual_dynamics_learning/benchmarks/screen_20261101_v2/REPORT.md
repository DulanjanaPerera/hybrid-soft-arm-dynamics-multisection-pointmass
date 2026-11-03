# Step3 preliminary offline benchmarks

2026-11-01. Canonical screen version2; no Simulink/hardware operation or
controller integration. Fitted models and transforms are in each task's
frozen_selection.mat. They are research predictors, not deployment artifacts.

## Protocol

Dataset two_day_v2, separate forward and inverse models. October19 development:
11,785 training anchors after deterministic stride3 per run; all13,977 validation
anchors. Both neural candidates use one seed73, CPU, at most20epochs, batch256,
Adam.001, L2=.0001, validation patience4. MLP32/16tanh, GRU16units, history10
frames with reset per window. No stateful long-history training or architecture
search. NN best epoch uses normalized pooled validation loss; final candidate
selection uses mean per-run vector RMSE. Ridge lambda grid[1e-6,.001,.1],
selected by validation;1e-6 chosen for both tasks. Stable double augmented-QR
solve exactly reproduces stored weights. Normalizers/targets/physics-context
transforms fitted only on training data.

Physics-assisted forward baseline is the exact local affine transition with
M,C,K,G frozen at current estimated state, command held over observed dt. It
is NOT exact nonlinear rollout. Learned output is next-position residual inmm;
velocity is filtered/inferred, not true measured mechanical velocity. Inverse
learns delta_d inbar; reconstruct three pressures with prescribed p0. Inputs
include nominal d and motion preview. Pressure history ends atk-1. At execution
preview must be planned motion and p0 externally given. Current experiment's
actual common pressure is used to specify that conditional offline task.

Each main task's weights/selection saved before loading October27 sequences.
All24,079 eligible external anchors evaluated with frozen candidates. Small
external step46/combined26 subsets remain selection-biased; macro metrics give
those tiny runs equal weight, so also inspect pooled and per-run CSVs. No
controller tracking scores or measured-force accuracy are reported here.

Screen v1 was superseded after a conditioning warning exposed single-precision
normal equations in ridge fitting. Same protocol rerun with double QR, without
external-informed hyperparameter changes. Validation-motivated shortcut
ablations added with their own saved protocol/frozen models before their
external evaluation. October27 has now been evaluated and must not be described
as a permanently untouched confirmation set for future iterative changes.

## Forward next-position prediction

Mean per-run vector RMSE, mm. One step is one native sensor interval, not
necessarily exactly50ms physical time.

| Predictor | Oct19 validation | Oct27 external |
|---|---:|---:|
| Frozen nominal mechanics | 4.98454 | 5.42676 |
| Position persistence | .10352 | .08861 |
| Constant velocity | .06128 | .05193 |
| Nominal + linear history residual | .05802 | .08939 |
| Nominal + MLP residual | .41469 | 2.47972 |
| Nominal + GRU residual | .30156 | 2.80554 |
| Kinematic baseline + ridge ablation | .05807 | .04887 |

Linear physics residual wins initial validation among main learned models,
but its cross-day aggregate is worse than constant velocity. It improves the
larger circle/helix external sets slightly and degrades pressure excitations,
especially where physical dt differs. Kinematic ablation is nearly identical
on validation and better externally: usefulness of nominal mechanics as a
learning prior is not established by this screen. Comparing only against the
poor nominal baseline would overstate learned dynamics improvement.

Free-running validation:20initial windows per run,100per model,10steps. Initial
q/v history observed; later q/v updated from prediction only. Pressure sequence
and query times are known. No true future state resets or per-step teacher
forcing. Gap/split guards retained. Horizon about.5s on day1, not a full circle.

| Predictor | Macro path vector RMSE mm | Macro final vector RMSE mm |
|---|---:|---:|
| Nominal | 7.42 | 8.20 |
| Existing full nonlinear RHS (ode15s) | 6.69 | 6.90 |
| Persistence | .50 | .76 |
| Constant velocity | .23 | .37 |
| Linear residual | .21 | .33 |
| MLP residual | 2.35 | 3.10 |
| GRU residual | 1.32 | 1.59 |

No trajectory exceeded the offline.05m coordinate-domain guard in these short
rollouts. This is not a stability/safety guarantee. Complete per-run figures
and any failure counts are in forward/validation_rollouts.csv.

The full nonlinear baseline separately integrates armS_single_dynamics with
ode15s, retains its integrated dq, and never resets from observed future states.
Same100initial windows, inputs and10step horizons; all completed without domain
guard failure. See nonlinear_nominal_rollouts.csv. This checks the actual
recursive mechanical model in addition to the frozen-coefficient predictor.
The learned forward networks here correct POSITION transitions, not physical
acceleration/force; do not inject their outputs directly into the ODE RHS.

## Conditional inverse pressure prediction

Mean per-run three-pressure vector RMSE, bar; p0 externally conditioned.

| Predictor | Oct19 validation | Oct27 external |
|---|---:|---:|
| Nominal inverse | .50588 | .46741 |
| Previous differential-pressure baseline | .00985 | .00906 |
| Nominal + full ridge residual | .000614 | .000626 |
| Nominal + MLP residual | .03613 | .02779 |
| Nominal + GRU residual | .03294 | .02764 |
| Pressure-history + common-pressure AR | .000621 | .000424 |
| Motion/preview ridge, no pressure history | .32888 | .40082 |

The AR ablation sees past pressure and prescribed common pressure but no motion
or nominal dynamics. It matches/beats full inverse ridge. Removing pressure
history makes performance much worse. Thus tiny inverse pressure errors do
NOT establish learning of missing mechanical dynamics or successful inverse
control of a new desired motion. Smooth/repetitive pressure schedules and
the supplied common-pressure condition permit an easy prediction shortcut.

No output clipping used to conceal infeasibility. Inverse ridge worst external
bound violation.00166bar; none exceed.01bar. MLP/GRU worst.07308/.08071bar,
with19.86%/28.00% external examples violating bounds by>.01bar. The initial
summary's strict1e-8 violation fraction includes tiny negative roundoff-like
predictions nearzero; inverse_constraints.csv adds magnitude-aware figures.

## Conclusion and next action

Preliminary comparison completed; no universal best architecture identified.
Do not deploy or claim controller improvement from these models. First-screen
failure does not prove recurrent networks unsuitable: memory length, capacity,
loss, seeds, initialization and data coverage have not been broadly tuned.

Next useful research work:

1. Revisit nominal calibration and a constrained actuator-memory baseline;
   avoid asking the residual to cancel most of an inaccurate prior.
2. Specify a motion-conditioned inverse validation that requires novel desired
   motion, not just extrapolation of smooth command traces. Retain AR baselines.
3. Address physical timing differences and get suitably varied excitation/
   speeds before longer model-memory/architecture searches if needed.
4. Use fresh independent confirmation data before any performance claim;
   later arm tests remain operator-run and require a qualified frozen model.

Artifacts: protocol.mat, fitted task banks, training histories/times, per-run
validation/external CSVs, rollout logs, ablation weights/reports,
ridge_reproduction.csv, inverse_constraints.csv and source_snapshot.
validateResidualScreening passed train-only statistics, target-free inference,
exact ridge reproduction and pressure-bound audit. No commits/push performed.
