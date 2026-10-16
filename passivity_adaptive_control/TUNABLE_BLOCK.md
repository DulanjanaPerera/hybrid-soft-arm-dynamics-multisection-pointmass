# Tunable PAC with switchable common pressure

Use **pacAdaptiveTunableBlock.m** in a duplicate MATLAB Function block/model.
The existing pacAdaptiveCommonPressureBlock and its ports remain usable.
No SLX file was edited. Set the new block explicitly to discrete **0.05 s**.

The common-pressure switch controls bounded allocation, not a gradient parameter
estimator. It chooses a mean pressure near the nominal 1.6 bar as needed to apply
the requested two-coordinate force. It does not add mechanical coordinates.

## Setup

Run once in the MATLAB base workspace (add this folder to the path first):

```matlab
[adaptationGamma,parameterScale,parameterMaxRate,parameterLower, ...
 parameterUpper,parameterInitial,adaptationOptions,historyOptions, ...
 pressureOptions,supervisorOptions] = pacTunableDefaults();
```

Paste pacAdaptiveTunableBlock.m into the new MATLAB Function block. Retain its
existing first 11 input connections, including lambda [40;40], dampingGain
[20;20], and the Boolean enableAdaptation input. Add the following ports in order.
Use separate Constant blocks whose values are the workspace variable names;
they can later be replaced with tunable sources. Vectors are double columns.

| Input | Size | Meaning |
|---|---|---|
| enableCommonPressureAdaptation | Boolean scalar | true: variable common allocation; false: current fixed-mean allocator |
| adaptationGamma | 10x1 | Parameter adaptation rates |
| parameterScale | 10x1 | Parameter normalization scales |
| parameterMaxRate | 10x1 | Absolute parameter change caps per second |
| parameterLower, parameterUpper | each 10x1 | Projection bounds |
| parameterInitial | 10x1 | Initial estimates, read only on a fresh run |
| adaptationOptions | 3x1 | [s-norm deadband (m/s); force-residual tolerance (N); maximum update interval (s)] |
| historyOptions | 4x1 | [history length (m); beta; gamma; state bound] |
| pressureOptions | 10x1 | [nominal common; s trim; a trim; minimum; maximum; delta-s cap; delta-a cap; slew (bar/s); activation; force-per-bar (N/bar)] |
| supervisorOptions | 7x1 | [pressure start; feedback start; blend duration; invalid grace; maximum sample age; anchor speed (m/s); stop time] (other entries seconds) |

Parameter order throughout: **k1 k2 k3 d2 d3 h1 h2 h3 b2 b3**.
Defaults reproduce the current configuration. Gamma acts together with scale
as `scale.^2 .* gamma`; changing both multiplies their effects.
The history alpha remains fixed at 1 in the normalized state law; the three
history force weights adapt. Mass, geometry and gravity model coefficients stay
in the mechanical model; they are not feedback tuning gains.

Configure structural settings before a run: initial values, scales, bounds,
history law, calibration/pressure-map constants and supervisor timings. Changing
them mid-run can invalidate the retained history or trim calibration. Bounds
are not an online reset input. Restart the simulation to reinitialize estimates.
Lambda, Ks, gamma, rate caps, deadband/tolerance and the two Boolean switches
are external for later controlled tuning. Abrupt live gain changes can cause
force jumps; the existing pressure slew limiter remains the final command cap.

## Outputs and logging

First **26 outputs are unchanged**. New outputs:

27. `commonPressure_bar`: actual mean of the three command pressures.
28. `commonPressureShift_bar`: command mean minus nominal common pressure.

Log these, pressureCommand_bar, tauRequested_N, tauEquivalent_N,
allocationResidual_N, limited, adaptationActive, parameterEstimate and both enable
inputs. Log the external tuning inputs too for reproducible comparisons.
status remains uint8; validity/limited/adaptation/bound outputs remain Boolean.

The common switch can change during control. OFF returns allocation to the
configured fixed mean; it does not freeze the last common value. Switching does
not reset parameter estimates or recapture the trim force offset. The supervisor
slews all three commands at the existing configured rate, so transient residual
can remain and will gate parameter adaptation. During startup, OFF/trim/fault
behavior is unchanged; variable common pressure operates only on control targets.

## Allocation and limits

ON uses the all-active nominal branch `p >= max(minimum, activation)`. For equal
force coefficients, pressures differences span a convex polygon defined by rails
and the trim-relative symmetric/antisymmetric correction caps. The allocator
projects the requested force onto that polygon (minimum Euclidean force error),
then chooses a feasible mean nearest the nominal mean. It uses bounded edge and
vertex enumeration, not an Optimization Toolbox solver. No pressure rail or
correction cap is relaxed. It requires no integral common-pressure state.

This deliberately excludes inactive pressures below activation in ON mode. Such
pressures supply no additional force range in the nominal map; the all-active
restriction may be conservative when interpreting correction caps in pressure
coordinates. OFF mode retains the existing activation-aware fixed-mean path.
The same 0.8-bar activation assumption is used for allocation and evaluation.
It is not a chamber-pressure measurement or a newly identified actuator law.

Changing common pressure can change physical stiffness/preload, which the ideal
nullspace map does not predict. Offline exact force balance is not a hardware
tracking or passivity guarantee. Keep lambda=40, Ks=20 and the previous 30-second
circle initially; compare common OFF and ON with the same initialization/reference.

## Offline validation

`validatePacTunable` checks feasible/infeasible forces, correction/pressure bounds,
projection optimality against feasible vertices, OFF-mode equivalence and command
slew during switching. Also run validatePacPressure and validatePacAdaptive.

Alternative to ten Constant blocks: paste pacTuningBlock.m into a single MATLAB Function block and connect its ten named outputs to the matching PAC ports. Edit its explicit defaults before a run. No workspace setup is required for this source block.
