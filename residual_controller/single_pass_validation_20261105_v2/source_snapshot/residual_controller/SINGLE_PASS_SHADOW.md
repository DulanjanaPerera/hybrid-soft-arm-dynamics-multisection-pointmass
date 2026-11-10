# Single-pass residual subtraction: passive inverse observation

Paste-ready source: **rlResidualSinglePassShadowBlock.m**. Add a separate
MATLAB Function block manually in your existing `residual_FWD_INV_validation`
model. Keep the forward observer and existing iterative inverse for comparison.
No SLX is edited by these helpers; no model-modification script is supplied.

This block has no persistent state, inverse iteration or dynamics integration.
It calls the frozen GRU once at actual current pressure, discards its returned
hidden state, and subtracts its two generalized forces from the nominal demand:

```
rFrozen = GRU(qMeasured, vMeasured, pActual, previousDt, historyPre)
target  = M(qMeasured)*ddqDesired + (C(qMeasured,vMeasured)+D)*vMeasured
          + K*qMeasured + boundForce(qMeasured) + G(qMeasured)
correctedRequest = target - rFrozen
pBeforeSlew = boundedPressureAllocation(correctedRequest)
pPrediction = clip(pBeforeSlew, pActual +/- slewRate*horizon)
```

The frozen residual has the sign used in training: nominal EoM LHS = B*p+r.
The nominal comparator uses the same allocation with r=0. Fixed/variable common
pressure, rails, correction caps and slew are retained. The residual is not
necessarily an additive positive pressure on all three actuators.

**This uses measured q/v plus desired acceleration, not desired q/dq throughout.**
It matches the existing iterative shadow block's force target. No desired-length
error feedback, integral action, PAC adaptation, gravity rotation or payload
estimator is introduced. This is an approximate inverse-model observation,
not a complete controller or a newly trained direct pressure-output network.

## Setup while the model is stopped

Keep the existing assignments from `setupResidualShadow`. If starting a fresh
MATLAB workspace, from the project root:

```matlab
addpath(fullfile(pwd,'residual_controller'));
[shadowOptions,ResidualShadowDiagnosticsBus,ResidualShadowContextBus, ...
 shadowPressureOptions,shadowInverseOptions,ResidualShadowInverseDiagnosticsBus] = ...
    setupResidualShadow(pwd);
```

Then add the new options and bus (also works in an already configured workspace):

```matlab
[singlePassOptions,ResidualSinglePassDiagnosticsBus] = ...
    setupResidualSinglePassShadow(pwd);
singlePassOptions(1) = shadowOptions(1); % Same frozen network as forward observer
```

Assign that **exact bus variable name** in the MATLAB base workspace. The new
setup verifies packet hashes, adds paths and creates the bus. It does not
overwrite your existing base-workspace options or reset a live observer. The
new two-entry vector is `[128;0.02]`: network width and diagnostic frozen-force
balance tolerance in N. The tolerance is not a pressure limit or solver setting.
Do not use eight-entry shadowInverseOptions at this block's last port.

## Manual ports and wiring

Paste the complete `rlResidualSinglePassShadowBlock.m` file into the new MATLAB
Function block. Helper files stay on the MATLAB path. Set Update method to
**Discrete**, sample time **0.05 s**, matching the existing forward observer.

| Input | Name | Fixed type/size | Connect to |
|---|---|---|---|
| 1 | inverseContext | Bus: ResidualShadowContextBus, scalar | Branch forward observer output5, same tick as existing inverse |
| 2 | desiredAcceleration_mPerS2 | double 2x1 | Same desired reduced `[ddl2;ddl3]` m/s^2 line as existing inverse; saved collection source is MATLAB Function2 output3 (`dd_length`) |
| 3 | shadowPressureOptions | double 10x1 | Branch the existing inverse's pressure-options Constant |
| 4 | variableCommon | Boolean scalar | Branch the existing inverse's common-pressure Boolean |
| 5 | singlePassOptions | double 2x1 | New Constant with value `singlePassOptions` |

The context already contains measured `[l2;l3;dl2;dl3]`, current final pressure
commands, physical previous dt, prediction horizon and **pre-update** GRU
history. Do not add measured-length, pressure or Clock ports. Do not delay the
context or relabel held measurements as fresh. The forward observer owns all
validity, missing/repeated-sample handling, warmup and neural state commits.
The stateless inverse cannot independently detect a deliberately delayed bus.

| Output | Name | Fixed type/size | Connect to |
|---|---|---|---|
| 1 | singlePassPressurePrediction_bar | double 3x1 | Scope/logging only, channels P1/P2/P3 |
| 2 | nominalPressurePrediction_bar | double 3x1 | Scope/logging only; name logged line `singlePassNominalPressure_bar` to distinguish duplicate nominal outputs |
| 3 | singlePassDiagnostic | Bus: ResidualSinglePassDiagnosticsBus, scalar | Bus Selector and logging |

In the block Symbols/Property Inspector, set the bus input/output to their exact
`Bus: ...` types, fixed size1. Set numeric ports to the listed fixed sizes.
Leave existing forward, commands, trajectory, acquisition, feedback, filters and
physical Kill wiring in place. **Neither pressure output may drive valves or a
controller input in this observation test.** Pre-slew diagnostic outputs are
hypothetical and likewise must never drive valves.

## Bus Selector: fields to display and log

First compare these three 3x1 signals on a Scope:

1. `correctedPressureBeforeSlew_bar`;
2. `nominalPressureBeforeSlew_bar`;
3. `residualInputPressure_bar` (the actual commands used as neural input).

Then compare output1/output2 on another Scope. This separates the learned
correction from its masking by slew clipping. Pre-slew values still respect
allocation rails/correction caps; they are not unconstrained inverse solutions.
The slew box is reanchored to actual current pressure each tick, as before;
it does not accumulate the preceding hypothetical prediction.

Log the complete new bus and outputs, plus existing observer flags/context,
raw NDI diagnostics, actual final pressure commands, trajectory and the actual
option vectors. Important fields:

| Fields | Meaning |
|---|---|
| valid | Finite output from a ready current observer context; not physical accuracy |
| networkEvaluations | 1 after the single GRU call; 0 for inputs rejected before that call |
| targetForce_N | Nominal EoM demand before correction |
| residualForce_N | GRU output at actual current pressure |
| correctedForceRequest_N | target minus residualForce |
| allocationLimited / slewLimited | Separate reasons for corrected-pressure limiting |
| nominalAllocationLimited / nominalSlewLimited | Same separation for comparator |
| pressureDifferenceBeforeSlew_bar / pressureDifference_bar | Corrected minus nominal, before/after slew |
| allocationError_N | B*pBeforeSlew minus corrected force request |
| slewError_N | B*(pAfterSlew-pBeforeSlew) |
| frozenBalanceError_N | B*pAfterSlew+rFrozen-target; equals allocationError+slewError |
| frozenBalanceSatisfied | Norm of frozenBalanceError <= tolerance; **not** nonlinear inverse convergence |
| nominalBalanceError_N | B*pNominalAfterSlew-target |
| frozenTotalForce_N | B*pAfterSlew+rFrozen, not measured force |
| sampleTime_host_s / horizon_s | Copied physical observation time and hypothetical slew horizon |

No `converged` field is supplied because there is no nonlinear inverse solve.
The network depends on pressure: r(pActual) can differ from r(pPrediction).
Even zero frozenBalanceError does not prove
`B*pPrediction+r(q,v,pPrediction,historyPre)=target`. That pressure-dependence
approximation is checked separately offline, with an additional neural call
outside this real-time block. No new interpolation or readiness bypass is used.
Invalid/cold contexts produce NaN pressures and valid=false, never zero commands.

Save to a **new filename**, for example after your manually started trial:

```matlab
save('circlesweep_singlePass_shadow_r01.mat','out','shadowOptions', ...
     'shadowPressureOptions','shadowInverseOptions','singlePassOptions','-v7.3');
```

## Offline verification

`validateResidualSinglePassShadow` checks both widths, known feasible force
targets/subtraction sign, deterministic repeated calls, immutable contexts,
invalid handling, rails/slew and error decomposition. It then replays the saved
circle/helix contexts, compares targets/residuals/nominal outputs with the old
logs, and quantifies frozen versus candidate-pressure model mismatch separately.
It saves per-sample CSV/MAT, summaries and pre/post-slew PNG/PDF/FIG plots.
`checkResidualSinglePassShadowCodegen` attempts fixed-size C source emission only.
Neither helper starts a Simulink model or connects a device.

Canonical results and exact verification limitations are recorded in
`single_pass_validation_20261105_v2/README.md`. MATLAB numeric/source-generation
results do not establish successful compilation in your manually edited SLX,
full real-time timing, physical inverse accuracy or improved closed-loop tracking.
