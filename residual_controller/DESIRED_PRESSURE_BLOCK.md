# Desired-motion pressure prediction with nominal + residual model

For smooth command selection after the desired-pressure block, use the new
residualCommandValid.m pressure-handoff function and
[COMMAND_HANDOFF.md](COMMAND_HANDOFF.md). It replaces the Boolean helper
and source Switch while retaining the external command Unit Delay/Kill path.

For the operator's later direct-actuation wiring, see
[ACTUATION_CAUSALITY.md](ACTUATION_CAUSALITY.md). The raw prediction-to-valves
connection closes an algebraic loop through actual-command feedback and the
persistent observer. A validity/startup selector plus explicit0.05s command
delay is required for the described causal plumbing; this document's original
installation remains observation-only and does not qualify closed-loop control.

2026-11-09. Paste-ready source: **rlResidualDesiredPressureBlock.m**.
This replaces the observed-motion target when comparing pressure for a desired
trajectory. Existing historical block sources are retained for log replay.
No SLX is edited by these functions or instructions.

The two branches now share the same desired-state nominal force demand:

```
desiredMotion = [l2d; l3d; dl2d; dl3d; ddl2d; ddl3d];
tauDesired = M0(qd)*ad + (C0(qd,vd)+D0)*vd + K0*qd + bound(qd) + G0(qd);
nominalPressure = allocate(tauDesired);
correctedPressure = allocate(tauDesired - residualActual);
```

The learned model supplies residualActual from actual measured q/v, actual
commanded pressure, previous physical sample interval and actual causal GRU
history. No desired values are substituted into that recurrent observation
stream. One GRU evaluation is made; no history is committed by this block.

This is desired-trajectory feedforward with a frozen observed residual. The
residual also depends on pressure, so this is a single-pass approximation to
the learned inverse, not a pressure-consistent nonlinear inverse solve or a
guarantee that the arm reaches the requested lengths. PAC integration later
requires its complete feedback/adaptive force demand, rather than substituting
this collection feedforward block for PAC.

## Manual replacement in residual_FWD_INV_validation

Stop the model first. Copy the **entire** rlResidualDesiredPressureBlock.m,
including its local minimumPressure function, into the existing single-pass
MATLAB Function block (MATLAB Function8 in test4). Set it to discrete0.05s.
Keep all outputs on Scopes/logging for this installation. Preserve command,
reference, filter, acquisition and physical Kill wiring.

The number of ports remains5inputs/3outputs. **Input2 changes from2x1 desired
acceleration to6x1 desiredMotion.** Do not leave the old acceleration-only
connection or supply three chamber lengths as the first three elements.

| Port | Name | Size/type | Connection |
|---|---|---|---|
| In1 | inverseContext | scalar Bus: ResidualShadowContextBus | Same forward observer output5; no added delay |
| In2 | desiredMotion | double[6 1] | Branch the existing6-vector entering collection inverse MATLAB Function3 |
| In3 | shadowPressureOptions | double[10 1] | Existing pressure-options Constant |
| In4 | variableCommon | Boolean scalar | Existing connection; retained but ignored by minimum-pressure policy |
| In5 | singlePassOptions | double[2 1] | Existing Constant; width must match forward observer |
| Out1 | desiredPressurePrediction_bar | double[3 1] | Corrected desired pressure AFTER pressure/slew limits; Scope/log only |
| Out2 | nominalPressurePrediction_bar | double[3 1] | Nominal desired pressure AFTER identical limits; Scope/log only |
| Out3 | desiredPressureDiagnostic | scalar Bus: ResidualSinglePassDiagnosticsBus | Existing diagnostic Selector/Terminator/logging |

In test4's saved model, **Mux6** combines the outputs of MATLAB Function2 and
feeds MATLAB Function3. Branch that existing Mux6 output to new input2.
It already contains [length(2);d_length(2);dd_length(2)], so no additional
derivative/filter block is needed. Select the wire, branch it, and connect
the branch to input2 after disconnecting the old acceleration-only branch.
Do not disconnect the original Mux6-to-collection-inverse wire.

If assembling the input manually, use a3-input Mux with these vector inputs,
in this exact order:

1. Desired [l2;l3],2x1,m.
2. Desired [dl2;dl3],2x1,m/s.
3. Desired [ddl2;ddl3],2x1,m/s^2.

Use the same reference signals and sample instant as the nominal collection
inverse. Feeding a differently filtered reference will change its force
demand and prevent the nominal pressure comparison from matching.

In Symbols/Edit Data, set desiredMotion size[6 1], pressure outputs[3 1],
and diagnostic size1/type **Bus: ResidualSinglePassDiagnosticsBus**.
The existing diagnostic bus definition is unchanged. Existing bus/options
can be reused. To recreate them from the project root while stopped:

```matlab
addpath(fullfile(pwd,'residual_controller'));
[shadowOptions,ResidualShadowDiagnosticsBus,ResidualShadowContextBus, ...
 shadowPressureOptions,shadowInverseOptions,ResidualShadowInverseDiagnosticsBus] = ...
    setupResidualShadow(pwd);
[singlePassOptions,ResidualSinglePassDiagnosticsBus] = ...
    setupResidualSinglePassShadow(pwd);
singlePassOptions(1)=shadowOptions(1);
```

These setup calls produce default options. Preserve any intentionally changed
existing options rather than resetting them inadvertently. Verify pressure
slots4=0,9=0 for the as-trained collection convention. Slot5 is maximum
pressure(default3bar);slot8 is slew(default6bar/s);slot10 is force/bar.
Common preload/trim and PAC caps(slots1:3,6:7) remain unused here.

## Scope and logging

Plot pressureRequested_bar, new nominalPressurePrediction_bar and new
desiredPressurePrediction_bar together. The new nominal should agree with
the collection nominal when they receive identical desired motion and use
the same nominal physics/map. A separate shadow slew limit can change it
when the actual context pressure differs substantially from that request.

Select and log these diagnostic fields:

- valid,networkEvaluations: false/zero during invalid or cold contexts;
  true/one for a valid evaluation.
- targetForce_N: the desired-state nominal force, comparable to tauDesired_N.
- residualForce_N,correctedForceRequest_N: subtraction constituents.
- nominalPressureBeforeSlew_bar,correctedPressureBeforeSlew_bar: pressures
  after allocation/max-pressure scaling but BEFORE slew. They are not
  unconstrained pressure requests when allocationLimited=true.
- allocationLimited,slewLimited,nominalAllocationLimited,nominalSlewLimited.
- frozenBalanceError_N,frozenBalanceSatisfied: balance with the frozen
  observed residual, not convergence of the pressure-dependent learned inverse.

The slew box is anchored to actual pressure at each context and allows
slew*horizon per channel(.30bar for6bar/s and.05s). There is no separately
accumulating hypothetical pressure trajectory while this block is shadowed.
Consequently a large hypothetical demand may remain limited on every sample
if the arm continues receiving the original nominal commands.

The block consumes the observer's existing context.valid gate. Invalid/cold
contexts return NaNs and valid=false, with no GRU evaluation. No gate is
bypassed to make plots continuous. Both blocks use the same simulation tick;
context featureDt/sampleTime/history remain governed by the forward observer's
physical observation timestamps and missing/repeated sample policy. Desired
reference derivatives use their existing simulation-time trajectory contract.
No recurrent update occurs from a clock tick, pressure candidate or repeat
call to this pressure function.

Save a fresh operator recording with the actual options, for example:

```matlab
save('circlesweep_desiredPressure_r01.mat','out','shadowOptions', ...
     'shadowPressureOptions','singlePassOptions','-v7.3');
```

## Offline evidence

See desired_pressure_validation_20261109_v1/README.md and summary.json.
Both network widths passed desired-target/sign/limits/invalid/immutability
checks. All5433 saved test4 rows were replayed,5118valid and315invalid.
Nominal pressure matches recorded collection request4.45e-16bar; actual
residual matches prior force1.07e-14N. First70s corrected peaks are
[.523,.468,.484]bar, versus recordednominal[.242,.240,.242]bar.
This is saved-data evidence, not a new hardware tracking test or actual
user-model compilation/timing qualification. No retraining was performed.
