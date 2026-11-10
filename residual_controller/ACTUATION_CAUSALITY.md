# Breaking the pressure/observer algebraic loop

Updated code is now available in residualCommandValid.m; follow
[COMMAND_HANDOFF.md](COMMAND_HANDOFF.md). It outputs the pressure command
directly and bypasses the old Switch. The Boolean-only helper/sketch below
is historical and does not provide the new smooth handoff.

Operator test5 update: the Switch/delay topology broke the loop, but the
simple Switch makes abrupt fallback pressure steps. See
Experiments/Model_Val/desired_pressure_test5_review_20261109_v2/README.md.
It is not a smooth handoff supervisor: the original fallback bypasses the
corrected branch's slew cap. Next work should add bounded grace/transition
logic and a shared final slew limiter while retaining Kill/stop priority;
do not interpret the wiring sketch below as complete closed-loop qualification.

2026-11-09. Operator screenshot from residual_FWD_INV_validation_v2 identifies
MATLAB Function4(persistent observer), MATLAB Function8 input1(context), and
valve switches/saturations/ZOH/RateTransition in an algebraic loop.
Read-only ZIP/XML inspection confirms the newly connected path:

```
desiredPressurePrediction
 -> Pneumatic Valves
 -> sampled actual command + commandTime
 -> observer inputs6/7
 -> inverseContext
 -> desiredPressurePrediction
```

The observer updates persistent state while producing its outputs. This path
cannot be solved as a simultaneous algebraic equation. A ZOH or Rate Transition
does not guarantee a one-tick delay in this loop. Do not suppress the diagnostic
or change observer history updates to work around it.

## Manual causal command path

Insert the following path manually, without changing the physical Kill or
existing valve saturation/output wiring:

```
desiredPressurePrediction_bar ----+
                                 | true
                            command Switch -> Unit Delay -> Pneumatic Valves
                                 | false
existing nominal/startup command-+
```

Switch condition: desiredPressureDiagnostic.valid AND all three candidate
pressures are finite. A small combinational MATLAB Function can compute the
condition if needed:

```matlab
function useResidual = residualCommandValid(p, valid)
%#codegen
useResidual = valid && all(isfinite(p));
end
```

The false branch must remain the existing protected nominal/startup command,
such as the original pressureRequested_bar branch through its existing
enable/stop logic. Do not wire the new NaN-producing prediction directly to
the Unit Delay. Keep original stop/pressure-enable and physical Kill priority;
valid is not a substitute for those controls or a complete dropout supervisor.
Do not run hardware until the startup, invalid-data and stop behavior is
explicitly resolved and checked by the operator.

Unit Delay parameters:

- Sample time:0.05s, matching the desired-pressure block/observer.
- Initial condition:[0;0;0]bar, a3x1 column.
- Input/output:double3x1.

With this placement, at tickk the actual command from the delay is the selected
candidate computed at tickk-1. The observer sees that actual command at tickk,
advances its existing history once for a new valid observation, and computes
the next candidate. There is a deliberate50ms command delay. This is causal
plumbing, not a proof of closed-loop stability or tracking performance.
The original desired feedforward reference is also effectively delayed50ms;
any later reference preview/control redesign must account for this explicitly.

Keep the existing feedback sampling packet intact:

```
final valve commands(3) + command-source simulation timestamp(1)
 -> existing coherent sampling/rate transition -> observer pressure/time
```

Do not delay only pressure or only timestamp. Do not insert a feedback-only
delay casually: the forward observer forecasts under the pressure it receives;
a delayed feedback value could differ from the new command actually applied
during that forecast. The command-path delay avoids that mismatch while the
new command is held over its control interval. Existing physical Kill changes
and within-interval pressure changes still require their usual diagnostics.

The desired-pressure block and forward observer code need no edits solely to
break this loop. Update Diagram manually after wiring; no model was loaded,
edited, compiled or run during this investigation. Hardware actuation remains
operator work; the offline/native tests were component tests, not qualification
of this newly closed loop.
