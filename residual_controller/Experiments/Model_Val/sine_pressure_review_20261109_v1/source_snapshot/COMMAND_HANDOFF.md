# Updated residualCommandValid: pressure handoff instead of Boolean selector

Operator sine test1 integration finding: the old Switch remained and treated
the3x1 handoff pressure as a Boolean control, injecting435NaN command samples.
Do not retain that Switch in the command path. See
Experiments/Model_Val/sine_pressure_review_20261109_v1/README.md for exact
recorded signal reconstruction. Output1 must drive Unit Delay directly.

2026-11-09. Paste-ready **residualCommandValid.m** replaces the original tiny
validity helper. It now outputs a selected/smoothed3x1 pressure command, not a
Boolean Switch control. Existing observer/GRU/readiness rules are unchanged.
No SLX is edited automatically. Install only after the current operator run
has stopped; do not change the ongoing sinewave experiment.

## Manual ports and wiring

Paste the entire file into MATLAB Function9. Set **Discrete0.05s** explicitly.
It now has7inputs and4outputs. Use fixed-size double3x1 pressure ports and
scalar Boolean flags. The time input must be the same simulation clock used
by the observer's t_sim_s port, not the host pollStart clock.

| Input | Symbol | Connection |
|---|---|---|
|1|pResidual_bar,3x1|Desired-pressure block output1, as before|
|2|residualValid,Boolean|Bus Selector **valid** from desiredPressureDiagnostic; replace the old A5 base-valid wire|
|3|pNominal_bar,3x1|Existing protected original command, currently the Switch's lower input|
|4|pApplied_bar,3x1|Actual NI-command feedback A8/pressure_bar, after physical Kill/saturation and coherent sampling; branch the observer pressure input6 wire|
|5|t_sim_s,scalar|Same sampled simulation clock as observer input8|
|6|enable,Boolean|Existing pressure-enable/Kill-permit signal,1when commands permitted; not NN readiness|
|7|reset,Boolean|Constantfalse during a run; reset only while pressure is disabled|

If the experiment uses a different nominal command source, such as pressure
sine excitation, input3 must be that experiment's existing protected fallback,
not an unrelated inactive circle reference branch.

**Bypass the old Switch.** Connect output1 directly to the **existing Unit
Delay input**. Retain delayTs0.05,IC[0;0;0]. Keep all valve/physicalKill/stop/
pressure-enable wiring downstream. Do not remove this delay or put the new
persistent block in an instantaneous feedback loop.

```
corrected candidate, valid, original fallback,
actual applied-command feedback, clock, enable, reset
                    |
          residualCommandValid
                    | output1: pressure3x1
              existing Unit Delay
                    |
          existing protected valve subsystem
```

Input4 must be actual command feedback, **not** inverseContext.pressure_bar
(which is unset while context is invalid), and not this block's instantaneous
output. Actual commands have already passed through the Unit Delay and
downstream Kill/saturation. This anchors the limiter to the command actually
being applied and avoids retaining a stale high-pressure command after Kill.

| Output | Symbol/type | Log meaning |
|---|---|---|
|1|pressureCommand_bar,double3x1|Selected/rate-limited next command, before external delay|
|2|mode,uint8scalar|State code below|
|3|invalidAge_s,doublescalar|Time since loss of a previously used residual candidate;0during startup/ready|
|4|slewLimited,Booleanscalar|Shared limiter is modifying this tick's target|

No new bus or setup function is needed. Connect outputs2:4 to logging/Scopes
or Terminators. A Boolean output formerly connected to the Switch has changed
to a pressure vector; leaving that old Switch connection is incorrect.

## Defaults editable near the top of the source

| Setting | Default | Effect |
|---|---|---|
|Ts|.05s|Must match actual block and command-delay periods|
|grace_s|.10s|Brief hold of actual command when residual readiness is lost|
|fallbackBlend_s|.50s|Blend from held command to current original fallback|
|reentryBlend_s|.50s|Blend from actual command to newly valid residual candidate|
|slew_barPerS|6bar/s|Shared cap for every ordinary source transition|
|maximum_bar|3bar|Same collection pressure rail|

At.05s, the maximum ordinary per-channel step is.30bar relative to actual
feedback. This applies to fallback as well as residual commands. Grace is
not extended over the full1.6s observer warm-up. If the observer remains
unready, fade toward nominal; do not retain an old rotating-trajectory
direction indefinitely. Returning residual candidates blend in rather than
switching sources abruptly. Fades track the current fallback/candidate values.

Modes:0disabled/reset,1nominal,2residualre-entry,3residual,4invalidgracehold,
5fallbackfade,6inputfault,7timingrestart.

Repeated same-clock evaluations return the cached command, not a second
state/timer advance. A backwards clock or a call interval other than.05s
restarts handoff at actual pressure and reports mode7; there is no large-dt
catch-up step. Set the MATLAB Function explicitly discrete rather than
relying on inherited timing from its continuous nominal input.

The block uses simulation-time hold/fade durations at the existing paced
execution period. It is not an independent physical-time watchdog: a stalled
MATLAB process cannot execute its timeout logic. Existing hardware/pressure
stop protections remain necessary and authoritative.

Disable/reset outputs zero and clears handoff state. Nonfinite time/applied
feedback outputs zero/mode6; these exception paths bypass ordinary smoothing.
Reset only with pressure disabled. A nonfinite fallback is never blended as
NaN: the fallback target becomes zero under the limiter and mode6 is logged.
Neither valid=true nor a finite candidate can override enable=false.
Physical Kill/stop must act immediately downstream and must not wait for
the smoothing ramp or external Unit Delay.

## What this change does not establish

GRU invalidation/reset and observation gates remain untouched. A pressure
handoff can reduce switching jerks; it does not restore missing NDI frames
or guarantee phi remains constant during extended warm-up. This is not a
complete hardware watchdog or a new stability proof. Source implementation
and offline command-path tests do not qualify actual SLX compilation/timing.

Verification: command_handoff_validation_20261109_v1. Synthetic startup,
dropout, re-entry, disable/reset, duplicate ticks, bad inputs and timing
restart cases pass. The saved test5 candidate/validity/nominal signals drove
a simulated command-delay path with no arm/NN-state recomputation:7501ticks,
maxordinarynewstep.300bar versus recorded1.514bar,0newsteps>.3bar versus30
recorded. This counterfactual verifies the command algorithm only; its future
measured phi, candidate forces and validity would change on the real arm.
Keep the current sinewave recording as its own baseline; use a new filename
for a later manually run handoff-enabled experiment.
