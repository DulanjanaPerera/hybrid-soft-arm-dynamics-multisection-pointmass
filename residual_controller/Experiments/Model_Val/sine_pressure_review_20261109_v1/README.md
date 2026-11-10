# Sine desired-pressure test1: old Switch bypasses the new pressure handoff

2026-11-09. Offline sinewave_FWD_INV_test_1_fwdPressure.mat only. Model1.59,
ReachedStopTime212s/4241observer ticks, only out saved. Operator execution
contains the new four-output pressure handoff. The separately saved SLX still
contains the earlier Boolean-helper revision; no model load/edit/run was used
to infer the executed command path. Numeric identities below establish it.

## Primary finding: pressure output is still used as Switch control

The logged signal named **useResidual** is now a3x1 pressure vector from the
new helper, not aBoolean. All4241handoff output vectors are finite. A complete
sequential replay of the new source using corrected candidate, diagnostic
validity, collection-inverse fallback and actual NI feedback matches:

-pressure outputs1.49e-14bar,
-all mode codes exactly,
-all invalid ages exactly,
-all limiter flags exactly.

However, selected_pressure does NOT equal the preceding handoff output;
maximum discrepancy is~3bar. It instead matches this reconstruction EXACTLY:

```
for each pressure channel:
    if handoff_output(channel) > 0:
        switch_output(channel) = raw_residual_candidate(channel)
    else:
        switch_output(channel) = raw_sine_command(channel)
selected_pressure(k) = switch_output(k-1)
```

Finite values match0bar maximum error, and the NaN pattern matches exactly.
This is the retained old vector-controlled Switch, with raw sine command on
its lower data input. The smoothing block is computing pressure correctly,
but its output is being interpreted as per-channel Boolean selection rather
than being used as the actuator command. The old signal name concealed the
interface change. The requested removal/bypass of that Switch was not realized
in the executed signal path.

## Fault sequence and oscillation

| Simulation time | Logged evidence |
|---|---|
|145.85s|Observer ready; finite residual candidate|
|145.90s|NDI tip missing; candidate becomesNaN; handoff holds actual finite pressure in mode4|
|145.95s|OldSwitch selects NaN candidate because held pressure>0; delayed NI-command log containsNaNs; handoff reports mode6|
|146.00s|Zero handoff control selects finite raw sine fallback; observer starts a new segment|
|146.05s|Finite nonzero handoff pressure selects still-cold NaN candidate again|
|Following ticks|Alternation repeats, blocking normal observer re-warm|

435observer-grid samples have at least one nonfinite NI command. The handoff
has435input-fault mode6samples between145.95and198.60s, but **its own pressure
outputs remain finite throughout**. The NaNs are injected downstream by the
old Switch. The recorded NI-command values are not measured chamber pressure;
the actual analog-output hardware's handling of NaNs is not established here.

Only5orientation-invalid samples are actual tipMISSING; base alwaysvalid.
434observer invalid samples occur with both sensors valid, consistent with
its nonfinite-pressure rejection.2physical sample intervals exceed.1s.
Command-age logs remain0. Thus these large observer invalidations cannot be
attributed to hundreds of physical sensor dropouts. The wiring creates a
separate invalid-command feedback cycle after an initial real interruption.

Modes:537nominal,27re-entry,3240residual,2gracehold,0fallbackfade,435inputfault.
No handoff slew-limit flags occur: normal raw candidate changes are already
within its shared step allowance, and the new handoff output is not driving
the actuator path anyway. A finite-adjacent NI-command step maximum.3bar does
not make this a valid/safe command stream: the NaN samples must be counted.

## Required manual wiring correction before another actuation trial

```
MATLAB Function9 output1 (pressureCommand_bar / old wire useResidual)
      -> existing Unit Delay input
      -> existing protected Pneumatic Valves
```

Bypass the old Switch completely. Do not connect pressure output to its middle
control port, and do not use that Switch's output as the Unit Delay input.
KeepUnitDelayTs.05s/IC[0;0;0], physicalKill/stoppriority downstream.

Handoff input3 should use the intended protected original command, currently
the raw sine generator pressureRequested_bar from /sineWave. Executed replay
shows it instead used /MATLAB Function3's reconstructed nominal inverse,
differing up to.09135bar. That small fallback-source difference is secondary
to the retained Switch, but should be made consistent with the experiment.
Use input2desiredPressureDiagnostic.valid,input4actualA8command feedback.
Rename the old useResidual wire to pressureCommand_bar to clarify its type.
Mode/invalid-age/slew outputs are diagnostics, not Boolean source selectors.

Keep the bad run as evidence. A new manually installed/build/run recording
is required before assessing this handoff on hardware. No autonomous rerun,
SLX edit, pressure action or runtime-code change was performed in this audit.

## Other limits, not solved by this correction

868/3267valid inverse rows hit3bar peak-pressure scaling,3hit candidate slew.
Desired force matches recorded tauDesired1.43e-14N. P1/P2 active length RMSE
is[2.153,5.084]/[5.335,1.938]mm; corrupted P3command phase gives
[9.440,18.379]mm. These are descriptive, not evidence against a correctly
connected handoff or a matched controller performance comparison.

A supplementary frozen-history pressure-map test at654valid contexts found
local spectral radius median/P95/max.188/.237/.257,0above1. This is only the
instantaneous single-pass allocation map with q/v/history/target frozen,
NOT the full recurrent/plant/control Jacobian. It neither proves stability
nor rules out remaining model/constraint-driven oscillations. The proven
wiring fault must be corrected before interpreting the remaining oscillation.
Candidate-pressure model-balance norm median/P95/max.131/16.847/21.272N
reflects pressure dependence/limits; not measured physical force errors.

## Retained evidence

summary.json,retained_switch_check.json,handoff_source_replays.csv,
handoff_modes.csv,channel_metrics.csv,trace.csv,analysis.mat;
overview/P1_detail/P3_detail PNG/PDF/FIG,source/raw hashes/source_snapshot.
Reproduce via analyzeSinePressureRun(root,freshFolder). Actual in-model runtime
source/options were not saved with out; source equivalence is established
numerically, not claimed from the stale on-disk SLX.
