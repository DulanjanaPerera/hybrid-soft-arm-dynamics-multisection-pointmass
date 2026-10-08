# Current status reference

Read [project state](../project_docs/PROJECT_STATE.md) and
[open issues](../project_docs/OPEN_ISSUES.md) before using these experiment
instructions. The October 6 recording exposed an invalid-input early return
that bypasses the feedback inverse's normal pressure slew limiter. The nominal
checks described below do not cover that fault transition. The saved
FBL_commonP.slx now contains both supplied feedback scripts; wiring, current
Kill state and equivalence to the recorded run remain unverified.
# Common-pressure experiment

## Paste-ready offline circle report

`analyzeCommonPressureCircle.m` reads a saved `out.logsout` MAT file and prints
a short report to the MATLAB Command Window. It detects the logged phi target,
circle count and timing from the desired signals, so the filename and any
intended settings cannot silently substitute for what was actually logged.
It summarizes the five seconds before circle motion (which may be a phi ramp),
measured versus desired l2/l3 and phi,
logged Kp/Kd gains (including whether they varied), force components,
validity, pressure limiting, common/differential pressure, command versus
post-Kill NI-path agreement, and each circle's tracking. It does not open a model,
run hardware, draw plots or write files.

From the repository root in MATLAB:

```matlab
addpath('common_pressure_testing')
analyzeCommonPressureCircle('FB_circle_phi_0.2rad.mat');
% Or use an absolute MAT path. With no argument, it reads the newest
% FB_circle*.mat beside the analyzer.
```

Copy the printed `COMMON-PRESSURE CIRCLE REPORT` when sharing new results.
NI-path equality compares logged software commands; actual chamber pressure
is not measured. Measured lengths are inferred from NDI orientation, not
independent muscle-strain sensors. The existing
For plots of the current common-pressure logs, use
`plotCommonPressureTracking.m`:

```matlab
addpath('common_pressure_testing')
plotCommonPressureTracking(out)             % current SimulationOutput
% plotCommonPressureTracking(out.logsout) % dataset directly
% plotCommonPressureTracking('FB_circle_phi_1rad_trial6_good_10s_0.5bar.mat')
```

It opens six figures: desired/measured theta and phi, wrapped theta and phi
errors, desired/measured l1/l2/l3, the two inverse-dynamics force channels,
three commanded pressures, and desired/measured tip-sensor XYZ. The function
uses `meas_length` and NDI `position` from the current logged model, hides
invalid orientation/position samples, and displays desired and measured theta
wrapped to [-pi,pi] with gaps at wrap crossings. It does not save files,
load a model, or access hardware.
Force plots show `M*v`, `(C+D)*dq`, `G+Keff*q`, and requested total; pressure
plots are software commands, not independently measured chamber pressure.

## Next feedback test: paste-ready V2 fault handling

`commonPressureFeedbackInverseBlockV2.m` is the revised inverse MATLAB Function
block code. Paste its entire contents over the inverse block code in the user's
test model; the four inputs and eight outputs retain the V1 order. The saved
`FBL_commonP.slx` was **not** edited by Codex. Keep
`commonPressureFeedbackReferenceBlock.m` and its 0.01 rad U/V references as-is.
The gain blocks are outside the inverse function: for one cautious comparison,
set both Kp entries from 1100 to 1400 and both Kd entries from 100 to 120.
These values are a test proposal, not experimentally validated tuning.

V2 holds the last pressure command across a brief invalid model/NDI sample.
If valid data has not returned within 0.10 s of the last valid sample, it
latches a fault and slews all three commands toward zero at 0.8 bar/s; restart
the simulation to clear the latch. It still commands zero immediately if
Clock time itself is invalid, since a time-based slew is then unavailable.
The separate hardware Kill path remains necessary. Log `inverseModelValid`,
NDI validity, controller `pressure_bar`, and post-Kill regulator commands if
available. Monitor the previously observed local buckle. The source does not
add integral or adaptive trim correction.

Offline MATLAB fault injection with
`testCommonPressureFeedbackInverseBlockV2.m` passed single-frame hold,
recovery, sustained-fault slew/zero, and no automatic re-pressurization after
latch. This did not compile or run the Simulink model or access hardware.
The subsequent `feedback_testing_results_V2.mat` retained Kp=[1100;1100]
and Kd=[100;100] and omitted all eight formerly logged inverse-block outputs,
including `pressureCommand_bar` and `inverseModelValid`. Restore those signals
and post-Kill regulator commands before using a new recording to judge the
physical invalid-frame pressure response. See the current validation results
in `project_docs/VALIDATION_RESULTS.md`.

## First bounded feedback test (paste-ready blocks; no model edited here)

Use a **copy** of `feedback_linearization_trackingWorks.slx`. Paste
`commonPressureFeedbackInverseBlock.m` into its root `MATLAB Function1`
inverse-model block. Its first three inputs and eight outputs preserve the
existing order; connect the new fourth input `t` to `Reference clock`.
The block retains the saved dynamic model and uses the manual trial-4 trim
as a 1.6 bar mean, all-active pressure allocation. At 50 s it anchors the
nominal static model force at the measured pose; it then applies incremental
feedback corrections, clamped to +/-0.04 bar in both `s` and `a`. It limits
each pressure output to 0.8 bar/s internally. The old minimum-pressure
allocator and its per-active-channel 0.8 bar offset are replaced; do not add
that offset again. `tauDesired_N` remains the absolute model force, while
`tauAvailable_N` reports the ideal force of the actual rate-limited pressure
command, so they need not be numerically equal at the empirical trim.

For the first test, paste `commonPressureFeedbackReferenceBlock.m` into root
`MATLAB Function7` in the same model copy. Keep `Reference clock` connected
to its time input. Connect its `theta` output to the existing theta branches
of the desired-angle Mux and NDI comparison input, and its new `phi` output
to their phi branches in place of the old `Switch`/Ramp signal. This removes
the original continuously rotating circle and 1.5 rad bend reference from
the test. Set Stop Time to **240 s**. The protocol is zero pressure for 10 s;
trial-4 trim [1.6656,1.6492,1.4852] bar until 50 s; feedback toward
geometric zero from 50 s; a smooth 0.01 rad U-direction reference step from
90-130 s and V-direction step from 160-200 s; and a return to zero pressure
from 220 s. Visually monitor the known local buckle and retain the physical
Kill control. Log desired/measured lengths, U/V or theta/phi, tracking error,
`pressureCommand_bar`, post-Kill regulator commands, `pressureLimited`,
`inverseModelValid`, and NDI validity/timestamps.

The block pair was evaluated as pure MATLAB code offline on a 0.05 s grid.
The mock run satisfied finite 0-3 bar outputs, the initial/final zero
windows, +/-0.8 bar/s pressure slew, and the expected reference phases.
A representative measured V=-0.027 rad produced pressure changes in the
direction needed to increase V. These checks do not compile the pasted
Simulink model or validate NI/NDI hardware behavior. The user's working
analysis assumption is command pressure equals chamber pressure.

## Local calibration around manual trial 4

Paste the entire `commonPressureLocalCalibrationBlock.m` function into one
Simulink MATLAB Function block. Connect a Clock to its `t` input and route
the three pressure outputs through the existing Mux, downstream Rate Limiter,
and valve Kill path. Log `phaseID` along with the post-Rate-Limiter pressure,
post-Kill regulator-path commands, NDI lengths/orientation/XYZ, and validity.
Set Stop Time to **395 s**. Positive phase IDs label holds; negative IDs label
five-second transitions. The block stays at zero pressure after 395 s.

It starts at zero, holds the manual trial-4 trim
`[1.6656,1.6492,1.4852]` bar, probes `s_bar` by +/-0.04 bar one direction
at a time, then probes `a_bar` by +/-0.04 bar, always returning to the trim
between probes, and finishes at zero. Active plateau commands span
1.4452-1.6923 bar. The transition is a five-second smoothstep, with a
maximum generated slew of 0.50 bar/s; retain the existing downstream Rate
Limiter. All pressure-to-motion interpretation uses the project's working
assumption that commanded pressure represents chamber pressure. The schedule
is open-loop and the manual visual stop condition for any growing/sudden
local buckle remains important because NDI tip logging does not measure
the distributed arm shape.

`common_pressure_testing.slx` is the user-created forward/inverse comparison
model. `commonPressureSequenceBlock.m` is the 300 s open-loop P1/P2/P3
command schedule. Neither analysis script nor MATLAB command below loads or
runs the model or accesses NI/NDI hardware.

From the repository root in MATLAB:

```matlab
addpath('common_pressure_testing')
results = analyzeCommonPressureTest;
% Or: results = analyzeCommonPressureTest('C:\full\path\to\another_run.mat');
```

The default input is `common_pressure_testing/test_300s.mat`, which must
contain `out.logsout`. The analysis writes to
`common_pressure_testing/results/test_300s/`:

- `phase_summary.csv`: the final five simulation seconds of each logged
  `phaseID`, including commanded common and differential pressure, validity,
  lengths, signed bending, forward-model prediction, NDI XYZ, and within-tail
  range/slope. A final-five-second average is not proof of equilibrium.
- `condition_summary.csv`: the neutral, positive and negative bend response
  at each common level; signed U/V gains and return-to-neutral shifts.
- `report.txt`: a compact numerical summary and interpretation limits.
- `pressure_trajectories.png`, `lengths.png`, `bending.png`,
  `condition_comparison.png`, `tip_position.png`: offline figures.
- `analysis_results.mat`: aligned data and both tables for further work.

The script uses `meas_length` from NDI and selects the forward model's
`comp_lengths` by block path, since the log has duplicate `comp_lengths`
names. It aligns model and command signals to the 0.05 s NDI time grid and
excludes invalid orientation, sensor-pair, and configuration samples.
`U=(l2+l3)/0.013` is signed bending toward the nominal theta=0 direction;
`V=(l3-l2)/(sqrt(3)*0.013)` measures cross-plane bending. These preserve
direction where the nonnegative bend magnitude `phi` cannot. Tip XYZ changes
can also result from bending at fixed backbone length.

`des_pressure`, `pressureUsed_bar`, and `P1/P2/P3 Regulator` are **logged
software signals**, not measured chamber pressure. The latter three describe
the regulator command path after the model's gating/saturation. Matching them
to `des_pressure` does not establish actual pressure delivery. The apparent
scalar stiffness estimate in the condition table assumes equal nominal
pressure-to-force gains, static equilibrium, and
`K=k*[2 1;1 2]`; unequal PMAs, valve dynamics, hysteresis, drift, and
cross-plane motion can all change that estimate. Treat it as a diagnostic
command-to-motion conversion, not an identified material parameter.
