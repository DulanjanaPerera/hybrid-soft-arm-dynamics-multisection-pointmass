# Verification and experimental results

Results below summarize existing recorded evidence; migration did not rerun
MATLAB tests or operate hardware. Distinguish calculation checks, model smoke
tests and measured behavior. Full settings/provenance remain in linked local
READMEs, original reports and dated log entries.

## NDI / estimator evidence — 2026-10-01

Straight calibration: 100/100 accepted pairs, 2.30 mm XYZ reference error;
separate preflight 98/100 pairs, median interval 0.0502 s. These checks do not
establish global orientation/strain accuracy after mounting changes.

Six 30 s sine trials, 601 samples each, ten invalid orientation samples among
3606. Filtered difference versus four-state CV Kalman used 0.05 mm measurement
SD and 0.003 m/s² per-step acceleration SD. At tau=0.12 s the mean prediction
proxy was about 0.106 versus 0.110 mm; no accuracy winner established.

| Filter tau (s) | Low-command velocity RMS (mm/s) | Prediction RMS (mm) | Velocity-change RMS (mm/s) |
| ---: | ---: | ---: | ---: |
| 0.12 | 0.0633 | 0.1061 | 0.6135 |
| 0.20 | 0.0534 | 0.1027 | 0.4064 |
| 0.30 | 0.0464 | 0.1011 | 0.2868 |
| 0.50 | 0.0381 | 0.1009 | 0.1827 |

Status: offline comparison completed; true velocity and controller phase lag
unvalidated. [Method and saved outputs](../estimator_comparison/README.md).
Source recording hashes and per-trial baseline metrics remain in the October 1 log.

## Early feedback implementation — 2026-10-01

Synthetic checks passed missing-frame hold, sine derivatives, error sign,
nonnegative/bounded pressures, force reconstruction, saturation and invalid-frame
zero command. Hardware-free copy with NI/NDI removed compiled and ran 0.5 s,
pressure 0–0.8002 bar. Status: narrow offline smoke pass; precedes later routing
and dashboard revisions and does not certify the current model or plant.

## Inverse and stiffness evidence — backfilled 2026-10-05

- Inverse recording `simulationresulta_20261001_1936.mat`: 2659 samples,
  167.836 s physical time, 99.85% valid, 7.48% limited, three qualifying holds,
  two stationary tails, one target maintained under 1 mm / two-second rule.
  Status: analysis completed, partial tracking success. Zero-rate holds cannot
  validate dynamic terms. [Report](../inverDynamic_test/results/report.txt),
  [method](../inverDynamic_test/README.md).
- Stiffness sweep `Comparision_results.mat`: grid 1350 down to 200 N/m;
  best tested 650, reduced-coordinate RMS 3.55231 mm over 600 valid samples;
  nominal baseline replay discrepancy <=4.61e-7 m. Status: replay consistency
  pass and conditional same-recording fit, no independent physical stiffness
  identification. [Report](../FWD_INV_COmparision/stiffness_sweep_results/report.txt).
- P1 static fit: means about 4132 N/m at phi=0.09425 to 953.333 at 2.05594 rad;
  low-angle SDs 1843 and 2597 N/m. Zero deadzone, pooled loading/unloading,
  command pressure assumed. User reported poor unloading agreement using phi
  alone. Status: empirical interpolation available, general law unvalidated.
  [Method and outputs](../planar_stiffness_test/README.md).

## Existing tracking / plotting — 2026-10-04–05

`testing_20261004_1311.mat`, 95–125 s: desired/measured phi=1.768/1.049 rad,
mean reduced errors [+9.37,-4.01] mm, requested force [16.01,0.39] N,
command [0,2.006,0.829] bar, limiting false. The inspected inverse used K=250,
not later K=650. Status: substantial steady error; modeled force consistency
does not demonstrate delivered force accuracy.

`plotCircleTracking.m` produced seven figure groups offline for three saved circle
recordings and an out-object input. Status: plotting execution pass, no tracking
performance claim; orientation-derived versus measured XYZ may differ.

## Common-pressure 300 s — 2026-10-05

6001/6001 valid; physical interval 0.0508 s, span 304.58 s. Final-five-second
tails at common means 1.2/1.6/2.0 bar gave neutral U=0.1065/0.1702/0.1958 and
V=0.0717/0.1444/0.1523 rad. Measured differential U gains
0.3096/0.4649/0.3897 rad/bar versus modeled 0.9641. Largest tail U range
0.00477 rad. Status: ideal equal-pressure neutrality disagrees with measured
pose; software commands matched regulator-path logs, chamber sensing absent.
Apparent K=2199/1464/1747 N/m is a conditional diagnostic.
[Report](../common_pressure_testing/results/test_300s/report.txt).

## Neutral-trim 387 s / validation 199 s — 2026-10-05

387 s: 7741/7741 valid; reference U,V=(0.03051,0.01501) rad.
Local Jacobian condition numbers 1.57/1.39/1.47 support distinct command
directions at the three means. Proposed triplets [1.328,1.292,0.981],
[1.766,1.766,1.268], [2.236,2.172,1.592] bar extrapolate beyond tested probes.
Status: candidate estimation only. [Report](../common_pressure_testing/results/trim_test_387s/report.txt).

199 s: 3980/3981 valid. The mean-1.6 candidate overshot reference with U/V
residuals (-0.10136,-0.15204) and (-0.09697,-0.14795) rad on two approaches.
Intermediate [1.7,1.7,1.4] residual (+0.01795,-0.01143), XYZ norm 6.57 mm
versus candidate ~26 mm. Status: extrapolated candidate failed neutral target;
repeat results do not prove equilibrium. [Report](../common_pressure_testing/results/trim_test_199s/report.txt).

## Manual trial-4 cycling — 2026-10-05

635.25 s, 12706/12706 valid, 17 on/17 zero episodes, four triplets.
Best of four was [1.6656,1.6492,1.4852] bar. Five later repeats:
U=0.02293 (SD 0.00103), V=0.00690 (SD 0.00247), magnitude=0.02404
(SD 0.00107) rad. Relative bend change 0.0047–0.0080 rad and XYZ shift
3.8–4.3 mm. Status: within-run local candidate, not exact straightness;
first approaches differ and short holds may still drift.
Source: `common_pressure_testing/manual_test_results.mat`; dated log entry.

## Local calibration 395 s — 2026-10-06

7900/7901 valid. +/-0.04 bar s/a probes around trial 4 gave approximate
`J=[[0.397,0.069],[-0.079,0.490]]` rad/bar after neighboring-trim subtraction.
Five identical-command trim holds spanned U/V=0.0096/0.0138 rad and differed
from manual-run means; tail slopes up to ~0.00035 rad/s. a+ achieved
(U,V)=(0.0067,-0.0007), but at 237.30 s an unexplained ~0.0063 rad
orientation-vector step occurred without a corresponding large XYZ/command
change. Status: useful local directions, equilibrium/repeatability incomplete.
Source: `common_pressure_testing/local_calibration_395s.mat`; dated log entry.

## Bounded feedback numerical checks / first recording — 2026-10-06

Pure-source checks at 0.05 s passed nominal bounds, zero windows, normal
0.8 bar/s slew and reference phases; extreme-acceleration case reached
[1.6656,1.6892,1.4452] bar, mean 1.6, limited=true. Status: nominal numerical
pass only; does not cover the subsequently observed invalid-sample jump.

Physical recording: 240 s, 4800/4801 valid; Kp=[1100;1100], Kd=[100;100].

| Window (s) | Desired U,V (rad) | Measured U,V (rad) | Mean length-error norm (mm) |
| --- | --- | --- | ---: |
| 115–125 | 0.0100,0 | 0.0054,-0.0065 | 0.112 |
| 185–195 | 0,0.0100 | 0.0016,-0.0032 | 0.211 |

Relative U response +0.0016, V response +0.0036 rad; small final slopes
~<=0.00006 rad/s support persistent residuals. Mean pressure stayed 1.6;
ds range -0.00749..+0.00288, da +0.00390..+0.00973 bar; no limiting in
3300 feedback samples. Initial trim drift confounds claimed zero-target improvement.
At 14.40 s invalid NDI produced an immediate zero command, ~33.3 bar/s;
recovery 0.8 bar/s. Status: intended response direction, incomplete tracking;
fault-path slew requirement violated. No post-Kill logs in this MAT.
Source: `common_pressure_testing/feedback_testing_results.mat`; latest log entry.

V2 source-only fault injection on October 6 passed: a single 0.05 s invalid
frame held the prior pressure, valid input resumed without a pressure step,
a sustained invalid interval slewed toward zero at no more than 0.8 bar/s,
and a latched fault did not automatically re-pressurize. Executed with MATLAB
R2025a using `testCommonPressureFeedbackInverseBlockV2.m`. Status: offline
source check only; no Simulink compile or NI/NDI run. The saved model still
contained V1 at the time of that source check; later saved-model inspection
finds V2 in the main inverse chart.

## V2 physical recording — 2026-10-06

`common_pressure_testing/feedback_testing_results_V2.mat` spans 240 s with
4800/4801 valid NDI samples; the one invalid sample was at 12.65 s. Both gain
logs stayed at Kp=[1100;1100], Kd=[100;100], so the proposed gain increase
was not part of this test. The desired U/V steps remained 0.01 rad.

| Hold (s) | V1 measured U,V (rad) | V2 measured U,V (rad) | V1/V2 length-error norm (mm) |
| --- | --- | --- | ---: |
| U, 115–125 | 0.00543,-0.00649 | 0.00415,-0.00760 | 0.112/0.133 |
| V, 185–195 | 0.00156,-0.00324 | -0.00040,-0.00555 | 0.211/0.248 |

The V2 zero-target 80–90 s pose was U,V=(0.00168,-0.00555), versus
(0.00379,-0.00530) in V1. Relative to that local baseline, the U step moved
U by +0.00247 rad in V2 versus +0.00164 in V1. Relative to the 150–160 s
zero-target hold, the V step moved V by +0.00221 in V2 versus +0.00360 in V1.
Both runs retained offset; different trim poses/history prevent attributing
the larger V2 absolute errors to V2 fault handling.

V2 `logsout` has 48 elements versus 56 in V1. It omits
`pressureCommand_bar`, `pressureLimited`, `inverseModelValid`,
`tauAvailable_N`, `tauDesired_N`, and the three inverse force outputs.
No post-Kill pressure command is logged. Thus the MAT
cannot establish whether pressure held through the invalid 12.65 s frame or
whether the local correction bounds were active. Read-only inspection of
`FBL_commonP.slx`, saved before the MAT file, found V2 in the main
`MATLAB Function1` chart and V1 in `MATLAB Function8`. This supports the
intended source version but is not a measured pressure-path check.

## Higher-gain V2 recording near straight — 2026-10-07

`common_pressure_testing/feedback_testing_results_V2_newGains.mat` spans
240 s at 0.05 s logging intervals. It records Kp=[1400;1400] and
Kd=[120;120], constant, and the same 0.01 rad U/V references. NDI and inverse
model were valid for 4800/4801 samples. At the single invalid sample, 22.00 s,
controller and logged post-Kill NI commands held exactly
[1.6656,1.64918,1.48522] bar from the prior sample and resumed at the same
value on the next valid sample. All three logged NI-command channels matched
controller pressure throughout. The Kill signal was logged enabled at its two
event samples. Actual chamber pressure was not measured.

| Hold (s) | Desired U,V (rad) | Measured U,V (rad) | Mean reduced length-error norm (mm) | Mean controller P1,P2,P3 (bar) |
| --- | --- | --- | ---: | --- |
| Zero, 80–90 | 0,0 | 0.00102,0.00108 | 0.020 | 1.6758,1.6459,1.4783 |
| U, 115–125 | 0.0100,0 | 0.00422,0.00449 | 0.090 | 1.6680,1.6502,1.4818 |
| Zero, 150–160 | 0,0 | 0.00143,0.00152 | 0.028 | 1.6759,1.6460,1.4782 |
| V, 185–195 | 0,0.0100 | 0.00157,0.00167 | 0.134 | 1.6759,1.6389,1.4853 |
| Zero, 210–220 | 0,0 | 0.00139,0.00148 | 0.027 | 1.6759,1.6460,1.4782 |

Relative to each preceding zero hold, U advanced +0.00320 rad for the U
reference, while V advanced only +0.00015 rad for the V reference. Absolute
errors are below both earlier recordings, but the V step was not independently
tracked. Feedback pressure maintained a 1.6000 bar common mean. Per-channel
command ranges during 55–220 s were [0.00920,0.01147,0.00857] bar, maximum
slew was [0.00803,0.01223,0.00709] bar/s, and pressureLimited was false for
all 3300 feedback samples. Overall maximum pressure rate, including the
start/end ramp, was 0.8 bar/s. NDI tip XYZ per-sample changes during feedback
had median 0.035 mm, 99th percentile 0.100 mm and maximum 0.239 mm. No tip
jump or command spike is evident in this recorded interval; distributed local
buckling cannot be assessed from tip sensing.

The near-straight interpretation is limited: thetaObservable=false during
3298/3300 feedback samples, yet inverseModelValid=true. Median feedback phi
was 0.00270 rad; its maximum 0.00897 rad was near the calibrated 0.5-degree
(0.00873 rad) observability threshold. During all five steady holds in the
table, thetaObservable was false and theta was held at 0.8166 rad by the NDI
orientation estimator. The logged measured lengths are inferred from phi and
this held theta, so U/V cannot resolve a change in bend direction there.
The stronger gains produced smooth bounded commands near straight in this
single run; they do not establish directional tracking or general stability
at the singular pose. Source: `ndiSupportFnc.orientationFromPose` and
`AuroraReader_s2_Position_only` for the held-angle behavior.

## Common-pressure circle recordings — 2026-10-07

Read-only analysis of `common_pressure_testing/FB_circle_phi_0.2rad.mat` and
`FB_circle_phi_1rad.mat` (both 240 s). The second filename is misleading:
its logged desired phi is **0.1 rad**, not 1 rad. The first run ramped phi
over 5–20 s; the second over 5–10 s. Both held theta=0 through 35 s, then
requested one smooth revolution over 35–155 s and a phi return over
155–175 s. The saved `FBL_commonP_circle.slx` inspected after the runs has
one-input `thetaCircle(t)` with `des_phi=0.1` internally; it does not yet
contain the later proposed `des_phi`/`circle_time_s` input interface. This
saved source should not be assumed to be the exact source for the earlier
0.2 rad recording. The user subsequently chose internal editable phi,
per-lap duration and circle-count variables in the one-input reference;
that revised code is also not installed in the saved model.

| Recording / logged target | Mean measured phi, 30–35 s | Desired l2,l3 (mm) | Measured l2,l3 (mm) | Mean controller P1,P2,P3 (bar) |
| --- | ---: | --- | --- | --- |
| 0.2 rad | 0.1213 rad | 1.30, 1.30 | 1.1625, 0.3410 | 1.53485, 1.72173, 1.54342 |
| Named 1 rad; logged 0.1 rad | 0.0613 rad | 0.65, 0.65 | 0.6463, 0.0808 | 1.60849, 1.68573, 1.50578 |

The corresponding l3 errors are 0.9590 and 0.5692 mm, whereas l2 errors
are only 0.1375 and 0.0037 mm. Across the circle window (35–155 s), mean
measured phi is 0.1156 rad versus 0.2 desired, and 0.0509 rad versus 0.1
desired. Logged theta and tip XYZ show substantial rotation/motion, so the
failure is bend-magnitude/directional tracking rather than no motion.

All three logged post-Kill NI-path pressure-equivalent signals equal
`pressureCommand_bar` at every recorded sample in both files (maximum
difference zero). `inverseModelValid` and NDI `orientationValid` are true
throughout the 10–220 s analysis window; `pressureLimited` is false there.
At 32 s in the 0.2 rad run, logged acceleration correction is about
[1.04,1.72] N and gravity/stiffness force about [1.88,1.30] N. The
controller therefore supplies nonzero correction while a persistent error
remains. The result is consistent with finite PD error under static
model/actuator mismatch, particularly in l3; it does not identify the
specific cause. Post-Kill NI-path agreement verifies the software command
route, not chamber pressure, which was not sensed. No model, controller,
hardware or analysis script was changed for this inspection.

The later `analyzeCommonPressureCircle.m` utility was run offline against
both MAT files and reproduced these hold values, the logged 0.1 rad target
in the second file, and zero controller-to-NI-path command discrepancy. It
prints a per-lap report for future recordings; no hardware/model validation
is implied by that offline check.
The expanded report also reads sparse Kp/Kd change-event samples correctly:
both existing runs report Kp=[1400,1400], Kd=[120,120]. It prints separate
validity counts, pressure decomposition and force summaries. This was
checked against both files without modifying them.

## One-radian circle trials 3 and 4 — 2026-10-07

`common_pressure_testing/FB_circle_phi_1rad_trial3.mat` logs one 25.01 s
circle over about 9.99–35 s; `FB_circle_phi_1rad_trial4.mat` logs five 30 s
circles over 10–160 s. Both actually request phi=1.0000 rad and log
Kp=[3000,3000], Kd=[120,120]. The 5–10 s window previously labelled
"pre-circle hold" by the report is a phi **ramp** (mean desired phi 0.495),
not a settled fixed-theta 1 rad hold. Trial 4's 20 s phi return starts at
160 s and finishes before the 220 s feedback cutoff.

| Recording | Valid circle samples | Mean measured phi (rad) | Phi MAE (rad) | Theta wrapped MAE (rad) | Measured turns | `pressureLimited` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Trial 3, one circle | 499/500 | 0.7693 | 0.2379 | 0.3026 | 0.976 | 468/500 |
| Trial 4, five circles | 2990/3000 | 0.7425 | 0.2575 | 0.2655 | 4.977 | 2894/3000 |

Trial 4 mean phi by lap is 0.7540, 0.7434, 0.7407, 0.7386, 0.7360 rad:
the persistent bend-magnitude deficit is repeatable across laps with a small
downward drift. Orientation-binned means range from about 0.926 rad at
desired theta 45–90 degrees to 0.566 rad at 180–225 degrees. This supports
direction-dependent tracking; it does not by itself identify gravity,
actuator response, or stiffness as the cause. In trial 4, l2/l3 circle RMSE
is 2.6058/3.5583 mm. Actual phi reaches up to 1.0099 rad in trial 4 and
1.1161 rad in trial 3 at some instants, so the recorded deficit is not an
absolute inability to attain 1 rad at every orientation.

Read-only inspection of raw MAT traces shows at least one channel at the
1.2 bar lower rail in 463/500 trial-3 and 2881/3000 trial-4 circle samples;
no logged channel reaches the 3 bar upper rail. Trial-4 lower-rail counts by
P1/P2/P3 are 1422/675/811. Thus lower-rail clipping explains nearly all
`pressureLimited` samples, although that flag also covers correction clamps
and slew. The saved inverse source currently has +/-1 bar s/a correction
bounds and 2 bar/s slew. Reconstructing its unclipped target from logged
force and initial static anchor matches every valid trial-4 circle pressure
command to numerical precision. On that reconstruction, symmetric `ds`
hits its +/-1 bar clamp in 1835/2990 valid trial-4 samples; antisymmetric
`da` never hits its clamp. The raw target falls below 1.2 bar in 2873/2990
valid samples, but never exceeds 3 bar. A uniform upward pressure shift
could preserve that raw differential within 1.2–3 bar in 2583/2990 valid
samples; this is an *offline allocation feasibility calculation*, not a
validated common-pressure setpoint or physical-response prediction. In the
remaining samples, the requested pressure span can exceed the 1.8 bar
available window (maximum reconstructed span about 2.10 bar).

Trial 4 has ten isolated invalid NDI/model frames during the circles and
trial 3 has one; logged controller pressure holds through inspected invalid
frames. Controller and post-Kill NI-path pressure-equivalent logs match
exactly in both runs. Actual chamber pressure was not measured, and the
known local buckle is not captured by tip-only NDI sensing. No Simulink,
controller or analysis code was modified for this result inspection.

## Five-second one-radian circle trial 6 — 2026-10-07

`FB_circle_phi_1rad_trial6_good_5s_0.5bar.mat` logs five 5 s circles from
10–35 s, desired phi=1 rad, Kp=[3000,3000], Kd=[22,22]. The user reports
lowering the active command floor to 0.5 bar. Valid/model/NDI orientation
samples are 489/500. Mean measured phi=0.9883 rad, phi MAE=0.0485 rad,
range 0.8864–1.1545 rad, and five-lap mean phi=0.9918, 0.9902, 0.9873,
0.9865, 0.9854 rad. Measured theta turns=4.989, but wrapped theta
MAE=0.3173 rad; l2/l3 RMSE=2.7649/3.0615 mm. The phi improvement relative
to 30 s trial 4 is substantial, while directional and individual-length
tracking are still imperfect.

`pressureLimited` is true in 65/500 circle samples. Per-channel pressure
ranges are [0.8057,0.7521,0.6052] to [2.6656,2.5299,2.2885] bar.
At least one command is below 0.8 bar in 93/500 samples (P2 20, P3 73),
but none reaches 0.5 bar. The observed per-channel command slew reaches
2.0 bar/s. Controller and post-Kill NI-path commands match exactly;
chamber pressure is not sensed. This run also changes circle time and Kd,
so it does not isolate the effect of the lower bound. The analyzer's
5–10 s "pre-circle hold" window is a phi ramp, not a steady hold.

Retrospective lag analysis of trial 6 finds a 0.24 s best shift between
desired and measured theta; per-lap fits are 0.25/0.24/0.23/0.23/0.23 s.
Comparing desired theta at `t-0.24 s` lowers wrapped theta MAE from
0.3173 to 0.0692 rad. The analogous offline length comparison lowers
l2/l3 RMSE from 2.765/3.062 to 0.820/1.110 mm. This quantifies an
effective phase offset, not a proven physical transport delay or a
modified controller result. Commands reach approximately 2.0 bar/s in
65/499 adjacent circle intervals.

The later `FB_circle_phi_1rad_trial6_good_10s_0.5bar.mat` logs five 10 s
circles at phi target 1 rad, Kp=3000, Kd=22. It has 990/1000 valid circle
samples, measured mean phi 0.9634 rad, phi MAE 0.0485 rad, theta MAE
0.2713 rad, 4.978 measured turns, l2/l3 RMSE 2.3119/2.5840 mm, and
pressureLimited 5/1000. A pure time shift fits 0.40 s here versus 0.24 s
at 5 s/circle. An offline linear fit of theta error against desired
angular speed gives intercept/slope 0.196 rad/0.121 s for the 10 s run,
versus 0.150 rad/0.133 s for the 5 s run. Thus the signed theta error
contains both an approximately speed-independent offset and a speed-related
part; no unique physical cause is identified. The saved model's current
`des_phi=1.5` setting was written after this 1 rad MAT and should not be
mistaken for its logged target.

## Documentation migration checks — 2026-10-06

Historical log copied before editing with matching SHA-256
`31A09742DA0F9EDEF769C9733C263897B4CA066C44B0778180429C8B5A8C1EA7`.
Read-only FBL_commonP archive inspection finds both supplied feedback scripts;
this updates the earlier old-script observation, without verifying wiring,
compilation or recording equivalence. No experiment results were regenerated.
