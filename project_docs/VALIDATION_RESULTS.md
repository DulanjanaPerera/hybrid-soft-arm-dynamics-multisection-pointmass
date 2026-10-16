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

## Dynamic Ki switch during 1.5 rad circles — 2026-10-08

`StaticTest/changeKi_dynamically.mat` logs five 10 s circles at desired
phi=1.5 rad with Kp=[3000,3000] and Kd=[42,42]. Ki changes from
[700,700] to [0,0] at 43.31 s, inside lap 3. The base remains nearly
vertical throughout. For a phase-matched comparison, complete laps 1-2
(20-40 s, Ki on) and laps 4-5 (50-70 s, Ki off) were used; valid NDI/model
counts were 389/400 and 395/400. Lengths are orientation-inferred. The
angle error is desired minus measured, with theta wrapped to [-pi,pi].

| Metric | Ki 700, 20-40 s | Ki 0, 50-70 s |
| --- | ---: | ---: |
| Phi signed error | +0.010 rad | +0.089 rad |
| Phi MAE | 0.0917 rad | 0.1080 rad |
| Phi error standard deviation | 0.111 rad | 0.083 rad |
| Wrapped theta MAE | 0.2208 rad | 0.2194 rad |
| l2 MAE | 2.391 mm | 2.391 mm |
| l3 MAE | 3.030 mm | 2.976 mm |
| l3 signed error | +0.478 mm | +1.002 mm |
| `pressureLimited` | 1623/2000 (81%) | 115/2000 (6%) |

Ki=700 substantially reduces persistent phi underbend but does not
materially improve theta or absolute length tracking. Phi variation and
pressure limiting increase with Ki on. Nominal pressure rails are touched
for about 14.5% of command samples with Ki on and 5% with Ki off; the
limited flag also includes internal correction and slew limits, so it is
not equivalent to rail saturation. The integral state was not logged. The
post-circle zero-reference hold occurs only with Ki off and is not a
matched comparison. This single on-to-off trial does not establish a
general stability margin or isolate slow drift. No model or hardware was
run by Codex.

## Integral-run fault at 123 s — 2026-10-08

`integral_term/plot_data_110s.mat` records through 140.05 s and identifies
model `FBL_commonP_circle_BRot`. The desired phi remains 1.5 rad and the
theta reference continues past 123 s. The tip 0B validity briefly drops
at 123.00-123.05 s, recovers, then drops again for 123.30-123.70 s;
base 0A and base-orientation validity remain good. InverseModelValid is
last true at 123.29 s. The controller holds [1.891,0.646,2.262] bar until
123.39 s, begins a 2 bar/s ramp toward zero at 123.40 s, reaches all-zero
at 124.50 s, and remains latched after NDI validity returns at 123.75 s.
Controller and NI-path command logs match (max difference zero); the
logged Kill state is 1. This validates the sustained-invalid command-path
latch and recovery lockout for this user-recorded run. It does not measure
actual chamber pressure or identify why tip tracking failed.

The saved SLX StopTime is 250 s, while the inverse and integral controller
time bounds are 220 s. No 120 s MATLAB Function time gate was found;
the apparent stop at 123 s was the fault response.

## Two-run Ki angle comparison script — 2026-10-08

`integral_term/compareKiAngleErrors.m` was validated offline with
`without_i_1rad.mat` and `withi_1rad.mat`. It generated the requested
theta, phi, and angle-error figures on their shared 0-80 s interval.
The logged desired references match. Kp differs (1800 versus 3000),
while Ki differs (0 versus 300) and Kd is 42 in both runs, so these
recordings cannot isolate the effect of Ki.

## 2026-10-08 - Nominal passivity controller numerical verification

`passivity_adaptive_control/validatePacFixed.m` passed MATLAB R2025a
hardware-free checks: nine-point +/-8 mm inertia positivity/symmetry,
Mdot-2C skew symmetry (1.9626e-17 maximum residual), finite-difference dM
agreement (2.5032e-10 relative error), force summation, matched-plant
closed-loop/energy identity (1.7101e-15 N balance residual), invalid-input
and base-validity gating, and gravity-direction response. Four-second
matched-model ode15s tracking ended with 3.8647e-7 m error norm.
These checks do not cover model mismatch, pressure limits/deadzone/slew,
filtered/discrete measurement dynamics, hardware, or SLX/codegen compilation.

## 2026-10-08 - PAC pressure allocation and fault supervisor

`validatePacPressure` passed MATLAB R2025a offline tests of trim-force bias,
all-active unconstrained allocation (maximum residual 2.2888e-15 N), pressure
rails/correction clipping, nominal activation crossing, 2.0 bar/s slew,
stationary-anchor condition, changed gravity without anchor recapture,
brief-invalid hold/recovery, sustained-invalid/stale-data latch, no spontaneous
fault recovery, disable/re-enable, Clock faults and combined-wrapper execution.
Stage-1 mechanical checks also passed unchanged. Code Analyzer reported zero
messages. MATLAB Coder C source generation passed for fixed input dimensions
and logical validity/enable flags. No binary build, SLX/model execution,
hardware access, chamber-pressure measurement or sampled closed-loop physical
validation was performed. Nominal below-threshold force estimates and the
constant trim calibration require experimental evaluation.


## 2026-10-08 - PAC.slx saved wiring review (no execution)

Read-only ZIP/XML inspection of the user-created PAC.slx confirms measured
[q;dq] (4), desired [qd;dqd;ddqd] (6), gravity_arm, orientationValid_control,
baseOrientationValid, conditional measurement timestamp, Clock and boolean
true enable are connected to the eight pacCommonPressureBlock inputs in order.
The sensor uses AuroraReader_s2_BaseRotation, COM12/921600, 0.05 s, the existing
20260930 calibration (file exists), and reference gravity [0;0;-9.81].
Pressure connects directly to the retained Kill switches, 0-3 bar saturation,
10/9 voltage conversion and NI outputs. Kill has no saved Value override
(default Constant value 1, which selects pressure through these switches).

Outstanding integration findings: controller and validMeasurementTime have no
saved explicit discrete rate and receive continuous Clock; configure both
Discrete 0.05 s to avoid continuous persistent-state execution and false
timestamp freshness between polls. Solver is VariableStepAuto. Active circle
reference is phi=1.5 rad, 20 laps at 5 s/lap; phi ramp begins at 15 s and
circles at 20 s, overlapping the earliest 15-20 s feedback blend. Hold the
reference through trim/engagement (status=3) before starting its ramp.
The XYZ scope desired input has no source. Controller diagnostics are named
and terminated, but no DataLogging enable entries were found in saved system
XML; naming and global logsout enable do not select individual signals.
Enable signal logging explicitly for controller outputs, inputs, validity,
gravity/timestamp and post-Kill pressures. Ensure root and isolated PAC folder
are on MATLAB path. No SLX compile/run or hardware was attempted; no model
or controller source changed. This remains the fixed-parameter passivity
baseline, with adaptation pending.


## 2026-10-09 - First PAC circle recording: feedback authority too low

Offline MATLAB read of PAC_circle/1rad_10s_1cycle.mat (78 logged signals,
0-61.1 s). Supervisor trim starts 10 s, anchor captures 15 s, blend finishes
20 s; one invalid hold 30.60-30.65 s, no latched fault. Circle 25-35 s:
phiOrientation mean 0.016283 rad (range 0.011401-0.019246), desired phi 1 rad.
All pressure commands remain above activation; limited false throughout circle.
Controller pressure equals post-Kill NI-path commands; ranges P1 1.6614-1.6794,
P2 1.6306-1.6817, P3 1.4452-1.4974 bar. Valid model 199/200 circle samples.
Feedback force components bounded approximately +/-0.295 N; inertial terms
approximately +/-0.028 N. Valid allocation residuals are roundoff (~4e-15 N).
Therefore no evidence of command gating or saturation as cause of poor bending.

Source nominal gains lambda=2 /s and Ks=10 N s/m give restoring coefficient
Ks*lambda=20 N/m (plus M*lambda velocity feedback), much weaker than nominal
650*[2 1;1 2] N/m stiffness. At qd=[6.5;6.5] mm, q near zero, static feedback
is about [0.13;0.13] N versus nominal stiffness force [12.675;12.675] N at the
desired pose. The compensation uses measured q, so near straight it supplies
little stiffness compensation; finite force-feedback gains must drive initial
error reduction. Initial gains were offline matched-model starting values, not
experiment-tuned gains. Increase/tune restoring authority with damping and
allocation limits considered; adaptation is still absent. No source/model edits
or hardware/model execution were performed.


## 2026-10-09 - Reusable PAC static analysis and first static test

Added PAC_circle/analyzePacStatic.m and local README at user request. Offline
only: automatic longest nonzero constant-reference hold, latter-half settled
metrics, inferred Ks/lambda from recorded forces/sliding signal, engagement/
fault/validity/pressure rails/freshness, post-Kill command comparison, force
components/residual, tracking plots, text report and compact analysis MAT.
Explicit hold window handles repeat trials. Gains are identified under settled
desired-rate assumptions rather than read from current controller source.

Validated on static_phi_test1.mat: 0-99.25 s; plateau 30.05-60 s; settled
45.025-60 s. Desired/measured phi .2/.01539 rad; inferred Ks=[10 10],
lambda=[2 2]; feedback means [.02254 .02947] N. Mean length errors
[1.1277 1.4707] mm. Anchor 15 s, full control 20 s; no latched fault or
limiting in hold; one invalid hold sample. NI-path command difference zero;
residual roundoff, timestamp age 0-.05 s. Weak restoring authority persists.
No controller or SLX edits; no hardware/model execution.


## 2026-10-11 - PAC static gain-pair comparison (user-supplied reports)

Compared lambda=400,Ks=10 with lambda=200,Ks=20; both Ks*lambda=4000 N/m.
Reported settled phi .15444 vs .15214 rad, errors .04556 vs .04786 rad.
Length errors [.21436 .38031] vs [.23269 .39172] mm. New run has no hold
limiting/fault/invalidity and controller-to-NI command difference zero.
Settled inertial-force peaks fall [.11041 .085709] to [.021285 .016562] N
(about fivefold); lower measured motion also contributes, so this is not a
controlled proof of gain effect. Gravity and captured forceBias differ between
runs. Static error remains near 24 percent, supporting a persistent mismatch
rather than allocation clipping. Keep 200/20 as provisional baseline; repeat
static holds/directions to check consistency before adaptive compensation.
Evidence is pasted reports/screenshots, not independently reloaded raw MAT.
No controller or SLX edits.


## 2026-10-11 - Direction-dependent PAC static performance

User-supplied reports/screenshots for static_phi_lambda_200_damp_20_theta_Piby2.mat
and static_phi_lambda_200_damp_20_theta_Pi.mat; raw MAT not independently
reloaded for this comparison. All use phi=.2 rad, lambda=200 /s, Ks=20 N s/m,
Ks*lambda=4000 N/m, hold 30.05-60 s and settled latter half.

| Desired theta | Settled phi | Phi MAE rad | Theta wrapped MAE rad | Mean l2/l3 error mm |
| --- | --- | --- | --- | --- |
| 0 (previous run) | .15214 | .04786 | .04643 | .23269 / .39172 |
| pi/2 | .19065 | .00935 | .16734 | -.34171 / -.071021 |
| pi | .16713 | .03287 | .03166 | -.15466 / -.27375 |

Both new holds have full engagement/anchor, 600/600 valid samples, zero limited
samples, all commands above activation and zero controller/NI command mismatch.
Gravity is similar to theta=0 trial. pi/2 feedback mean [-1.3668 -.28411] N;
pi feedback mean [-.61858 -1.0949] N. Static allocation residuals are roundoff;
pi/2 ramp/return shows oscillations and a transient residual outside the hold.
Good pi/2 phi alone hides a ~9.59 degree direction error and unequal length
errors. pi theta error ~1.81 degrees. Captured forceBias differs by trial.

Interpretation: persistent mismatch is direction dependent, not explained by
hold pressure clipping. Mechanics/stiffness/coupling, pressure-force mapping,
trim bias and orientation geometry remain candidate causes; no unique stiffness
identification from these tests. Keep 200/20 provisional baseline; formulate
adaptation in both length coordinates. Distinguish additive force offset from
stiffness/coupling mismatch using repeats/multiple amplitudes before selecting
adaptive parameters. No adaptation implementation authorized by this discussion.

Prior gain discussion: lambda converts length error to reference velocity (1/s),
Ks converts sliding velocity error to force (N s/m); direct feedback has
position coefficient Ks*lambda and velocity coefficient Ks. Exposing gains as
2x1 block inputs is possible but was discussed only, not implemented by Codex.
Constant gains per run preferred initially; abrupt lambda changes alter
reference dynamics and constant-gain stability arguments do not automatically
apply. Physical params.D remains distinct from controller Ks.


## 2026-10-11 - Combined PAC static sweep implementation

Added PAC_circle/pacStaticSweepReferenceBlock.m: 9 direction/amplitude trials
plus theta=0/.2 repeat, start25 s, 5/30/5/10 s rise/hold/fall/zero-rest, end525 s,
Stop Time535 s. Readiness must reach status3 before start. Zero-reference rest
does not vent pressure or recapture trim. User pastes block; no SLX edits.
Added analyzePacStaticSweep.m: detects >=10 s logged constant nonzero holds,
reuses per-hold static analyzer without intermediate artifact overwrite, prints
compact tracking/gains/force/pressure/quality reports and saves summary CSV, TXT,
MAT. Exploratory unconstrained affine residual fit is not physical stiffness
identification. analyzePacStatic adds optional fourth writeArtifacts flag.
Validated analyzer on existing single-hold recording and all ten reference
plateau IDs/amplitudes in MATLAB offline. Full multi-hold physical validation
is pending; no hardware or Simulink run performed.


## 2026-10-11 - Ten-hold PAC sweep: repeat changes and affine fit limits

Evidence: user-pasted static_phi_lambda_200_damp_20_sweep.mat sweep report and
summary screenshot; raw MAT not reloaded. All settled holds valid/engaged, no
limiting, above activation; lambda200/Ks20. Theta0 phiActual .07259/.14517/.22800
for desired .1/.2/.3; pi/2 .09539/.18283/.27894; pi .04626/.12473/.21352.
Same theta0/.2 repeated at hold10 gives .10945 vs .14517 at hold2 (24.6%
decrease); feedback [2.2642 2.4454] vs [1.0893 1.7723] N, despite fixed
forceBias [-.69562 -2.9912] N. Changes demonstrate time/history dependence
or another changing condition; cannot isolate hysteresis, creep, air supply,
orientation drift or temperature from report. Direction and trial order confound.
pi/2 theta MAE .54541/.27384/.20097 rad at increasing phi, so phi magnitude
alone overstates accuracy. pi/.1 theta omitted by observable-angle threshold.

Affine fit deltaK=[1159.7 -22.415;1241 987.44] N/m, b=[-.045984;-.45184] N,
training RMSE [.68587 .65032] N, predictor condition2.338. Full rank does not
remove history confounding; strongly nonsymmetric fit is not a conservative
stiffness correction. Do not insert into EoM. Next repeat same theta0/.2
several times with unchanged anchor and track zero-reference baselines; separate
repeatability/history from fixed configuration dependence before selecting
stiffness versus adaptive bias parameters. No code/model changes.


## 2026-10-11 - Adaptive PAC force/pressure path implemented

User authorized implementation in passivity_adaptive_control only. Added
pacAdaptiveConfig/InitialState/History/Force/Update/CommonPressureBlock,
validatePacAdaptive, validatePacAdaptiveTracking, replayPacAdaptive,
analyzePacAdaptive and ADAPTIVE_INTEGRATION.md. Existing nominal sources/SLX
files preserved. New wrapper has previous 8 inputs plus 2x1 lambda, 2x1 Ks,
boolean enableAdaptation; previous17 outputs plus Khat,Dhat,hysteresis force,
bias,10parameter estimates,3history states,active/bound flags,adaptive force.

Khat is spring-structured SPD: [k1+k3 k3;k3 k2+k3], positive independent
springs/nonnegative shared coupling; Dhat diagonal positive. Parameter bounds,
rate caps/deadband and scale-adjusted gradient in config. Three normalized
Bouc-Wen feature states use physical muscle velocities T*dq, with fixed internal
coefficients. Their nonnegative force weights and two signed bias coefficients
adapt. Does NOT adapt alpha/beta/gamma internally, nor M/C/G. Paper online
adaptation is K/D, not inertia/hysteresis-internal coefficients.

Adaptive controller uses Dhat*v instead of old D*dq, giving a changed nominal
error response (initial (Dhat+Ks)*lambda=12000 N/m at lambda200/Ks20).
Adaptive on/off comparison must use same new wrapper. Updates apply NEXT
sample, once per fresh timestamp, only valid status3/no limiting/residual
<=.02 N; dt gaps>.075 s skip updates. Freeze retains estimates/history;
restart simulation resets. Physical Kill and pressure supervisor unchanged.
Projection bounds parameters, not physical closed-loop stability. Fixed-feature
history, saturation, sampling/rate caps/deadband mean ideal asymptotic proof
cannot be claimed for hardware. Static damping/bias/pressure mapping confound.

Verification: algebra/regressor/force sums, SPD projections, descent sign,
ideal unprojected energy identity, history reversal, startup/fresh-only update,
invalid hold/fault checks pass in MATLAB. Ideal continuous force mismatch
simulation at lambda200/Ks20 (30s) improves final10s RMSE from
[.10960 .12827] mm to [.05988 .07452] mm. This is a chosen representable
plant, not physical validation. Replay of existing 10701-sample sweep finite,
7633 update samples, zero fault samples, pressure maxima [3 2.1240 2.3221]
bar; record-fixed replay is not closed-loop prediction. C source generation
passed with fixed dimensions/logical flags; no binary/SLX build or hardware
run. Initial generation commands failed due literal temporary output path;
corrected functional codegen call passed and generated workspace dirs removed.
Parameter estimates are bounded effective coefficients, not unique physical
identification. Integration/model wiring and experiment checks remain pending.

Adaptive analyzer report and plots additionally verified with synthetic
Simulink signal Dataset logs generated by offline wrapper calls.


## 2026-10-11 - First adaptive-wrapper physical tests: oscillatory baseline

Read three saved PAC_phi_0.2 recordings offline (statis, enable_inMiddle,
enable_fromBegining). During45-60s desired.2, measured phi means
.4226/.4300/.4298 rad, ranges .058-.694/.033-.696/.052-.696. Limited
fractions .897/.934/.910. Dhat*v force peaks near60N, feedback near30N,
request peaks near85N; allocation residual peaks near90N. Learning disabled
already oscillates, so adaptive gain tuning cannot explain baseline behavior.
Middle enable yields zero updates; beginning enable yields112 updates mostly
20-26s and tiny changes, then frozen. Gates intentionally inhibit learning
when limited. Limiting includes correction clipping/slew; pressure below3bar
does not mean force request realizable.

Important implementation/tuning correction: Dhat*v adds Dhat*lambda position
feedback, giving initial12000N/m at lambda200/Ks20, vs old direct4000N/m.
Offline continuous-force simulation did not validate this setting against
50ms sensor/filter/pressure delays. Nominal gain carryover unsuitable.
Proposed next new-wrapper fixed-gain baseline lambda50/Ks40 gives initial
(Dhat+Ks)*lambda=4000N/m and less lambda acceleration amplification.
Check learning-disabled response before adaptation; do not remove freeze gate
or increase adaptation rates/pressure limits to address oscillatory baseline.
No code/SLX edits made for this diagnosis.


## 2026-10-12 - Reduced-gain adaptive comparison

Offline read PAC_phi_0.2_static_reducedGains_engagedMiddle.mat and engagedFalse.mat.
Gains lambda50/Ks40. During45-60s, desired phi.2, adaptive/fixed actual
.16658/.15206 rad; MAE .03342/.04794 (~30.3% reduction). Mean length errors
[.14369 .29258]/[.22465 .40123] mm (~36%/~27% reduction). Feedback means
[.28708 .58477]/[.44933 .80235] N. No hold limiting in either; full run16
limited samples each, zero faults. Adaptive phi range .15499-.17547 versus
fixed .14910-.15659; learning continues, not fully converged. User reports751
updates; final spring coefficients [677.45 699.57 810.52], damping40.012/40.024,
small history weights, bias[-.00927 .12572] after return. Final values after
return differ from hold values; not physical parameter identification.
Keep gains/rates unchanged initially; longer hold or full sweep now useful
for directional/history tracking evaluation. No code/SLX changes.


## 2026-10-12 - Adaptive 10-hold sweep: good .3 rad, limited transfer

User supplied PAC_sweep.mat reports/screenshots, 10701 samples0-535s,
lambda50/Ks40,9601 adaptive updates,20 limited samples total,no faults.
All reported hold latter halves valid/engaged/no limiting/above activation.
At desiredphi.3 theta0/pi2/pi actual .30369/.29674/.29744; phiMAE
.00445/.00326/.00497 rad; thetaMAE .01817/.02628/.01833 rad. At .2 first
theta0/pi2/pi actual .18019/.19137/.18679. Finaltheta0/.2 actual .13643
vs initial .18019 (24.3% decrease), error .06357 vs .01980. This is
limited transfer/readaptation/history performance, not global convergence.
At pi2/.1 thetaMAE .46708rad despite good phi. Effective spring/history/bias
estimates bounded; D estimates only40.232/40.336 versus initial40, so no
evidence of uniquely identified physical damping. Final parameters are
compensation estimates; all nominal-baseline comparisons across gain/form
changes are confounded. Proposed application-relevant next amplitude .5rad
with unchanged gains/rates,10s ramp,60s hold,include return to earlier .3/.2
to evaluate transfer. Keep constraints and observe direction/length errors.
No controller/SLX/source changes.

## 2026-10-12 - High-phi adaptive sweep: unloading onset and slew-limited oscillation

Offline raw-data review of passivity_adaptive_control/PAC_sweep_high_phis_unstable.mat,
10701 samples, 0-535 s. Reference amplitudes .5/.8/1.5 rad. The onset at110 s
is the theta0/.8 unloading ramp110-115, not the 1.5 command (starts125).
During95-109 phi .804-.812 (mean.808), all orientations valid, no limiting,
force allocation residual roundoff. During110-115 adaptation active94%;
spring estimates [568.1 804.9 824.8] at110 become [467.7 682.6 397.9] at115;
k3 falls52%. K remains positive definite (eigenvalues678/2344 to561/1385 N/m).
First limited113.950 s, first invalid orientation115.700 s, first below
activation128.500 s. During115-125 desiredphi0, actual .0186-.462 rad,
SD.1008 rad, approximate dominant frequency1.9 Hz; limited and 2bar/s slew
hit86.5%, valid92%, adaptation active5.5%. Pressures minima[1.637 1.037 .9104]
bar: above .8 activation during this first oscillatory zero interval.
Feedback force swings about -12.7..11.9 N; requests -20.2..19 N and allocation
residual about +/-20 N. Before onset feedback magnitude was below.311 N.
Later pi2/.5 hold195-210 meanphi.515 conceals range.114-.889, limited94%.
Estimates remain bounded, damping changes small, zero latched faults. Later
activation residual can freeze learning even limited=false; these flags are
not equivalent. User-reported1537 updates/7232 limited samples.
Leading hypothesis: unloading triggers delayed force-feedback/pressure dynamics
interacting with slew limiting; rapid transient adaptation may contribute.
Timing excludes invalid orientation or activation crossing as initial trigger,
but does not prove adaptation causal or distinguish pneumatic/sensor delays.
Next diagnostic priority is transient plots/time-window metrics and a matched
frozen-estimate unloading comparison; consider slower10-20 s returns and learning
gating during return/recovery after review. Do not increase slew or remove gates
on this evidence alone. No controller or SLX changes, no hardware run.

## 2026-10-12 - High-phi repeat trial confirms unloading failure

Offline comparison PAC_sweep_high_phis_unstable_trial2.mat (0-193.7 s)
with trial1. Before110, trial2 phi .80527-.81167, SD.00102, no limiting.
Return starts110 s; first limited113.800 versus113.950 trial1; first
invalid orientation114.700 versus115.700; below activation128.300 versus128.500.
Oscillation grows approximately113.6-115 s; no mathematically unique onset
sample is asserted. During115-125 phi .01599-.45883, SD.09137, limited88%,
valid90.5%, residual peak20.347 N; trial1 corresponding phi .01865-.46166,
SD.10079, limited86.5%. Strong repeatability ties failure to .8 unloading,
not random startup or initial deadzone crossing. Slew limiting remains a
constraint/possible amplifier, not proven original cause. Adaptation-versus
sampled feedback/pneumatic delay causality remains unresolved. Hold4 meanphi
.5248 conceals MAE.2455. No controller/SLX changes or hardware execution.

## 2026-10-12 - Freeze-adaptation recordings
Offline changing_enabling_adaptive_tr1/tr2.mat reviewed. Both adapt: last active
104.75/102.90 s, last parameter change104.76/102.91 s. Trial1 is NOT an
adaptation-off-from-start baseline (1106 updates and changed parameters).
Trial2 quiet105-110 phi .8077-.8099, no limiting; frozen parameters then
unloading110-115 followed115-125 phi .0606-.4638, SD.0848, limited87%.
Thus continued parameter updates during unloading are not necessary for failure.
This does not exclude learned model/history dynamics as contributors.
Both logs .01 s vs earlier .05, measurement timestamp unique increments.05 s.
Possible multirate supervisor/slew effect: initial .5 hold limited~62-65% but
residual only~.005 N, so raw limited counts are not comparable across rates.
Need verify controller execution/sample holds before retuning. Trial1 already
oscillates95-105 phi .8324-1.045, whereas trial2 quiet. Trial2 1.5 hold mean
near desired conceals135-150 range1.438-1.625 and88% limiting. No source/SLX edits.

## 2026-10-12 - Next isolation experiment and verified constraints
Source review: pressure rails .2-3 bar, s/a correction +/-2.5 bar, slew2bar/s,
nominal activation .8 bar (minimum below activation), fixed mean1.6. No timeout.
Invalid grace.10s, max sample age.075s. Allocation assumes linear force above
activation; crossing activation yields residual without necessarily limited=true.
At static target tau includes Dhat*v-Ks*s = (Dhat+Ks)*lambda*error-Ks*dq.
Current D40/Ks40/lambda50 gives4000N/m restoring and40Ns/m explicit measured
velocity coefficient, plus M*a/C*v. Frozen adaptation does not freeze history z.
Hysteresis is plausible but not established. Next verify controller sample time
.05s consistent with fresh sensor data; then short identical .8 ramp/hold/return,
frozen estimates before return, compare Ks40 versus20 with lambda50 unchanged.
This reduces both feedback stiffness and velocity authority; cannot isolate pure
velocity causality. Retain limits; log rate/rail/correction reasons separately
when implementing diagnostics. No code changes in this discussion.

## 2026-10-13 - PAC .05s execution confirmed; severe wall-time slowdown
User sample-time overlay confirms PAC outputs .05s, NI path .01s. Pacing off;
user reports .05 simulation seconds per wall second. Profile attempt aborts in
AuroraReader_s2_BaseRotation.setupImpl preflight:100/100 valid paired,88 overruns,
median interval1.3434s, reason polling slower than configured period. Screenshot
PAC call count0: profiler cannot assess PAC runtime from this run. Reader source
preflight occurs before stepImpl; physical polling paced via tic/toc in stepImpl.
Thus this failure establishes slow acquisition/setup under profiling, not PAC
computation cost. Profiler overhead versus serial/host delay unresolved. Next
check preflight without profiler, preferably standalone zero-pressure acquisition;
retain .05 period/thresholds and do not bypass preflight. No code/model edits.

## 2026-10-13 - Restart resolves NDI slowdown; new frozen high-phi onset
User restarted MATLAB and acquisition now works; PAC inherited rate retained.
Offline static_phi_test1_unstable.mat: last parameter change102.710s, despite
reported disable110s (updates may naturally gate earlier).105-110 phi.8088-.8108
quiet no limiting.115-125 desiredzero actualphi.1398-.2195,SD.0185,limited2.5%,
valid100%: residual bend rather than previous large return oscillation.
Next1.5 ramp125-130: first limiting127.600, invalid128.000, activation crossing
128.150.128-130 limited77.5%, residualpeak20.48N;130-140 phi1.3395-1.5985,
limited85.5%,allAboveActivation0. New onset tied to larger loading ramp, not
continued parameter updates. Deadzone-map/allocation interaction plausible,
not proven initial cause (limiting/invalid precede crossing). Theta nearphi0
ill-conditioned; length/component tracking required to establish true lateral
oscillation. Prior pure unloading mechanism not universal across runs. Next
isolate amplitude at .8 vs1.0 with identical frozen setup and rates, avoiding
simultaneous gain/ramp changes; exact source limits remain .2-3bar/slew2.
No code/model changes.

## 2026-10-13 - Commit preparation and offline validation
User requested committing pending updates on codex/passivity-adaptive-controller.
Included PAC implementation/models, experiment recordings/reports, plotting updates
and project documentation. validatePacFixed, validatePacPressure, validatePacAdaptive
and validatePacAdaptiveTracking passed offline. No hardware models executed; these
checks do not establish physical sampled-loop stability. No remote push requested.

## 2026-10-13 - 1.0-rad target does not prevent frozen-loop oscillation
Offline static_phi_test1_unstable_tr2.mat reviewed. Last parameter change99.160s.
105-110 quietphi.8054-.8078, no limiting.110-115 unload limited15%, residual
peak10.657N;115-125 desiredzero phi.0618-.4539,SD.0892,limited88.5%,valid91%.
All commands above nominal activation in these intervals (min .9215bar during
115-125). Failure already present before third target1.0 starts125, so reduced
third amplitude cannot explain onset.130-140 phi.9784-1.2321,limited86%,
allAboveActivation100%,residual17.89N. Continued updates and nominal activation
crossing are not required. Previous run had much quieter return despite similar
frozen estimates; calibration bias/history/physical state differ. Leading issue
remains sampled feedback/pressure dynamics interacting constraints, not proven
hysteresis. Next restore verified explicit .05 PAC execution after MATLAB restart
and confirm wall time; then compare Ks20 vs40,lambda50,identical reference/freeze.
No source or model changes.

## 2026-10-13 - Trial3 .05s recording: timing alone insufficient
Offline static_phi_test1_unstable_tr3.mat,3668samples0-183.35 implies.05s
logging; report lambda50/Ks40 unchanged. Last parameter change102.65s.
105-110 quietphi.8088-.8102,no limiting.115-125 desiredzero phi.1069-.4612,
SD.0866,limited85%,valid90.5%,aboveActivation100%,residualpeak18.30N.
130-140 actualphi.9731-1.2578,limited90%,aboveActivation100%. Thirdtarget1.0
holdMAE.1653 vs trial2.1603; failure persists at .05 logging, controller rate
requires compiled verification. Quiet static micro-limiting flags disappear,
but actual transient force deficits persist. Ks20 comparison still not done.
Next change Ks40 to20 only with verified.05 timing and identical reference/freeze;
not simply more adaptive updates or removing pressure constraints. No code edits.

## 2026-10-13 - Ks20 removes observed oscillation, residual tracking errors remain
Offline static_phi_test1_20damping_tr4.mat:lambda50/Ks20,.05s logs,last parameter
change103.70s. Compared tr3Ks40,115-125 phi now.2033-.2573,SD.0119 vs.0866,
no limiting/invalidity and residual roundoff vs85% limiting/18.30N residual.
Return remains biased ~.2rad for desired0; suppression is not zero tracking error.
130-140 phi.8898-.9429,SD.0103,no limiting. Report1rad holdMAE.04165 vs.16528
tr3 (~75% reduction);thetaMAE.00186 vs.19117. Only16limited fullrun, allholds0.
pi2/.5 holdthetaMAE.51148,lengtherrors[-1.9865 -1.5063]mm despitephiMAE.04893:
directional compensation remains poor with frozen parameters. Source force law
Dhat*v-Ks*s changes both velocity coefficient (40->20Ns/m) and static corrective
coefficient (Dhat+Ks)*lambda ~4000->3000N/m; improved physical response does not
prove pure damping or identify hysteresis. Supports excessive sampled feedback
and slew interaction. Hold Ks20/lambda50 as provisional baseline; repeat same
run before adaptation-on full directional sweep or faster trajectories. No code
changes; results not a physical stability proof.

## 2026-10-13 - Full adaptive Ks20 sweep: accurate holds, recovered transient
static_phi_test1_20damping_full_tr5.mat,.05 logs0-535s,lambda50/Ks20,
7051updates288limited,no faults. All latter-half holdvalid100%,limiting0;
phiMAE .00208-.02109 and thetaMAE .00177-.03823. Directionpi2/.5 now
phiMAE.00208/theta.03758 versus frozen trial4 .04893/.51148 (separate-run,
not isolated proof). Repeattheta0/.5MAE.01139 vs.00867 initial.
Raw260-265 pi2/.8 unloading,265-275zero,275-2801radloading.245-260quiet.
265-275phi.0091-.3703,SD.0892,limited79.5%,valid95%,aboveactivation100%,
residualpeak12.059N.275-280limited79%,residual15.583N.280-285limited2%,
295-310none,phi.9774-.9900. Recovery after next nonzero target is not proof
zero-reference oscillation self-decays; reference changed275s.
k3 reaches0 lower projection bound during run, does not imply physical muscle
stiffness0; K can remain SPD via k1,k2. Nominal activation residual~.045-.048N
in1rad holds exceeds.02 adaptive gate and may freeze learning without limited.
Hold gains/rates; next identical fullrepeat or isolated extendedzerohold at265
before next ramp tests persistence/recovery. No faster circles or gain increase
until transient resolved. No code/model edits.

## 2026-10-13 - Raw k3/velocity/force analysis of full tr5
k3 first reaches0 at171.550s; stays0 throughout245-310 including quiet
245-260 and recovered295-310. Thus0 is not sufficient cause or sudden260trigger.
K eigenvalues at260 [462.25 1057.2],270 [285.48 826.74],280[225.31 809.65],
295[418.33 1029.7] N/m remainpositive. k2 drops462->225 overtransition.
245-260 velocityRMS[.01229 .01460]mm/s vs265-275[6.536 9.491]mm/s;
velocitypeaks~15/21mm/s.265-275 dampingForce(Dhat*v) min[-6.66 -9.503],
max[8.212 9.127]N; feedbackmin[-3.446 -5.069]max[4.223 4.881]; inertia
peaks<.784N; staticpeaks<3.54N. Force requestmin[-9.027 -14.98],
max[7.702 11.32].Large D*v is not directlyphysicalvelocity damping:
v=dqd+lambda*(qd-q),hence zero-referenceD*v=-D*lambda*q; sensitivity to
position error (reference velocity), feedback alsousesmeasuredvelocity.
Dominant force swings are damping-compensation plusfeedback, not inertia or
lossSPD. Measurements filtered/inferred cannot distinguish actual oscillatory
motion from derivative amplification/delay alone; length/pressure phase review
needed. Negativek3 unsupported; next isolate reference/velocityforce sensitivity
via gains/lambda comparison or simulation, do not blindly alterterms andclaim
passivity. No controller/model edits.

## 2026-10-13 - Lambda25/Ks20: oscillation suppressed, adaptation settling slower
Read full pasted report and raw static_phi_test_20d_25Lmabda.mat0-512.85,
17limited6503updatesno faults; all evaluated holdsvalid100%,limiting0.
Reference changed vslambda50:rise/fall1s,42s trialspacing vs5s/50s before,
calibrationbias also different; cannot attribute every improvement solelylambda.
pi2/.8 unload224-225,zero225-235 phi.1162-.2113 SD.0158,no limiting,
compared previouszero265-275SD.0892/79.5%limiting.1radpi2 hold251-266
phi.9503-.9626,MAE.04293; allocationresidual.0407N exceeds.02 adaptationgate.
Theta0/1rad residualpeak.8635N,P1mean.7403bar implies nominal activationmapping
limits force fidelity despite limitedfalse. Repeattheta0/.5 hold404-434 actual
starts.051-.085 during404-409,419-434range.174-.508 andMAE.17078, not just
fullysettled offset: directionalreadaptation/transfer poor. Last updates508.15.
Lambda reduction lowers(D+Ks)*lambda approx3000->1500N/m; measuredvelocity
coefficientKs20 unchanged. Conservative baseline25/20; next matchedreference
comparison lambda30 or35 only, not higher adaptation gains simultaneously.
Settling report must distinguish convergence to target vsdrift/model mismatch.
No source/model changes.

## 2026-10-13 - Lambda40/Ks20 improved tradeoff, brief return transient
Read report/raw static_phi_test_20d_40Lambda.mat,.05logs0-484.65,same1s ramps/
42s trialspacing aslambda25.58limited6036updatesnofaults;all evaluatedholds
valid100%,limiting0.225-235zeroafterpi2/.8 return phi.014-.2092,SD.0213,
limited2%,valid99%,residualpeak.465N: brief transient rather than sustained
lambda50 returnoscillation.1radtheta0MAE.0150 vslambda25.0401;
pi2MAE.0253 vs.0429. Repeattheta0/.5MAE.0183 vs.1708;
419-434phi.496-.5256,SD.0059 vslambda25.174-.508,SD.108.
Fullresidualpeak2.424N vs15.583lambda50full, .8635lambda25.
Initialstates/calibrationdifferent;singletrialnotstabilitymargin proof.
Provisional workingbaseline40/20 withsame.05timing;repeatbeforefastertrajectories,
retainbounds andrates. Activationmapresidual persists1radholds(.269/.055N),
canblockadaptationdespite nolimiting. No code/model edits.

## 2026-10-13 - Adaptive circle1rad30s per lap
Raw PAC_circle/AdaptiveController/1rad_60s.mat reviewed; filename60s but logged
reference theta speed.20944rad/s confirms30s period,7complete laps20-230 and
partial8 through235.7. phi1rad. Lap1MAEphi.0453 theta.2348(lengthRMSE2.079/
2.272mm) tolap7phi.0305 theta.1954(RMSE1.888/1.800mm). Meanphi~.984;
valid>99%,laplimited~.7-2.5%,no faults perreport. Allaboveactivation only65.7-
73.9%: nominal mapallocation errors persistwithoutpressureclamping and canfreeze
learning. k2floor50 for1023samples (~51.15s),k3zero52samples;coefficients cycle
withdirection andhistoryweightsdrift upward,not physicalparameter convergence.
Ks inferred20; staticlambda inference is invalid oncircle (nonzerodesiredrates),
do not report approximate40.8 as verifiedinput. AngleMAE.195rad~11.2deg with
phierror~.03; direction dominates lengtherror. Effective lag estimate fromMAE
alone not uniquelyidentified. No source/model edits. Keep30s baseline; inspect
signedphaseerror/request-vsallocatedforcebytheta before speed increase or
adaptationrate tuning. Need matchedPAC-vsFBL referencecomparison later.

## 2026-10-13 - Repeat10circle1rad30s: persistent phase error and allocation mismatch
Found PAC_circle/AdaptiveController/1rad_30s.mat matching6720samples0-335.95.
10complete circles20-320,period30s. Lap2thetaMAE.2890 tolap10.2752rad;
phiMAE.0487->.0579,lengthRMSE2.784/2.895->2.756/2.601mm; littleconvergence.
Lap3signedtheta.28584 equalsMAE:dominantlylag;equivalent1.36s at.20944rad/s,
not uniquelypuretransportdelay. Aboveactivationfractionlap2.355->lap10.278,
residualpeaks up to[3.647 5.561]N withalmostnolimitedflags. Lap3residualgate
passes~32.9% samples. k2floor50 on1555samples(~77.75s),k3zero47.
Adaptationupdates2261whilehistoryweights/damping/biascontinuechanging;notmodel
convergence. Allocationactive-map mismatchprioritybeforegainincrease: existing
fixedmeanlinearallocation allowscommandsbelowactivation thenforceevaluationclips,
nominalresidualcanblocklearning. Needpressure-vsangle/signedphaseandfeasibility
review;activationmapassumednotchambermeasurement. Endthetaerrorwhilephifalls
excludeasunobservable. No code changes.

## 2026-10-13 - Activation-consistent fixed-mean PAC allocator
User authorized pacAllocatePressure change. Exact activation-map inverse now
uses requested+forceBias to derive effective pressure differences. Preserves1.6
bar mean; when effective muscle pressure iszero, its inactive commandedpressure
can decrease below.8 to satisfyfixedmean. Candidate accepted only withinexisting
.2-3bar rails ands/a +/-2.5trimcorrectionclamps. Infeasiblecandidate retains bounded
legacyfallback,nowlimited=true fornonzeroallocationresidual; no hiddenmeanincrease
or relaxedconstraints. Supervisor2bar/slew remains andmayproduce transientresidual.
Nominal activationthreshold.8 remains assumption; doesnot identifyphysicalpressure.
validatePacPressure includes feasibleinactiveexactforceand infeasibleflag checks;
validatePacAdaptive passes. No SLX edited byagent. Repeat1rad30s baseline,
lambda40/Ks20 andrecordallocationresidual/limits/updates beforehighergainorspeed.

## 2026-10-13 - Fixed-mean activation inverse circle did not improve tracking
Raw1rad_30s_tr1.mat newallocator:lap10thetaMAE.3094 vsprior.2752,phiMAE.0888
vs.0579,RMSE3.141/2.821 vs2.756/2.601mm. Newlimited semantics includes
forceinfeasibility:cannotcompare3184vs25 asphysicalclampingcounts. Lap10limited
56.2%,activationabove20.2%. Offlinecandidateconstraints50-320 valid5384samples:
minimumviolations2244,maximum1486,scorrection496,acorrection0(overlapcounts).
Therefore fixedmean1.6 cannot exactlydeliver manyrequestswithinallbounds;
activation consistency improvement didnotresolvephysicaltracking. Needbounded
variablecommonpressure allocation/closestfeasibleforce solver to avoid legacy
fallback mismatch. Do notrelax rails or assumehighercommonalwayshelps. Proposed
next allocator changes separateauthorization;keepgains/rates forcomparison.
No controller edits in this review.

## 2026-10-16 - Documentation clarification and commit preparation
Pressure nullspace doesnot require3mechanicalcoordinates:retain2lengthcoordinates
and3pressureinputs;optional pneumaticpressurestates modelvalvedynamics. Fixedmean
allocation contributesforce mismatch;variablemeanproposal remainsunimplemented.
Wrappedtheta error blocksnippet supplied inchat atan2(sin(des-meas),cos(des-meas));
not added asstandalonefile;theta nearstraight poorlyobservable.
User requests allpendingupdates committed. Includesworkspaceanalyzers,PACplotting,
activationallocator,models,setup scripts,recordings/reportsandlogs. Offline
validatePacPressure/validatePacAdaptive passed; no hardwaremodel executed.
Allocator circle limitation retained/documented,not presented asperformancefix.
