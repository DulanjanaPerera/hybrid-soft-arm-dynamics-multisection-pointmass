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

## Documentation migration checks — 2026-10-06

Historical log copied before editing with matching SHA-256
`31A09742DA0F9EDEF769C9733C263897B4CA066C44B0778180429C8B5A8C1EA7`.
Read-only FBL_commonP archive inspection finds both supplied feedback scripts;
this updates the earlier old-script observation, without verifying wiring,
compilation or recording equivalence. No experiment results were regenerated.
