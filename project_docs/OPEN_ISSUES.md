# Active unresolved issues

Current as of 2026-10-08. Resolutions belong in the log and relevant results or
decisions; remove resolved items here. Prior code edits require explicit direction.

## OI-01 — Invalid-input pressure discontinuity

- Evidence: feedback recording at 14.40 s, one invalid 0.05 s NDI pair;
  approximately 1.666 bar to zero immediately, about 33.3 bar/s.
- Cause: invalid/nonfinite-input early returns reset `lastPressure_bar` and
  bypass the normal 0.8 bar/s limiter. Present in the source and the saved
  FBL_commonP inverse script inspected during migration.
- Attempted: nominal-grid V1 checks missed this transition. Paste-ready V2
  holds one brief invalid frame, latches after 0.10 s without valid data, and
  slews to zero at 0.8 bar/s. Offline MATLAB fault injection passed.
- Status/action: the October 7 higher-gain run logged one invalid frame at
  22.00 s; controller pressure and post-Kill NI command both held the previous
  trim, and the next valid frame resumed without a step. The October 8
  `integral_term/plot_data_110s.mat` run additionally logged a 0.45 s tip-0B
  dropout beginning 123.30 s. After the 0.10 s grace period, the controller
  latched the fault, slewed to zero at 2 bar/s, and stayed off when tracking
  returned. The transient hold and sustained-invalid command-path latch
  are now observed in user-recorded runs. The intermittent tip-tracking
  cause, unusable Clock exception, and actual chamber-pressure response
  remain unresolved; retain the physical Kill path.

## OI-02 — Persistent tracking error and imperfect trim

- Evidence: latest U/V step holds retain error despite pressureLimited=false;
  trim was still drifting at feedback engagement.
- Suspected cause: model/actuator/static-bias mismatch and PD-only correction;
  the recording does not isolate a cause or identify new gains.
- Attempted: stiffness sweeps, empirical phi stiffness, neutral-trim estimates,
  manual trim, local calibration, bounded incremental feedback.
- V2 recording retained Kp=[1100;1100], Kd=[100;100] rather than the proposed
  increase. U/V hold length-error norms were 0.133/0.248 mm versus
  0.112/0.211 mm in V1; starting pose/history differed, so this does not
  isolate a controller-code effect.
- The October 7 test did use Kp=[1400;1400], Kd=[120;120]. Zero-target
  length-error norm fell to 0.020–0.028 mm; U/V hold norms were 0.090/0.134 mm.
  Relative U response was +0.00320 rad, relative V only +0.00015 rad.
- Later 0.2/0.1 rad circle recordings give mean measured phi 0.1213/0.0613
  rad during the fixed-theta 30–35 s holds. At 0.2 rad desired l2=l3=1.30
  mm but measured l2/l3=1.1625/0.3410 mm; at 0.1 rad desired 0.65/0.65
  mm but measured 0.6463/0.0808 mm. The l3 error dominates both holds.
  Acceleration correction is nonzero; active-window pressure limiting is
  absent, and logged NI-path commands match controller commands.
- At desired phi=1 rad and Kp=3000, trial 4 makes 4.977 measured turns over
  five requested circles, but mean measured phi is 0.7425 rad. Limiting is
  active on 2894/3000 circle samples; at least one channel hits the 1.2 bar
  floor on 2881/3000, none hits 3 bar, and reconstructed symmetric `ds`
  reaches its +/-1 bar clamp on 1835/2990 valid samples. Bend magnitude is
  direction dependent (about 0.57–0.93 rad across desired-theta sectors).
- Later trial 6 with a user-reported 0.5 bar command floor, 5 s/circle and
  Kd=22 improves mean phi to 0.9883 rad and phi MAE to 0.0485 rad; limiting
  falls to 65/500 samples. Theta MAE remains 0.3173 rad and individual
  lengths retain roughly 2.8/3.1 mm RMSE. Because floor, speed and Kd all
  changed, the floor's isolated contribution is unknown. Commands dip below
  0.8 bar in 93/500 circle samples, but never reach 0.5 bar.
- Offline trial-6 theta alignment identifies a repeatable 0.23–0.25 s
  effective lag across all five laps. A 0.24 s comparison shift reduces
  theta MAE to 0.069 rad and l2/l3 RMSE to 0.82/1.11 mm. This does not
  isolate actuator, slew, filter or logging-time cause; the 2 bar/s command
  slew is reached in 65/499 adjacent circle intervals.
- The 10 s/circle run fits a larger 0.40 s pure time shift despite nearly
  the same phi MAE. Across both runs, theta error is better described by
  a 0.15–0.20 rad offset plus a 0.12–0.13 s speed coefficient; neither
  regression identifies the physical cause. Test timing lead separately
  from a later angular offset rather than assuming a fixed 0.24 s delay.
- In `StaticTest/changeKi_dynamically.mat`, Ki=[700,700] for complete
  20-40 s laps and Ki=[0,0] for complete 50-70 s laps, with fixed
  Kp=[3000,3000], Kd=[42,42] and nearly fixed vertical base. Ki reduced
  phi signed underbend from +0.089 to +0.010 rad, but theta and absolute
  l2/l3 errors changed little. Phi error variation and the pressure-limited
  fraction increased (81% on versus 6% off). The integral state was not
  logged; the zero-reference hold occurred only with Ki off.
- Status/action: open. The near-straight NDI angle was unobservable during
  nearly all small-step feedback samples, so inferred directional lengths
  there are limited. At 1 rad, further gain increases alone cannot resolve
  pressure clipping and may request still larger inaccessible corrections.
  Review allocator headroom offline and, if a physical test is later chosen,
  use fixed-theta holds to separate static capability from rotation lag.
  The user-observed local buckle limits assumptions about raising common
  pressure. For the Ki controller, log its state and determine whether
  internal correction bounds, pressure rails or slew cause frequent limiting;
  compare matched repeated trials and assess anti-windup before increasing
  Ki. Command-only data cannot isolate model stiffness from actuator
  asymmetry or actual chamber pressure.

## OI-03 — History-dependent operating pose and local buckle

- Evidence: manual later trial-4 mean U,V=(0.02293,0.00690) rad differs from
  subsequent local-calibration trim holds; identical trim commands gave U/V
  spans 0.0096/0.0138 rad. User observed a small local buckle.
- Suspected cause: unmodeled state/history/prestress; no unique physical cause
  established. NDI tip orientation cannot reconstruct distributed curvature.
- Attempted: repeated pressure cycles, +/-0.04 bar probes; one a+ hold had an
  unexplained 0.0063 rad orientation change near straight.
- Status/action: open; longer repeat holds from different histories, signed U/V
  and XYZ/timing checks, and visual buckle observation. No camera system exists;
  do not make it a prerequisite for all further discussion or analysis.

## OI-04 — Exact saved-model reproducibility and pressure logging

- Evidence: intermediate V2 MAT lacked controller and post-Kill logs; the
  October 7 new-gains MAT restores pressure, force, limiting, validity, Kill
  and all three NI-command logs. Read-only earlier saved-model inspection
  found V2 in main `MATLAB Function1` and V1 in `MATLAB Function8`.
- Attempted: read-only archive inspection; no compile, run or save here.
- Status/action: command-path logging is now available and exact controller/NI
  equality was observed in the new-gains and circle runs. The historical
  circle MAT named `FB_circle_phi_1rad.mat` logs a 0.1 rad reference; later
  `trial3`/`trial4` files do log 1 rad. The saved active one-input
  `thetaCircle(t)` now uses internal `des_phi=1`, 30 s/lap and five circles,
  while the trial-3 MAT describes one 25 s circle. The saved model was
  written after trial 4, so use each MAT's logged reference and gains to
  identify that run. Full hardware response remains unverified; chamber
  pressure is not sensed.

## OI-05 — Conditional parameters and model inconsistencies

- Evidence: nominal historical K=1350, prior recorded inverse K=250, later K=650;
  current offline runner D=20 versus common-pressure inverse D=40. Stiffness
  sweep selects 650 conditionally; phi-only fit disagreed especially on unloading.
- Attempted: scalar sweep, P1 static estimation/PCHIP, older coupled sysid fit.
- Status/action: open; compare on held-out data with consistent pressure,
  coordinate and timing conventions. Do not transplant the older fit blindly.
  Inv_Dyn_val last documented bound update was blocked at +/-30 mm; verify
  current contents before treating that historical pending item as resolved.

## OI-06 — Predictor readiness for nonlinear MPC

- Evidence: older separate scalar NARX networks trained on corresponding
  single-actuator chirps; validation uses measured past output, not autonomous
  recursive rollout. Centered label smoothing and normalization before split
  need leakage review; default predictor Ts=1 s is not a verified sample rate.
- Attempted: read-only factual audit, no current NMPC implementation.
- Status/action: open; held-out simultaneous-actuation recursive tests, correct
  sample alignment, enforce the length constraint, and audit preprocessing.

## OI-07 — Measurement and pressure-model uncertainty

- Evidence: orientation-inferred lengths, uncalibrated pollStart-to-device timing,
  provisional deadzone/damping and asymmetric equal-command response.
- Attempted: calibration/preflight, estimator sensitivity and pressure trials.
- Status/action: ongoing limitation; retain explicit validity and timestamps,
  verify mounting, and assess controller delay. Proceed under the user's chosen
  command-pressure assumption; pressure sensing is unavailable, not a blocker.

## OI-08 — Near-straight bending-plane observability

- Evidence: in the October 7 new-gains run, thetaObservable was false for
  3298/3300 feedback samples while inverseModelValid was true. During every
  steady target hold, theta stayed at 0.8166 rad. The NDI orientation estimator
  holds prior theta below the calibrated 0.5-degree straight-angle tolerance;
  measured lengths are computed using that held angle. Median feedback phi
  was 0.00270 rad, below the tolerance.
- Consequence: near-straight length-derived U/V can appear smooth but cannot
  resolve a new bending direction. This especially limits interpretation of
  the weak V response under the 0.01 rad reference.
- Status/action: open; use thetaObservable in tracking analysis, compare raw
  tip XYZ and orientation magnitude, and repeat at a bend large enough to
  make theta observable before drawing directional-control conclusions.

## OI-09 — Gravity direction when the whole arm is rotated

- Evidence: the NDI poll has an absolute base quaternion, but the current
  reader exports only base-relative bend; the common-pressure inverse fixes
  `params.g=[0;0;-9.81]`. The offline runner applies a base rotation, and
  `armS_single_core`/the RHS require careful gravity-sign interpretation.
- Status/action: open; establish physical vertical in the NDI frame, verify
  base-sensor mounting and gravity sign at known static orientations, expose
  and gate the calibrated base attitude, then feed the same arm-frame model
  gravity to all relevant dynamics paths. Do not treat recalibration in a
  rotated pose as a substitute for measuring base rotation.
