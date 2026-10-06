# Active unresolved issues

Current as of 2026-10-06. Resolutions belong in the log and relevant results or
decisions; remove resolved items here. Prior code edits require explicit direction.

## OI-01 — Invalid-input pressure discontinuity

- Evidence: feedback recording at 14.40 s, one invalid 0.05 s NDI pair;
  approximately 1.666 bar to zero immediately, about 33.3 bar/s.
- Cause: invalid/nonfinite-input early returns reset `lastPressure_bar` and
  bypass the normal 0.8 bar/s limiter. Present in the source and the saved
  FBL_commonP inverse script inspected during migration.
- Attempted: nominal-grid offline slew checks passed but did not cover this
  pressurized invalid-frame transition. No fix has been implemented.
- Status/action: open; agree on transient versus persistent invalid handling,
  implement when requested, and inject faults offline before another run.

## OI-02 — Persistent tracking error and imperfect trim

- Evidence: latest U/V step holds retain error despite pressureLimited=false;
  trim was still drifting at feedback engagement.
- Suspected cause: model/actuator/static-bias mismatch and PD-only correction;
  the recording does not isolate a cause or identify new gains.
- Attempted: stiffness sweeps, empirical phi stiffness, neutral-trim estimates,
  manual trim, local calibration, bounded incremental feedback.
- Status/action: open; assess settled repeatability and command sensitivity,
  then compare authorized gain/bounded bias approaches with rate and stability
  checks. The clamp did not explain the latest residual error.

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

- Evidence: latest MAT lacks post-Kill regulator logs. Earlier log inspection
  reported old scripts on disk; migration inspection now finds both supplied
  feedback scripts in FBL_commonP.slx.
- Attempted: read-only archive inspection; no compile, run or save here.
- Status/action: partially clarified, still open. Verify wiring, gains, reference,
  Kill, timing and logging; capture saved settings matching the next recording.
  No equivalence to the prior recorded run is established solely by script presence.

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
