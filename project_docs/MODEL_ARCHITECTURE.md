# Current model and software architecture

Current as of 2026-10-06; numeric settings are model-specific.

## Physical system and coordinates

One hanging module, stiff backbone, three PMAs spaced around radius 13 mm,
active arc length 174.11 mm, nominal section mass 0.1 kg. Independent bending
coordinates are length changes `q=[l2;l3]` (m), with `l1=-l2-l3`; state
`X=[q;dq]`. Use length-change notation consistently: these are not absolute
actuator lengths. Pressure area uses radius 6.5 mm, separate from PMA offset.
Signed bend `U=(l2+l3)/r`, `V=(l3-l2)/(sqrt(3)*r)` remains useful near straight
where bending-plane theta is poorly observable.

## Measurement subsystem

`AuroraDriver_2.m` and `PortHandle.m` handle serial acquisition;
`AuroraReader_s2_Position_only.m` provides Simulink System-block integration;
`ndiSupportFnc.m` owns geometry, calibration, preflight, orientation/length
conversion and validation. `f20260219_2_task2config_withLenExt.m` supplies
configuration conversion and `setup_ndiConfig.m` initializes configuration.

Sensor 0A is base, 0B tip; one BX poll provides the paired XYZ/quaternions.
COM12, 921600 baud, nominal 0.05 s polling. Arm axes satisfy
`x_arm=x_NDI, y_arm=z_NDI, z_arm=-y_NDI`. Base/tip offsets are 55/50 mm.
Calibrated base-relative orientation is the primary theta/phi/length signal;
XYZ inverse kinematics is an independent comparison route. Held display values
on invalid frames are stale; validity flags must gate estimation/control.
The BX poll includes the base quaternion, but the current System block does not
output absolute base attitude. Whole-arm tilt relative to gravity therefore
does not currently reach the dynamics. A proposed extension would use the
calibrated base-sensor-to-arm rotation and a known NDI-frame vertical to form
the arm-frame gravity vector, with invalid-pose gating and a checked sign
convention before closed-loop use.
An experimental `AuroraReader_s2_BaseRotation.m` class now preserves outputs
1:19 of the original reader and adds output 20 `R_NDI_from_arm` (3x3), output
21 `gravity_arm` (3x1), and output 22 `baseOrientationValid` (logical). The
default reference-model gravity `[0;0;-9.81]` reproduces the current inverse
convention at the original vertical calibration pose; its physical sign still
needs an independent static check. Invalid frames hold the last values while
clearing validity. The saved Simulink models still use the original reader and
do not yet consume these outputs.

## Mechanical dynamics

`armS_single_entry.m` -> `armS_single_dynamics.m` -> `armS_single_core.m` and
`integratedJacobianProduct_compact.m`, `integratedPosition_nume.m`,
`christoffelSymbol.m` compute section dynamics:

`M(q)*ddq + (C(q,dq)+D)*dq + G(q) + Keff(q)*q = tau`.

K is a coupled 2x2 matrix; smooth diagonal penalties bound reduced lengths.
The common-pressure inverse currently uses K=650*[2 1;1 2], D=40*I,
gravity [0;0;-9.81], mu=2000, bounds [-0.035;0.035;1e6].
In the saved `FBL_commonP_circle.slx`, the active V2 inverse hardcodes that
gravity vector, while its `forwardModel` block separately hardcodes the same
vector and still uses K=1350*[2 1;1 2]. Base rotation currently updates
neither path. The forward block's stiffness setting is separate from the
orientation change and should be accounted for when comparing predictions.
The offline runner uses K=650, D=20, signed two-input force and world-to-arm
gravity rotation; its conventions differ from the three-positive-valve path.

## Feedback data flow

Desired theta/phi -> desired reduced lengths -> causal desired derivatives.
Paired calibrated NDI -> measured reduced lengths/validity/poll time -> causal
measured derivatives. Errors are desired minus measured. Virtual acceleration
is `v=ddq_des+Kd.*(dq_des-dq)+Kp.*(q_des-q)`; inverse dynamics computes
absolute generalized force `M*v+(C+D)*dq+G+Keff*q`.

Historical baseline filter settings: desired velocity/acceleration tau=0.10/0.15 s;
measured tau=0.20/0.30 s, nominal 0.05 s update. Verify these in the actual
saved experiment model before relying on them.

## Pressure allocation and experiment interfaces

The older inverse selects a minimum nonnegative effective-pressure triplet,
adds 0.8 bar only to active channels, and scales for a 3 bar command ceiling.
The current bounded experiment replaces that allocation with a 1.6 bar mean,
unequal empirical trim and incremental s/a corrections. Source interfaces:

- `commonPressureFeedbackInverseBlock(measured, desiredAcceleration,
  orientationValid, t)` returns pressure, desired/available forces, limited and
  model-valid flags, mass/Coriolis/static-force diagnostics in existing order.
- `commonPressureFeedbackReferenceBlock(t)` returns theta and phi for the
  240 s small-step protocol. Saved FBL_commonP contains both scripts as of migration.
- The user-maintained `FBL_commonP_circle.slx` saves V2 in its active inverse
  chart with feedback from t=0, five-second blend, +/-1 bar s/a correction
  bounds, 1.2–3 bar active-pressure clipping, 2 bar/s internal slew, and vent
  from 220 s. These current saved values differ from the older standalone V2
  source. The saved active one-input `thetaCircle(t)` now sets internal
  `des_phi=1`, `circle_time_s=30`, `num_circles=5`, starts circles at 10 s,
  and returns phi over 20 s after the last lap. It therefore has no fixed-phi
  hold between the 5–10 s bend ramp and circle motion. The user superseded
  the previously suggested external-input reference. The saved model was
  written after the latest MAT recording; confirm logged `des_phi`, circle
  count, gains and pressure command behavior for each individual run.
- Open-loop sequence/trim/manual/local-calibration functions supply scalar
  P1/P2/P3 commands and phase annotations; manual block requires downstream slew.

Commands pass through the model's gating/saturation/Kill and NI analog outputs
to physical regulators. Software pressure logs do not measure chamber pressure;
analysis assumes command equals chamber pressure at the user's direction.
An ordinary hardware-bearing Simulink run may open both serial and NI hardware.

## Software / artifact organization

| Location | Responsibility |
| --- | --- |
| Root SLX models | NDI-only, pressure-only, dynamics validation, controller variants, FBL_commonP |
| `forwardModelValidation/` | Six pressure sine recordings |
| `estimator_comparison/` | Causal derivative estimator comparison and saved settings/results |
| `inverDynamic_test/` | Static inverse/tracking analysis and original report |
| `FWD_INV_COmparision/` | Offline command replay and stiffness grid sweep |
| `planar_stiffness_test/` | P1 static fits, PCHIP evaluator and paste-ready block |
| `feedbacklinearization_tests/` | Hardware/circle recordings and seven-group plot utility |
| `common_pressure_testing/` | Pressure protocols, trim/feedback sources, recordings and offline analyses, including `analyzeCommonPressureCircle.m` for compact reports and `plotCommonPressureTracking.m` for six logged-run figures |
| `ndi_orientation_recordings/` | Calibration raw poses and fitted references |
| `project_docs/` | Current state, history, decisions, evidence, issues and reusable theory |

MATLAB/Simulink and experiment-computer NI/SMC support are required for physical
runs; legacy serial is still used. Standalone offline analyses load MAT-files
and do not operate hardware. Animation uses `drawingArm_single_stationary_base.m`
and `HTM_nume.m`; no MEX is required. Existing local READMEs give execution
instructions and output paths; the older learned-predictor project is an external
research reference, not an installed dependency of this controller.
