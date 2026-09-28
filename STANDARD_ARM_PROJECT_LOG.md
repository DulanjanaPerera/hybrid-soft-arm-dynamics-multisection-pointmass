# Standard arm project log

## Goal and scope

Build a three-section standard distributed-mass arm model for residual-dynamics learning and orientation-invariant gravity-compensated control, then extend it to a moving floating base. The generalized coordinates remain [l12;l13;l22;l23;l32;l33]. The nominal model is the existing distributed-section model, not the point-mass beta model. Keep section masses, geometry, Taylor order, translational-only soft-material kinetic energy, stiffness, damping, and actuation law unchanged.

## Baseline restored on 2026-09-24

Worktree: codex/standard-arm-floating-base, starting at 0e659b4 (codex/hybrid-simulation-only). Git status was clean before changes. The existing armS_standard_core.m, compact integral functions, geometry/Jacobian functions, and Christoffel function were retained unchanged. The following four source files were copied byte-for-byte from compare-pointmass-standard (their origin is commit 943c1957015064d75c0b4b095d7688e1f3c534dd):

- armS_standard_dynamics.m: six-coordinate RHS and the existing stiffness, damping, and actuation law.
- armS_standard_entry.m: fixed three-section runtime-parameter entry.
- build_armS_standard_mex.m: source MEX build and three MATLAB/MEX RHS comparisons.
- validate_armS_standard.m: independent 16-node Gauss quadrature plus symmetry, positive-definiteness, derivative, and Christoffel checks.

Their Git blob hashes match the source branch exactly: dynamics 0133aa669c2015c8624d446a7e133f6d9995f754; entry ed16aa75b0783be18121d6096f608942e4fe08e8; build 14edb67c1790dffb30b79c1b267147cbf6794d25; validation 34d75ad1db7508d629c8573fa85b778c9aacb4d9. No point-mass beta comparison files were imported. The existing hybrid model was not altered.

## MATLAB validation before MEX

MATLAB R2025a was used with N=3, L=0.278 m, r=0.013 m, mi=[0.1;0.1;0.1] kg, and g=[0;0;-9.81] m/s^2. The validator used dq=[0.01;-0.02;0.015;0.003;-0.005;0.008], a central-difference step of 1e-7, and these three q poses: zero; [-0.001;-0.001;-0.001;-0.001;-1e-6;-1e-6]; and [-0.008;0.003;-0.004;-0.006;0.002;-0.005].

| Pose | Relative M symmetry | Relative dM/dq error | Relative skew residual for Mdot-2C | Relative quadrature M error | Quadrature G error | Minimum eigenvalue of sym(M) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0 | 0 | 0 | 1.816e-15 | 0 | 2.160036370791e-01 |
| 2 | 1.893e-16 | 2.070e-10 | 1.534e-16 | 1.752e-15 | 1.574e-15 | 2.163611124478e-01 |
| 3 | 7.952e-17 | 5.542e-11 | 1.146e-16 | 1.112e-12 | 1.550e-15 | 2.328300925023e-01 |

All three Cholesky positive-definiteness checks and validator assertions passed. These are sampled checks, not a global positive-definiteness proof. A bent-pose MATLAB ode15s smoke run from pose 3 with zero initial velocity, D=600 I, K=2200 I, zero tau, mu=2000, and lKbounds=[-0.02;0.02;1e6] reached 0.002 s in 92 output points, all states finite, with endpoint state change norm 1.403147577790e-02.

## Source MEX build

The build was invoked after validation as build_armS_standard_mex(fullfile(tempdir,'armS_standard_mex_build')). MATLAB Coder and the selected Microsoft C++ compiler completed code generation and compilation in 179.5 s. The build script's scaled maximum MATLAB/MEX RHS errors across its three cases were 0, 5.731e-15, and 2.341e-15 (worst 5.731e-15), below its 1e-8 threshold. It copied the validated 380,928-byte armS_standard_mex.mexw64 into this worktree. This is a Windows binary; another OS needs a native rebuild. No build limitation occurred on this host. MATLAB initially could not start in the restricted execution sandbox (File system inconsistency); validation and build succeeded with normal MATLAB filesystem access.

## Stationary orientation-aware base proposal

Keep the current standard entry and core intact. Add a thin runtime interface accepting R_world_from_arm (3x3 proper rotation) and g_world (3x1), in addition to the existing standard entry inputs. With the base stationary and orientation fixed during a solve, compute g_arm = R_world_from_arm.' * g_world and pass g_arm into armS_standard_entry. A constant base translation has no effect on the generalized gravity vector. No base angular velocity, angular acceleration, or translational acceleration terms belong in this stationary interface.

The current convention is G = integral(J_arm.' * g_arm dm), with g=[0;0;-9.81] in the reference arm frame, and the RHS subtracts G: M*qdd = tau - (C+D)*dq - G - K(q)*q. Preserve this convention and its sign; gravity compensation in this RHS uses tau_grav=G. Do not substitute the point-mass beta model.

Checks already made: R_world_from_arm=I reproduced the original full RHS exactly (maximum absolute difference 0). For a +90 degree roll about arm x, g_arm=[0;-9.81;0]; M and C differences were exactly 0 at pose 3, while the gravity-vector difference norm was 6.928982306785e+01. The existing three-pose validation also passed with this rotated gravity; quadrature G errors were 1.503e-15, 1.697e-15, and 1.426e-15.

Before moving-base dynamics, test the new interface at identity, yaw about world gravity (same G), +/-90 degree pitch and roll, and sampled proper rotations. For each, compare G with independent world-frame Jacobian quadrature using R_world_from_arm*J_arm, verify M and C are orientation-independent, verify the existing Mdot-2C identity and dM/dq, and compare MATLAB/MEX RHS with runtime orientations. Test gravity linearity and the existing sign convention in static compensation. Require orthonormal R and det(R)>0 at the interface. Then derive and validate moving-base velocity, acceleration, and coupling terms separately.

## Decision and next step

Use the standard distributed-mass baseline as nominal dynamics. Retain the inherited functions, physical parameters, and six-coordinate convention. Implement and test the stationary orientation wrapper next, including runtime MATLAB/MEX parity; proceed to floating-base dynamics only after orientation tests pass. Nothing was pushed, and neither Ozi nor compare-pointmass-standard was edited.

## Stationary orientation interface implemented on 2026-09-24

The checkout was confirmed as the `standard-arm-floating-base` linked Git worktree on `codex/standard-arm-floating-base` at `0e659b4`; `origin` is `https://github.com/DulanjanaPerera/hybrid-soft-arm-dynamics-multisection-pointmass.git`. MATLAB's Source Control menu did not identify the linked worktree's `.git` pointer file, but Git itself reported the worktree, branch, remote, and untracked baseline correctly. No Git initialization was needed. The four original untracked MATLAB source files retain the exact blob hashes listed above, and the original standard MEX was not replaced.

Added `armS_stationary_base_entry.m`: it accepts runtime `R_world_from_arm` and `g_world`, requires a finite 3x3 proper orthonormal rotation and finite 3x1 gravity, computes `g_arm=R_world_from_arm.'*g_world`, and calls the unchanged `armS_standard_entry`. This is a stationary base with orientation fixed during a solve. Its six coordinates, physical parameters, mass model, stiffness, damping, actuation, and gravity sign are unchanged. Added `stationaryBaseTestRotations.m`, `validate_armS_stationary_base.m`, and `build_armS_stationary_base_mex.m`. The new 385,024-byte Windows `armS_stationary_base_mex.mexw64` is separate from the original MEX.

MATLAB R2025a validation used the three poses and parameters above, with identity, yaw 60 degrees about world gravity, signed 90-degree roll and pitch, and two deterministic compound rotations (24 pose-orientation cases). The independent reference uses 16-node Gauss quadrature of `J_world=R_world_from_arm*J_arm`. Maximum results: identity RHS absolute error 0; relative world-quadrature M error 1.112e-12; world-quadrature G error 1.776e-15; orientation changes in M/C/dM exactly 0; relative dM finite-difference error 2.070e-10; relative skew residual for Mdot-2C 1.534e-16; gravity-linearity residual 8.335e-15; and static-compensation RHS norm 2.453e-15 using `tau=G+K(q)q`. Yaw preserved gravity in the arm frame. Reflection and nonorthonormal inputs were rejected. All assertions passed.

MATLAB Coder and Microsoft C++ built the separate MEX in 305.4 s. Its build script compared all eight orientations at runtime and reported a worst scaled MATLAB/MEX RHS error of 2.994e-15. The final validator compared 24 pose-orientation cases plus one changed runtime gravity vector and reported a worst scaled error of 9.728e-15 across 25 checks. MATLAB startup again failed only in the restricted filesystem sandbox; validation and build passed with normal MATLAB filesystem access. Nothing was committed or pushed; neither Ozi nor compare-pointmass-standard was edited.

Next: derive and validate moving-base velocity, acceleration, and coupling terms in a separate step. This stationary wrapper does not model moving-base dynamics.


## Stationary-orientation simulation runner on 2026-09-24

Added `runDynamicSimulation_armS_stationary_base.m` in the independent
OneDrive clone under `Lab/Research/float_dynamic/standard-arm-floating-base`.
It runs one fixed-orientation `ode15s` solve per editable
`[roll pitch yaw]` row, using `Rz*Ry*Rx`, the original standard distributed
arm RHS, and either the separate stationary-base MEX or MATLAB source.
The script stores each trajectory with its rotation and arm-frame gravity,
and plots coordinate histories and final world-frame backbone shapes.
It does not add point masses or moving-base inertial terms.

MATLAB R2025a smoke run used identity, +90 degree roll, and +90 degree
pitch with the script's default physical parameters and bent initial pose.
Each solve reached 0.1 s with 11 requested output points and all states
finite. The final state-change norms were 3.669358e-02, 3.261798e-02,
and 4.183200e-02, respectively. Plotting completed without error. The
new runner was not committed or pushed.

## Fixed-orientation standard-arm animation on 2026-09-24

Added `runDynamicSimulation_armS_stationary_base_custom.m` and
`drawingArms_stationary_base.m`. The editable runner accepts one fixed
base orientation as roll/pitch/yaw degrees or a measured 3x3
`R_world_from_arm` matrix. It uses the unchanged standard distributed
arm dynamics and existing pressure, stiffness, and damping laws. The
animator rotates the entire recursive backbone into world coordinates,
draws the arm-frame axes at the fixed base origin, and can record MP4
plus MAT output. Recordings are ignored by Git. The gravity default is
`[0;0;9.81]`, matching the user's current orientation-sweep runner;
the prior validation used `[0;0;-9.81]`. Gravity direction remains an
editable world-frame input.

MATLAB R2025a ran the default identity case for 5 s with 301 frames,
finite states, and final state-change norm 1.318550e-02. A temporary
matrix-input +90 degree roll test ran for 0.1 s with 3 frames; the
rendered section-1 endpoint matched `R_world_from_arm*p_tip` with
error 0. The MP4 and MAT were created at 56,446 and 1,937 bytes.

## Animation figure reuse and fixed limits on 2026-09-24

Updated `drawingArms_stationary_base.m` to clear and reuse Figures 1
and 2, following the original `drawingArms.m`. The arm axes now use
fixed X and Y limits [-1,1] m, Z limits [-1.5,1.1] m, and view
[11,13]. The world-frame backbone rotation remains unchanged.
A MATLAB R2025a test called the animator twice in one session
(identity then +90 degree roll): both figure handles were reused,
exactly two figures existed, and all three limits matched the
requested values. No dynamics files changed.

## Single-module branch and implementation on 2026-09-24

Created local branch `codex/single-module` in the independent OneDrive clone at `Lab/Research/float_dynamic/standard-arm-floating-base` from `99f281c`. Preserved the preexisting modified and untracked files; did not edit the three-section standard or stationary-base source. No commit or push was made.

Added `armS_single_core.m`, `armS_single_dynamics.m`, and `armS_single_entry.m` for one distributed section with `q=[l12;l13]` and `X=[q;dq]`. The core uses the existing compact integrated-position and integrated-Jacobian products with `N=1`. Mass, Christoffel, gravity, nonlinear stiffness, damping, and actuation laws follow the standard distributed-mass model. The thin `armS_single_stationary_base_entry.m` accepts runtime `R_world_from_arm` and `g_world`, validates a proper rotation, computes `g_arm=R_world_from_arm.'*g_world`, and calls the single-module entry. It introduces no moving-base terms.

Added `validate_armS_single.m`, `build_armS_single_stationary_base_mex.m`, `runDynamicSimulation_armS_single_custom.m`, and `drawingArm_single_stationary_base.m`. The runner accepts one fixed orientation as RPY degrees or a measured rotation matrix, an editable signed two-pressure vector, and source or MEX RHS. The animator transforms the one-section backbone into world coordinates, reuses Figures 1 and 2, and keeps X/Y [-1,1] m, Z [-1.5,1.1] m, view [11,13]. Optional recordings go into the Git-ignored `single_module_recordings/` folder.

MATLAB R2025a validation covered three poses and eight fixed orientations. All assertions passed. Maximum relative independent 16-node quadrature errors: M 2.305e-15 and G 2.359e-15. The isolated first section of the three-section model matched to 4.441e-16; the identity-orientation RHS matched exactly (absolute error 0). Maximum relative finite-difference dM error was 3.270e-10; the Mdot-2C skew residual was 1.145e-16; static-compensation RHS norm was 4.512e-16.

MATLAB Coder and the configured Microsoft C++ compiler built the separate `armS_single_stationary_base_mex.mexw64` in 82.0 s. It passed 24 pose-orientation MATLAB/MEX comparisons, including runtime mass, torque, and gravity changes, with worst scaled RHS error 0.000e+00. The default MEX runner reached 5.000 s with 301 frames, final state-change norm 1.414213e-02, and final q approximately [-6.87667e-10;-6.87667e-10]. A 0.1 s signed-pressure check gave pressure [-2;-2] bar, q=[-0.00360418;-0.00360418], tip X=-0.075119314 m; pressure [2;2] bar gave q=[0.00360418;0.00360418], tip X=0.075119314 m. A +90-degree roll animation check gave world-frame tip-geometry error 0, reused both figure handles, produced exactly two figures, and retained all requested axis limits.

This branch provides one stationary-orientation module for control experiments. Recursive multi-section and quadruped extensions remain future work; moving-base dynamics were not introduced.

## Git commit record for the single-module branch

The earlier standard-arm validation, Maple worksheet, MEX binaries, and stationary-orientation simulation and animation files were committed as `a020747` (`Add standard-arm validation and orientation simulations`). The single-module source, MEX, runner, animation, validation, README, and recording ignore rule are committed in the following commit (`Add single distributed module with stationary orientation`). The earlier statements that files were uncommitted describe their state when those entries were written. Both commits are local to `codex/single-module`; no push was made.

## Two-computer handoff and experimental direction on 2026-09-28

The active repository is the independent clone at `Lab/Research/float_dynamic/standard-arm-floating-base`, on `codex/single-module`. Its `origin` is `https://github.com/DulanjanaPerera/hybrid-soft-arm-dynamics-multisection-pointmass.git`. The earlier statements about uncommitted or unpushed work above are historical snapshots. The commits after the previously published `408f5bc` are:

| Commit | Change |
| --- | --- |
| `e81b6e1` | Add the two XYZ-based single-module inverse-kinematics functions. |
| `4c8d5c0` | Improve straight-pose handling, bounded bend estimation, and fit diagnostics. |
| `ad8f982` | Add the reviewed two-sensor Aurora acquisition driver, port-handle class, and MATLAB sampling script. |
| `9e0f684` | Add the Simulink NDI reader and model; connect the reader to `AuroraDriver_2`. |
| `675b503` | Retain the original `AuroraDriver.m` as a labeled legacy reference. |

The control experiment is for one stationary-base module with `q=[l12;l13]`. Base orientation may be measured and held fixed for gravity compensation; the model does not predict base motion. The intended sequence is (1) compare measured trajectories with the nominal distributed-mass model, (2) insert the already identified coupled-stiffness matrix and measure the change in prediction accuracy, then (3) learn remaining hysteresis, direction-dependent stiffness, and damping effects. These are planned experiments, not completed validations. A residual model and controller, possibly MPC, have not been implemented. The later quadruped will reuse four modules and the recursive dynamics structure. Do not interpret the stationary-orientation wrapper as moving-base dynamics.

The planned control sensor is the IMU; NDI XYZ is the reference for geometry reconstruction, system identification, and validation. Collect synchronized IMU, NDI, pressure, and timestamps when possible so that estimates used during training can be reproduced by the eventual IMU-based controller. The UKF for length-change estimation exists in another folder but has not been imported or integrated here. No residual-learning or coupled-stiffness code was added in the NDI commits.

`f20260219_2_task2config_withLenExt.m` fits bend direction `theta` and bend angle `phi` from the tip-sensor XYZ relative to the base sensor, expressed in the arm base frame and in meters. It accounts for base-sensor offset, module length, tip-sensor offset, and the tip's backbone fraction. It reports fit residual and diagnostics. Near straight poses, bending direction is unobservable and the function uses a caller-supplied prior direction. `f20260219_2_task2length_withLenExt.m` maps the fitted configuration to three PMA length changes; use its entries `l(2:3)` for the two dynamic coordinates. Verify sensor offsets and frame alignment on the physical rig before using these estimates as ground truth.

`AuroraDriver_2.m` and `PortHandle.m` are the active acquisition path. One BX poll supplies both sensors per sample. The driver rejects incomplete scalar reads, does not map an unknown handle onto sensor 1, and returns NaN positions for missing or unusable readings; validity and frame numbers are available as extra outputs. `m20251106_3_double_sensor_reading.m` records both sensors, host timestamps, frame numbers, and validity into `ndiSamples`, preserving raw positions in millimeters and a converted meter copy. Its defaults are COM12, 921600 baud, and 500 samples. Sensor order follows discovered handle order; the hardware must establish which ID is base and which is tip. Neither active nor legacy driver verifies BX CRCs. Both still use MATLAB's legacy `serial` interface, which R2025a warns will be removed in a future release.

`NDI_readingOnly.slx` calls `AuroraReader_s2_Position_only.m`, whose defaults are COM12 and 230400 baud and whose MATLAB System block uses interpreted execution. The model assumes sensor 1 is the base and sensor 2 is the tip. It forms sensor-2 minus sensor-1, applies the fixed rotation `[1 0 0; 0 0 1; 0 -1 0]`, converts millimeters to meters, and calls the XYZ-to-`theta`/`phi` function with `L=[0.05;0.178;0.05]` m and `xi=1`. The model's fixed step is 0.01 s. Its fixed frame rotation, sensor order, offsets, achievable sample timing, and both baud-rate settings need experimental verification. The Simulink model carries the previous `phi` in memory but does not feed a previous `theta` to the inverse-kinematics function, so straight-pose direction currently falls back to that function's default.

Offline MATLAB R2025a checks passed for the new classes, including missing-sensor invalidation, the reader's missing-second-sensor error path, and model diagram compilation. A mock two-sensor Simulink run produced position `[0.0100;0;0.2779]` m, `theta=0`, and `phi=0.0719` rad at the displayed precision. The standalone inverse-kinematics function returned `theta=0`, `phi=0.0719101904` rad for the same position. No physical NDI acquisition was performed on this computer: `serialportlist('all')` reported only COM1, while the NDI hardware is connected as COM12 to the separate experiment computer. `AuroraDriver.m` is an older class retained for reference and is not called by the Simulink reader; it can retain stale data and does not reject incomplete BX reads.

For collaboration, create a **separate clone of this same Git remote** on the experiment computer and a Codex/MATLAB project pointing to that clone. Check out `codex/single-module`; do not initialize a second independent repository or work from one shared OneDrive-synced checkout on both machines. Use the remote to exchange commits. The experiment computer should perform the COM12 hardware checks and commit their outcomes; the modeling computer can continue dynamics and analysis. If both computers edit concurrently, use separate branches and merge reviewed changes. Keep raw experiment recordings in an appropriate data location and commit code, configuration, concise results, and provenance to the repository.

Next on the experiment computer: confirm COM12 visibility and port ownership; start with a short two-sensor acquisition; record handle IDs, valid-frame counts, actual sample timing, and known-pose geometry; confirm base/tip assignment, sign and axis mapping, millimeter-to-meter conversion, and the measured offsets; then run the Simulink model briefly and compare its position and angles with the standalone functions. Hardware findings should be appended to this log. Only after the measurement chain is trustworthy should the nominal-versus-coupled-stiffness comparison and later residual identification proceed.
