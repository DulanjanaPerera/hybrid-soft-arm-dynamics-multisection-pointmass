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

