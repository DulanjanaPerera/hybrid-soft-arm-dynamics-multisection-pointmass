# Hybrid arm simulation

This branch contains the MATLAB sources needed to simulate a three-section
distributed-mass arm with independently attached backbone point masses. It
contains no beta-fitting, point-mass replacement comparison, or study outputs.

## Dynamics and parameters

The six coordinates are `[l12;l13;l22;l23;l32;l33]`. The model adds each
attached mass contribution to the standard distributed section mass:

```text
M = M_distributed + M_added
C = C_distributed + C_added
G = G_distributed + G_added
```

The model keeps translational kinetic energy, the original gravity convention,
stiffness law, damping, and pressure input. The added point masses use ordinary
translational kinetic energy (unit beta). In the editable runner,
`params.mi` gives the existing distributed masses in kg, while
`params.addedMass` gives additional masses in kg for sections 1–3.
`params.addedXi` gives their section backbone positions from 0 (base) to
1 (tip). All three masses and positions stay fixed during each ODE solve.

## Run

Open MATLAB in this folder and run:

```matlab
runDynamicSimulation_armS_hybrid_custom
```

Edit the configuration at the top of that script for physical parameters,
initial coordinates and velocities, duration, solver settings, and recording.
`useMex=true` uses a MEX built for the current computer; `useMex=false`
calls the MATLAB dynamics without a MEX. The script leaves `t`, `X`, `X0`, and `params` in the
workspace. When `recordVideo=true`, it saves the arm animation as an MP4
and the trajectory/parameters as a matching MAT file in
`hybrid_recordings/`. That output directory is ignored by Git. The
coordinate-history figure remains interactive and is not included in the MP4.

## Build on macOS or another computer

The included `armS_hybrid_mex.mexw64` was built for MATLAB R2025a on
Windows and does not run on macOS. The complete MATLAB source dependency
chain and `armS_hybrid_entry.m` are included. On the other computer, install
MATLAB Coder and a C++ compiler supported by that MATLAB release, then in
MATLAB select the compiler and build:

```matlab
mex -setup C++
build_armS_hybrid_mex
```

The build runs in a folder under `tempdir`, checks two runtime input cases
against the MATLAB entry function, and copies the validated native MEX next
to the scripts. Its extension comes from `mexext`, so the build script does
not hard-code a Windows or macOS binary name. If MATLAB Coder is unavailable,
set `useMex=false` in the editable runner to use the MATLAB implementation.
The native MEX build has been verified on Windows; the macOS build must be
run and checked on the recipient's Mac.

You can also choose another temporary build folder:

```matlab
build_armS_hybrid_mex(fullfile(tempdir,'my_hybrid_build'))
```

The build uses the same runtime interface; changing added masses or positions
does not require recompilation. `armS_hybrid_entry.m` is the MEX entry point.
The point-mass source core retains an internal beta argument because it is
shared source, but the hybrid model always supplies `ones(3,3)`; there are
no beta values to tune in this simulation branch.

The attached mass is modeled as a point on the section backbone. An offset
mount or a rigid payload with meaningful rotational inertia needs an
extended model.

## Standard distributed arm at stationary base orientations

Run `runDynamicSimulation_armS_stationary_base` from this folder to compare
fixed base orientations without attached point masses. Edit
`orientationRPYDeg` near the top of the script; each row is
`[roll pitch yaw]` in degrees, with
`R_world_from_arm = Rz(yaw)*Ry(pitch)*Rx(roll)`. The orientation is fixed
during each ODE solve. The script uses the same six coordinates and existing
standard distributed-mass dynamics, stiffness, damping, and pressure law.
It leaves `results`, `X0`, and `params` in the MATLAB workspace; each
`results(k)` contains the rotation, arm-frame gravity, time vector, state
history, and elapsed solve time. It plots coordinate histories and final
arm shapes in world coordinates when `showPlots=true`.

Set `useMex=false` to run MATLAB source. With `useMex=true`, the separate
`armS_stationary_base_mex` binary is required. Rebuild it on another
platform with `build_armS_stationary_base_mex`. This stationary-base model
accounts for orientation-dependent gravity; it assumes zero base velocity
and acceleration throughout each solve.

## Animate one fixed base orientation

Run `runDynamicSimulation_armS_stationary_base_custom` to simulate and
animate one stationary base orientation, using the standard distributed arm
without attached point masses. Edit `baseRPYDeg` with
`baseOrientationMode='rpy'`, or set `baseOrientationMode='matrix'` and
provide a measured `R_world_from_arm_input`. The matrix maps arm-frame
axes into world coordinates and is held fixed throughout the solve.
The script leaves `t`, `X`, `X0`, and `params` in the workspace.

The animation draws the backbone and base axes in world coordinates and
shows the six coordinate histories. Repeated runs reuse Figures 1 and 2.
The arm view uses the original fixed limits: X and Y [-1,1] m,
Z [-1.5,1.1] m, and view [11,13]. Set `showAnimation=false` for a
solve without figures. Set `recordVideo=true` to save an MP4 and MAT file
under `standard_orientation_recordings/`; this folder is ignored by Git.
Set `useMex=false` if the native stationary-base MEX is unavailable.
The editable gravity default follows the current orientation-sweep runner.
