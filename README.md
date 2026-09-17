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
