# P1 stepped stiffness experiment

`planar_stiffness_estimation.slx` drives P1 through 0, 0.5, 1, 1.5, 2, 2.5,
3, 2.5, 2, 1.5, 1, 0.5, 0 bar, each held for 10 simulation seconds.
It repeats three times, then commands zero. Stop time is 400 s: 390 s of
cycles plus a final zero hold. P2/P3 stay zero. The standalone
`stepped_ramp_wave.m` matches the MATLAB Function3 block for inspection.

Existing valve Kill/dashboard behavior is preserved; verify its live state.
The model uses NDI startup preflight. This setup has not been run on hardware.
Explicit logging includes commanded and post-Kill/saturation valve pressures,
NDI lengths, orientation angles, validity, and physical poll timestamps.
Initialize ndiConfig and save the run as a MAT-file containing `out`.

After the physical experiment, run from the project root:

```matlab
addpath('planar_stiffness_test');
results = estimatePlanarStiffness('FULL_PATH_TO_SAVED_RESULTS.mat');
```

The analysis searches k=100:10:10000 N/m in K=k*[2 1;1 2]. At each pressure
hold it matches nominal static force to force inferred from commanded pressure
at the mean measured pose. This is an equilibrium force-residual search,
not a time-domain trajectory fit. Gravity and the +/-35 mm bound penalty are
retained. Effective pressure radius is 6.5 mm; geometry radius is 13 mm.
Deadzone is set to zero for this first simple experiment, per the requested
assumption. Thus estimated stiffness includes pressure/deadzone/model effects.

Last three physical seconds of holds >=5 seconds are used. Low displacement
(<0.5 mm norm), zero pressure, invalid tails, and tails with >0.5 mm coordinate
range are excluded from means. These thresholds are editable in cfg. A range
check is a simple stationarity screen, not proof of equilibrium.
Loading/unloading and repeated cycles are pooled by pressure level. Outputs:
`results/hold_estimates.csv`, `results/stiffness_variation.csv`,
`results/stiffness_variation.mat`, and `results/stiffness_vs_phi_expanded.png`.
The table includes pressure, mean phi and lengths, mean/std stiffness, count,
force residual, and grid-edge flags. Inspect residuals and grid-edge flags
before using the table in a controller. NDI lengths are orientation-derived;
the fitted scalar does not independently identify physical material stiffness.
No controller automatically uses the output table.

The estimator uses `des_pressure`, or separate `P1`/`P2`/`P3` logs when
the combined signal is absent (input commands, not pressure sensors).
This assumes the Kill gate passed the input during the experiment; otherwise
commanded pressure does not represent what was requested from the valves.
`valvePressure_bar` is an optional post-Kill software signal, not measured
chamber pressure, and is no longer required by the estimator.

## Stiffness from phi

Run `fit = fitStiffnessCurve` after estimation, then
`[k, inRange] = stiffnessFromPhi(phi_rad, fit)` to evaluate stiffness in N/m.
The curve uses shape-preserving cubic interpolation (PCHIP), excluding pressure
levels for which any included hold reached the stiffness search boundary.
The usable range is saved in fit.PhiRange_rad; the evaluator returns
NaN outside it rather than extrapolating unverified stiffness. The knots
are interpolation data, not identification of a general nonlinear law.
Standard deviations across holds are shown as variability, not confidence bands.

Outputs: `stiffness_phi_fit.mat` (piecewise-polynomial coefficients and source
knots), `stiffness_phi_knots.csv`, `stiffness_phi_curve.csv`, and
`stiffness_phi_fit_expanded.png`. This P1 curve is conditional on the zero-deadzone
pressure assumption and pooled loading/unloading measurements. It is not
installed in a Simulink controller. Refit after regenerating the source table.

Expanded search update: all six pressure levels are now included in the PCHIP
curve, spanning approximately 0.09425 to 2.05594 rad. The evaluator still
returns NaN outside the fitted interval, including phi=0. Choose an explicit
out-of-range policy before a forward simulation. Large low-phi variation is
retained in the saved standard deviations. No forward model is changed by fitting.

## MATLAB Function block version

Paste `stiffnessFromPhiBlock.m` into a MATLAB Function block. It has a scalar
phi input (rad) and scalar k output (N/m), hardcoded PCHIP coefficients, and
no file loads. Use K=k*[2 1;1 2], retaining the separate bound penalty.
Below/above the measured range it holds the endpoint values 4132/953.333 N/m.
Nonfinite phi uses 4132; this fallback does not replace NDI validity gating.
The standalone `stiffnessFromPhi` evaluator still returns NaN outside its
measured range; endpoint extension is explicit in this block version only.

For a purely continuous forward model, phi computed from the Integrator's
current length state can feed this block directly; the state integrator
normally breaks direct feedthrough. If deliberately using previous phi, use
an explicit Memory block (previous major step) or a discrete Unit Delay with
a defined sample period and initial phi=0. These introduce different delays.
No delay, controller change, or wiring is installed automatically. This P1
empirical curve is not validated for other bending planes or hardware control.
Refitting the curve does not automatically update the hardcoded block.
