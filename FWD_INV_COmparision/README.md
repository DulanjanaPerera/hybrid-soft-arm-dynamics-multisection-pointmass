# Offline stiffness sweep

From the repository root, run:

```matlab
addpath('FWD_INV_COmparision');
results = sweepStiffness;
```

The default input is `Comparision_results.mat`. Pass another full MAT-file
path as the first argument to analyze a compatible log. All outputs go into
`stiffness_sweep_results/`; rerunning overwrites those analysis outputs.

This replays logged commands through the existing MATLAB dynamics, without
loading Simulink models or connecting hardware. The sweep is 1350, 1250, ...,
250, 200 N/m. Only the scalar in K=k*[2 1;1 2] changes. Geometry, mass,
damping, gravity, nonlinear bound penalty, zero logged initial state, pressure
area, and 0.8 bar deadzone stay fixed.

The 1350 replay is checked against the saved model length trajectory. Commands
are interpolated on simulation time to match that baseline; physical poll
duration is reported separately. A physical-time replay or delay analysis
would be a different experiment. Invalid NDI samples are omitted from errors.

`sweep_scores.csv` lists per-length RMS and the reduced-coordinate score used
to select the best tested value. l1=-l2-l3 is dependent and is not counted as
an independent coordinate in that selection. `sweep_results.mat` preserves
all predictions, settings, input arrays, and scores. PNG plots show all
predictions, the best versus baseline, and error versus stiffness.

The selected value is conditional on this recording and all fixed assumptions;
it is not independent identification of physical stiffness. Nonlinear bounds
may dominate softer simulations. NearBoundFraction flags sampled l2/l3 within
2 mm of the +/-35 mm bounds. Calibration, actuator asymmetry, pressure response,
and timing can also affect the fit.

2026-10-02: Modeled length bounds raised to +/-35 mm. Existing input recordings
retain their original model settings; baseline replay differences may therefore
include the bound change. Results are regenerated under the new bounds.
