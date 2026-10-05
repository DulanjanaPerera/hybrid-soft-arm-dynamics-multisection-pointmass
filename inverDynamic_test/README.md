# Inverse dynamics test analysis

Run in MATLAB from the repository root:

```matlab
addpath('inverDynamic_test');
results = analyzeInverseDynamics;
% Alternatively:
% results = analyzeInverseDynamics(fullfile(pwd,'inverDynamic_test','another_run.mat'));
```

The function reads saved `out.logsout` data only. It never loads or simulates
models and never opens NDI or NI hardware. It uses the existing repository
math functions for nominal static forces without modifying them. Outputs are
written to `inverDynamic_test/results/` and overwritten on a subsequent run:

- `hold_summary.csv`: one row per constant sampled-length target lasting at
  least five physical seconds.
- `analysis_results.mat`: settings, model parameters, summary, and aligned data.
- `tracking.png`, `errors_and_forces.png`: complete-run plots.
- `report.txt`: timing, validity, and interpretation limits.

Hold numbers refer to all reference segments, including short segments during
slider changes, so the retained hold numbers are not consecutive. The last
three physical seconds are evaluated. Stationarity requires maximum filtered
speed of 0.20 mm/s and maximum coordinate range of 0.50 mm. These are explicit
engineering thresholds, not identified noise limits. Change the `cfg` block
to assess sensitivity. The function requires logged desired dl and ddl to be
zero. Slowly changing targets and holds shorter than five seconds are plotted
but excluded from static tables.

Length errors are desired minus measured. Settling requires both coordinate
errors to stay within 1 mm until the end of the hold for at least two seconds.
A settling time of zero means the initial pose already met that tolerance;
it does not establish a transient response. Stationary tails with large
tracking errors have not settled to the requested target. Theta errors use
circular statistics and exclude unobservable or nearly straight poses.

Static required force is calculated at the mean measured tail pose, including
the nonlinear bound penalty. Its residual is required minus logged available
force. Use it only for stationary tails. Available force is the inverse's
estimate after pressure scaling, not measured force. The independent pressure
mapping check applies the same provisional 0.8 bar deadzone, with a separate
0.0065 m pressure radius; the geometry offset remains 0.013 m. A tiny mapping
or desired-force check verifies calculation consistency, not plant accuracy.

Actual chamber pressure and post-Kill voltage/pressure are not assumed to be
logged. Static residuals cannot independently identify stiffness, deadzone,
pressure response, calibration, or channel mapping. Transient plots do not
validate mass, Coriolis, or damping with these zero-rate inverse inputs.
NDI lengths are inferred from orientation, not independent muscle strain.
All plots use elapsed physical BX poll time; simulation time is retained in
the table for locating events in Simulink. Parameter logs are held at their
last recorded value. Invalid samples are excluded from tail metrics.
