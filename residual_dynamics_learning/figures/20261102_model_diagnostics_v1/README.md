# Retained model-diagnostic figures

Created 2026-11-02 from existing saved numerical results; no new model fitting
or hardware execution. See project_docs/FIGURE_INDEX.md for archive policy.
F01-F04 are saved as editable MATLAB FIG, vector PDF and 300 dpi PNG. F05 is
the original existing PNG, copied without alteration. source/ contains the
plotting function, source reports and the original example plotting code.
data/ includes original CSVs and plotted-value tables. source_manifest.csv
links copied inputs to their original paths; archive_hashes.csv fingerprints
the bundle excluding the hash manifest itself.

## Captions and limitations

**F01_prediction_horizon.** Mean vector path RMSE of nonlinear model predictions
over 10,20,40 sensor steps in 30 held-out October19 windows, six per waveform.
The x-axis uses the mean recorded physical duration at each step count. Errors
exclude the initial frame. Initial q/dq are observed/causally derived; subsequent
states are integrated without measured-state resets, using recorded commands.
Mechanics improves over the original, but persistence/constant velocity remain
better. Windows can overlap and are not independent repeats. Selected memory
uses validation-selected delay0.2s; it is not an identified physical delay.

**F02_static_validation.** Vector equilibrium-position error on each of four
held-out near-static pressure holds. Original and calibrated predictions are
paired by run/trial; lines connect hold categories and imply no interpolation.
Calibrated stiffness/bias used15 StepWave training holds. All static-equilibrium
solves start from zero. Small sample; no error bars or population-level claims.

**F03_damping_diagnostic.** Training output-error cost for shared diagonal
damping40,100,250,500,1000 Ns/m with calibrated stiffness/bias fixed and no lag.
Ten training windows; cost is the sum over windows of mean squared coordinate
errors in mm, i.e. mm^2, not mean vector RMSE. Dashed line marks the actual fit
upperbound250. Values above it are a diagnostic profile, not newly fitted
parameters or validation-selected bounds. Lower cost at high D can reflect
suppression of inaccurate motion; do not call it identified physical damping.

**F04_forward_screening.** Mean per-run vector next-position RMSE on October19
validation and October27 external data, with a logarithmic y-axis. Outputs are
one native sensor step, not uniformly50ms physical time. Frozen local nominal
transition differs from the fully nonlinear model in F01. The pressure-run
October27 subsets are sparse (StepWave46, combined26 anchors); cross-day timing
differs, so equal-weight macro averages have limited representativeness.
Kinematic ridge is a separately frozen ablation. Keep the source screening
report's leakage/shortcut, finite-budget and external-evaluation limitations.
This figure does not establish a learned physical residual or control benefit.

**F05_example_rollouts_original.** Preserved original nonlinear rollout plot:
first selected validation window per waveform, observed lengths, original model,
calibrated+selected memory and constant velocity. Raw plot is illustrative;
not all validation windows. Exact inputs/initial states are stored in the
calibration identification_data.mat and original source is retained. This
legacy PNG has no original FIG/vector export; do not claim otherwise.

## Reproduce in a new folder

From the repository root in MATLAB:

```matlab
addpath('residual_dynamics_learning');
archiveResidualFigures(pwd,fullfile(pwd,'residual_dynamics_learning', ...
    'figures','new_version'));
```

The archived source snapshot reproduces F01-F04 from canonical saved results;
copied data preserve inputs if canonical results are later superseded. No
simulation or Optimization Toolbox is required for this export function.
It refuses destinations containing prior outputs (a prewritten README is allowed).
The archive README/index/hash manifest are
curated records, separate from the plotting function's generated outputs.
