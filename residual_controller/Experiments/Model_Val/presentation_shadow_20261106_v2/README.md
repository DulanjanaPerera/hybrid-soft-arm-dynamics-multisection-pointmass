# Real-time residual-model motion and error presentation figures

Created2026-11-06. User explicitly selected the **later real-time shadow
recordings**, not the earlier November2 GRU64 free-run examples. These are
matched previously issued0.05s predictions from nominal and residual+nominal
models, compared with current valid NDI-derived measurements. No fitting,
interpolation, smoothing, model/SLX/hardware run, controller edits, commit or
push was performed. Raw recordings remain unchanged.

## Recommended presentation figures

Each stem has PNG4320x2430/300dpi, vectorPDF/SVG and editableMATLAB FIG.

| Experiment/view | Theta/phi trajectories + angular errors | Length trajectories + length errors |
|---|---|---|
| Complete circle, requestedphi1rad | circle_00g_phi1_complete_cycle_angles | circle_00g_phi1_complete_cycle_lengths |
| First helix loading/unloading episode | helix_00g_first_episode_angles | helix_00g_first_episode_lengths |
| Full sine experiment, no payload | sinewave_00g_full_angles | sinewave_00g_full_lengths |
| Full sine experiment,100g | sinewave_100g_full_angles | sinewave_100g_full_lengths |

Also retained: circle_00g_full_angles/lengths and helix_00g_full_angles/lengths,
covering the entire375s sweep and290s helix recording. There are12figures total,
not12new experiments. The two sine files each cover212s, all3pressure channels
and both repeats. No pressure-channel segments are stitched together.

Every figure has trajectories in the top row and corresponding errors below.
Black=measured, orange dashed=nominal, blue=residual+nominal. Only requested
predictors appear; the constant-velocity baseline is not plotted. Earlier
audits show it remains important for judging short-horizon accuracy.

## Source and window identity

| Stem | Raw source in parent folder |
|---|---|
| circle_00g | circlesweep_test.mat |
| helix_00g | helix_FWD_INV_test.mat |
| sinewave_00g | sinewave_test.mat |
| sinewave_100g | sinewave_test_100g.mat |

The circle close-up is159–169s, the first complete10s phase5cycle at requested
phi_reference1rad (metadata selection, not smallest error). Actual measured
bend is substantially below1rad. This is prediction of the recorded response,
not successful tracking of the requested circle. The helix close-up is5–65s,
first chronological trial including both phase6loading and7unloading. Full
records remain available. Payload labels follow user/filename; attachment and
calibration timing are not embedded, so no isolated causal payload claim.
All four model versions are residual_FWD_INV_validation1.51.

Close-up horizontal axes start at the selected source-window start; full
records start at0simulation seconds. Host and simulation time are kept distinct.
Physical0.05s horizon comes from the logged context; prospective scoring uses
the existing logged timing tolerance. Retimed/future-issued outputs are not
substituted into these comparisons.

## Exact prediction/measurement semantics

Nominal and augmented trajectories use **predictionAtMeasurement_m** columns2
and1, respectively. Both are measured-state-initialized one-step predictions
issued before the target observation, with identical score masks. They are
not a free nominal simulator versus a measurement-reset residual model. No
future target is used to construct the earlier pressure-conditioned prediction.
The plot's full duration does not make these full-duration free-run predictions.

The observed q=[l2;l3] is taken from meas_length and checked against the shadow
diagnostic qMeasured_m. It denotes inferred length changes, not independent
strain measurements. No l1 panel is added because l1=-l2-l3 in this model.
Angles are reconstructed consistently from each model/measurement's q:

```
u=(l2+l3)/0.013
v=(l3-l2)/(sqrt(3)*0.013)
theta=atan2(v,u)
phi=hypot(u,v)
```

Measured geometric angles reproduce the independently logged thetaOrientation/
phiOrientation within9e-16rad (theta modulo2pi on observable rows). The separate
XYZ-fit theta/phi signals and the separate old nominal simulator are not mixed
into the comparison. No reconstructed XYZ is shown in these figures.

All errors are **prediction minus measured**, not desired-tracking error.
Theta errors use atan2(sin(delta),cos(delta)); theta trajectories are wrapped
[-pi,pi], with gaps inserted at wrap discontinuities to avoid false vertical
lines. Theta error axes also show[-pi,pi]. Near-straight predicted direction
is undefined and masked when predictedphi<=0.5degree. The measured theta uses
its own logged thetaObservable flag. This threshold is a declared geometric
display/scoring convention, not calibrated model uncertainty.

Nominal and augmented theta MAE in prediction_metrics.csv use the **same
intersection of available angular rows**, because one model can predict a
near-straight configuration at a different time. The plotted theta trajectories
and errors retain each model's own legitimate availability; undefined theta
is not assigned zero error. Length/phi metrics use the common length-score mask.

## Validity and independent verification

Physical measurement validity requires orientationValid==1,basevalid==1 and
finite inferred lengths. Model curves additionally require scoreValid and
finite matched predictions. Measured curves remain visible through model
warmup when the NDI measurement is valid; unavailable model curves remain NaN.
Sensor gaps, observer resets and warmup are not bridged. Native signals use
bounded previous-sample alignment(maxage.075s) if their clocks differ, with no
extrapolation. Duplicate thetaObservable names are selected explicitly from
the NDI MATLAB System source. Full signal inventories retain paths.

For every scored row in allfour recordings, the matched prediction equals a
strictly earlier issued prediction EXACTLY. Independent length error
recalculation agrees with the logged diagnostic exactly, and raw q agrees
with diagnostic q exactly. These checks do not require a model/hardware run.

| Source | Physical valid rows | Matched scored rows |
|---|---:|---:|
| Circle | 7479/7501 | 6767 |
| Helix | 5773/5801 | 5030 |
| Sine,0g | 4229/4241 | 3866 |
| Sine,100g | 4231/4241 | 3984 |

Global coordinateRMSEmm, nominal versus augmented:
circle[4.296,4.714]vs[.748,1.220];helix[4.388,4.534]vs[1.022,1.180];
sine0g[3.696,2.897]vs[.482,.429];sine100g[4.505,3.733]vs[.729,.809].
These are scored-sample prediction errors, not physical force, inverse-pressure
accuracy, controller tracking improvement or closed-loop stability evidence.

## Retained/reproducible artifacts

- plotting_data.mat: full measurements/predictions/angles/errors/masks/metadata.
- *_plot_data.csv: all curves/errors/flags and both timebases. Matrix columns
  theta_rad/phi_rad are [measured,nominal,augmented]; error columns are
  [nominal,augmented]. Other two-column lengths are[l2,l3], inmm.
- prediction_metrics.csv, recording_summary.csv, selected_windows.csv,
  signal_inventory.csv and verification.json.
- Input/source SHA256 and artifact hashes; plotting-source snapshot.
- Adjacent ZIP contains the complete final bundle. v1 is retained as a
  preliminary export: curves unchanged; v2 fixes common-row theta metric
  comparison and standardizes circular-error axes. Usev2.

From project root, with a fresh destination, offline only:

```matlab
addpath(fullfile(pwd,'residual_controller'));
plotResidualShadowMotionPresentation(pwd,fullfile(pwd,'tmp','new_shadow_plots'));
```
