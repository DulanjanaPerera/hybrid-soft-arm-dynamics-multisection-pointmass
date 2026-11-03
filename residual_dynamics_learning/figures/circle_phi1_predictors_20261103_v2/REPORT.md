# Requested phi=1 circle: five existing forward predictors

2026-11-03. User requested a visual circle comparison for the five predictors
previously listed. No new predictors, fitting, optimization, hardware/Simulink
execution or controller experiments were added. Frozen files are hash-checked.

Source: 20261102/circleSweep.mat and its verified export. Selection uses logged
phaseId=5 and phi_reference_rad=1, identifying trial3, five circles between
simulation259.00 and358.95s,20s per requested circle. Observed phi across this
section is .05237-.23015rad, not1rad. This was recorded feedforward collection,
not five different controller tracking runs. The requested reference, measured
response and each model's prediction are explicitly separate.

## Full-circle figures

full_circle_xy and full_circle_angles use the FIRST chronological circle,
259-279s simulation,20.23s physical poll-start duration,401samples including the
endpoint. Initial10frame measured q/v history is supplied once. Subsequent
q/v come only from predicted motion and the original causal velocity update;
future measured states are used for scoring/plotting only. All models use the
same recorded pressure commands and physical query intervals. Inputs use the
original sensor-row ZOH approximation, not measured chamber pressures.

Models are exactly nominal legacy-bound local-affine transition, constant
velocity, nominal+ridge, kinematic+ridge and nominal+GRU. Legacy contextmode0
preserves frozen training semantics. No corrected-context learned model or
new architecture was substituted. This nominal propagation is not exact
nonlinear ODE integration. All five complete the circle within the inherited
numerical q-coordinate guard; this is not a stability/safety qualification.

XY uses the SAME PCC mapping on predicted and measured orientation-inferred
lengths: radius.013m, module length.17411m, u=(l2+l3)/r,
v=(l3-l2)/(sqrt(3)r), phi=hypot(u,v),theta=atan2(v,u),
rho=L*(1-cos(phi))/phi. It is not an independent raw NDI XYZ comparison.
The logged orientation versus inferred-length inverse map roundtrip passes
within1e-10rad. Theta curves are unwrapped and lifted to the reference's
initial2pi branch for display; reference angles are not prediction features.

| Predictor | q2/q3 vector path RMSE, mm | PCC XY vector RMSE, mm | phi MAE, rad |
|---|---:|---:|---:|
| Nominal | 11.489 | 70.727 | .8560 |
| Constant velocity | 2.253 | 15.490 | .0465 |
| Nominal + ridge | 3.825 | 22.577 | .2423 |
| Kinematic + ridge | 4.750 | 28.185 | .3074 |
| Nominal + GRU | 1.096 | 6.521 | .0552 |

These errors compare prediction with measured RESPONSE, not requested reference.
The nominal model follows the requested ideal circle but poorly predicts the
measured small-radius response. GRU is closest in length/XY on this one selected
circle. Constant velocity has smaller phi MAE but does not rotate through the
circle, illustrating why bend magnitude alone is insufficient. Ridge corrections
drift over20s more than the aggregate5s results suggested. This is one circle,
not a newly selected architecture winner or controller demonstration.

## One-step figure

one_step_five_circles_angles_labeled shows the EXISTING one-step predictions
across all five circles. Measured history is supplied again at EACH sample;
this is not an autonomous100s rollout. Keep this distinction explicit when
showing near-overlapping learned curves. The GRU's isolated one-step excursion
is retained, not removed or smoothed. The final labeled variant explicitly
allocates its legend outside the tiled axes; the original export is preserved.

## Artifacts and reproduction

Three plot sets (final one-step set uses _labeled suffix), editable FIG, vector
PDF and300dpiPNG. All prediction/measurement/reference q/angles/PCC coordinates
and source-row/time IDs are in CSV; MAT stores initial protocol and predictions.
sources.csv records raw/export/normalizer/weight/data SHA256, checked unchanged.
source_snapshot contains plotting/prediction/mapping sources, plus the rendering
finalizer. V1 is an incomplete table-shape attempt marked DO_NOT_USE.

From repository root, choose a NEW output directory:

    addpath(pwd,fullfile(pwd,'residual_dynamics_learning'));
    plotResidualCirclePredictors(pwd,newFolder);
    finalizeCirclePlotLayout(newFolder);

The frozen models remain next-position corrections. This figure request does
not implement a force residual or authorize running predictors on the arm.
