# Directional asymmetry in the recorded nominal-model bend prediction

2026-11-02. Read existing SineWave/SineWave2 exports from
audits/validity_20261102_v1, backed by the original20261102 MAT recordings.
No dynamics integration, parameter fitting, learned-model scoring or hardware
execution. No mechanical helper or SLX changed.

## Finding

The approximately3rad versus2.5rad peaks match the logged simulated des_phi,
not the orientation-derived measured phi. The two recordings have identical
nominal-model peaks on the sensor grid:

| Excited channel | Model peak phi (rad) | Measured peak phi, SineWave / SineWave2 (rad) |
|---|---:|---:|
| 1 | 3.005162 | 1.952151 / 1.950837 |
| 2 | 2.526838 | 1.935315 / 1.931885 |
| 3 | 2.526837 | 2.493593 / 2.494099 |

These maxima are not necessarily simultaneous measured/model instants. Separate
columns in peak_summary.csv retain model-peak time, lengths and bound forces.
The plot shows the first complete sine trial per channel, not every repeat.

## Identified model asymmetry

The current saved forwardModel sets lengthBounds=[-0.035;0.035;1e6] andmu2000,
calls armS_single_entry, reconstructs l1=-l2-l3, then computes

    phi = hypot((l2+l3)/r, (l3-l2)/(sqrt(3)*r)), r=0.013 m.

There is no phi clipping at2 or2.5rad in this forward calculation. The other
orientation reader computes a principal tangent angle in[0,pi]; the XYZ-IK
scalar fit uses[0,pi-1e-6]. Neither imposes a channel-specific2.5rad cap.

armS_single_dynamics applies its nonlinear bound force only to the two stored
coordinates l2 andl3. The reconstructed l1 receives no equivalent bound force.
At the recorded model peaks (near3bar):

| Channel | l1 (mm) | l2 (mm) | l3 (mm) | Existing additional bound force in q2/q3 (N) |
|---|---:|---:|---:|---|
| 1 | +39.067 | -19.534 | -19.534 | approximately[0,0] |
| 2 | -16.424 | +32.849 | -16.424 | approximately[6.020,0] |
| 3 | -16.424 | -16.424 | +32.849 | approximately[0,6.020] |

Thus channel1 can extend simulated l1 past+35mm without a penalty, whereas
channels2/3 encounter additional restoring force on their extending coordinate.
The tanh law is smooth and acts before the nominal35mm threshold;32.849mm is
not a hard length clamp. This identified artificial directional resistance
explains the simulated peak discrepancy; it is not evidence that the real arm
has exactly this directional stiffness or that the angle calculation is capped.

## Possible correction, not implemented

If the intended limits represent equal physical actuator-length bounds, apply
the same scalar restoring law F(l) to all three reconstructed lengths and map
forces into q=[l2;l3] using the length Jacobian:

    J = [-1 -1; 1 0; 0 1];
    f_bound_q = J.' * [F(l1);F(l2);F(l3)]
              = [F(l2)-F(l1);F(l3)-F(l1)].

The current term is only[F(l2);F(l3)]. A correction would require a versioned
nominal baseline and validation of the force/potential/coordinate conventions;
do not silently alter old fitted models, datasets or benchmark reports. This
report has not run a corrected-model counterfactual, so it does not claim an
exact corrected angle or that this alone fixes mechanical prediction error.
The known issue should be explicit before interpreting residual learning as
missing physical dynamics. No generic calibration detour is being restarted.

## Generator move and acquisition

Operator reports channel3 moves the tip toward the field generator and is
repeating with the generator moved back. This is a stated setup hypothesis/
change, not confirmed resolution of the earlier tip dropouts. Record the new
generator pose/setup as a separate run and repeat straight-pose calibration/
preflight: the current reader's gravity/reference calibration assumes a fixed
NDI frame. Preserve old runs and avoid treating setup changes as model learning.

## Artifacts

peak_summary.csv; per-run plotted-value CSVs; both three-by-two diagnostic
figures as FIG/vectorPDF/300dpiPNG; diagnoseRlBoundAsymmetry.m; saved forward
chart and current mechanical source snapshots/hashes. Inputs reference existing
quality-audit exports instead of duplicating original raw MATs. Diagnostic code
checks equal observation/model timelines and zero-sum reconstructed lengths.
Both plots distinguish prediction from measurement. Raw inputs are unchanged.
