# November3 short NDI gap audit

Offline feasibility preview only. Raw recordings, exports, validity flags,
training masks, all frozen weights and November2 evaluation are unchanged.

57 invalid NDI samples in 55 episodes.
53 one-frame episodes meet the conservative preview rule:
one missing sample, valid clean endpoints, positive same channel/trial,
base still valid, individual physical intervals <=.1s and bracketing
valid-to-valid span <=.12s. Two-frame episodes remain missing.
Preview interpolates q2/q3 on actual host timestamps; q1=-q2-q3.
Pressure commands are never interpolated. No held-invalid value is
silently reclassified as a real measurement. Each preview row has
originalNdiValid and imputed flags; it has NOT entered any model pipeline.

## Reconstruction proxy test

Real valid blocks were artificially hidden, including neighborhoods of
actual dropouts. Linear reconstruction uses bounding observations.
Errors below compare with the hidden real measurements. They do not
measure error at actual gaps, whose truth is unavailable. Overlapping
test blocks are correlated; this is not an independent validation set.

| Recording | Hidden frames | Windows | Vector RMSE (mm) | p95 norm error (mm) | Max norm error (mm) |
|---|---:|---:|---:|---:|---:|
| circleSweep | 1 | 604 | 0.0711 | 0.1602 | 0.4136 |
| circleSweep | 2 | 590 | 0.1080 | 0.2372 | 0.5811 |
| SineWave | 1 | 499 | 0.0677 | 0.1317 | 0.4855 |
| SineWave | 2 | 495 | 0.0805 | 0.1547 | 0.7591 |
| helix | 1 | 597 | 0.0980 | 0.2198 | 0.4585 |
| helix | 2 | 586 | 0.1438 | 0.3292 | 0.6633 |
| StepWave | 1 | 482 | 0.0500 | 0.1037 | 0.3719 |
| StepWave | 2 | 477 | 0.0712 | 0.1674 | 0.5227 |
| combinePressure | 1 | 524 | 0.0401 | 0.0856 | 0.1695 |
| combinePressure | 2 | 516 | 0.0633 | 0.1374 | 0.4498 |

## Recommendation and next experiment

Short-gap interpolation is plausible as a separate offline training
input experiment, not an unconditional replacement for validity.
Keep the original measured-only target mask; do not train force or
acceleration truth labels from filled lengths. Retain long-gap resets
and time/pressure/calibration/trial boundaries. Preserve unchanged
measured-only validation, and compare causal operation under simulated
dropouts before accepting a new recurrent-state policy.

Two-sided interpolation uses the next valid sample, so it is not
available at zero delay online. A deployed estimator must predict/hold
causally, expose uncertainty and missingness, or explicitly accept delay.
Masking only the filled targets does not remove future leakage if
imputed inputs already advanced the validation hidden state.

An alternative is gap-aware physics rollout with loss only on real
observations. Free predictions can already traverse a missing target
using known commands/time without imputing measured q. Measurement-
conditioned training needs an explicit tested missing-input/state
policy; the current frozen reset convention has NOT been changed.

Do not interpolate wrapped theta near straightness. If raw rotations
are later used, quaternion SLERP is preferable to blending Euler angles.
Smoother appearance alone is not evidence of correct derivatives or dynamics.

References: MATLAB fillmissing supports physical SamplePoints, MaxGap
and returned fill masks: https://www.mathworks.com/help/matlab/ref/fillmissing.html
Quaternion SLERP: https://www.mathworks.com/help/uav/ref/quaternion.slerp.html
