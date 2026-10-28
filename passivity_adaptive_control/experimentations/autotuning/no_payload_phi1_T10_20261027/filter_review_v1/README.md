# Measured derivative filter review and smoothness-first campaign preparation

Outcome update: actual_arm_comparison/README.md reviews the user's12circle
new-filter run at candidate08 gains. It worsened fixed-analysis roughness and
length error. Restore pacMeasuredDerivativesOriginal.m (.20/.30s); the experimental
replacement is not selected. Historical offline results below remain for audit.

2026-10-28. User prefers smoothness over accuracy and approved l2/l3 RMSE ceilings
of2mm each. The existing campaign_length_v1, its1mm targets, trials, artifacts and
ledger remain unchanged. No Simulink file was edited or hardware/model executed.

## What is actually in the saved model

Read source XML from PAC_autotuning_manual.slx without loading/executing the model:

- Measured filteredDerivatives block: backward difference of raw lengths, then
  first-order velocity tau=.20s; difference of filtered velocity then tau=.30s
  acceleration smoothing. Invalid/repeated frames hold states.
- Desired lengthDerivatives: velocity tau=.10s, acceleration tau=.15s.
- pacAdaptiveForce uses measured q/dq, desired q/dq/ddq. Measured acceleration
  is not used in its force law. It is still useful for diagnostics, but smoothing
  it alone cannot change the current PAC force. Noisy q/dq, desired derivatives,
  actuator dynamics and actual stop/start motion remain separate possible causes.

Saved measured source is retained as pacSavedMeasuredDerivatives.m. Replaying
baseline1 and candidate08 with recorded pollStart_s reproduces logged velocities
exactly (max difference0). Replacing pollStart_s with controller time does not;
do not confuse these clocks when differentiating the physical samples.

## Proposed measured-filter replacement

File: passivity_adaptive_control/pacSmoothMeasuredDerivatives.m (relative to repo).
Standalone codegen-oriented MATLAB Function source, same three inputs/outputs as
the measured block. Two equal velocity-filter poles, exact ZOH propagation,
tauV=.15s per stage; acceleration smoothing tauA=.35s. Raw length output remains
unchanged. Invalid/repeated/backward timestamps hold derivative states. Long
gaps>.15s re-anchor position and hold derivative states for that frame; the
existing supervisor remains responsible for freshness/fault behavior.

This is an offline-tested candidate, not a physically validated improvement.
Derivative filtering can add delay and change stability margins. It cannot prove
jerks are caused by noise or distinguish real high-frequency motion from noise.
No zero-phase/future-data filtering is used in the live candidate.

## Replay results, not new physical experiments

Late110-140s results on the same raw recordings, using a10s sin/cos/DC fit:

| Recording/metric | Existing | Candidate |
|---|---:|---:|
| Baseline1 velocity non-fundamental RMS(m/s) | .002418 | .002234 |
| Candidate08 velocity non-fundamental RMS(m/s) | .002740 | .002439 |
| Baseline1 acceleration non-fundamental RMS(m/s^2) | .005150 | .004148 |
| Candidate08 acceleration non-fundamental RMS(m/s^2) | .006286 | .004808 |
| Baseline1 velocity filter lag l2/l3(s) | .216/.223 | .289/.299 |
| Candidate08 velocity filter lag l2/l3(s) | .226/.221 | .298/.296 |

Non-fundamental content includes actual motion harmonics as well as noise.
Do not call its reduction a measured reduction in physical jerk. The initial
tauV=.10/tauA=.20 prototype increased non-fundamental content in these recordings;
its summary is preserved in comparison_initial_tau010_tau020.csv. Current outputs
comparison.csv/replay MAT/PNG files correspond to the .15/.35 candidate.

compareMeasuredFilters reproduces the replay; it overwrites only its own analysis
outputs, not source recordings/model files. validatePacSmoothMeasuredDerivatives
checks ramps, raw-position passthrough, timestamp/gap/invalid behavior and synthetic
6Hz attenuation. Real control-loop phase/performance remains to be tested.

## Fair smoothness scoring and broader search draft

pacTuneSmoothCampaignDraft returns a non-active purpose=draft configuration:

- Primary objective: velocityErrorRateRmse_mPerS2, normalized by .025.
- User-approved position constraints: l2 and l3 RMSE<=2mm each.
- Velocity<=.01m/s, zero pressure limiting; provisional existing roughness cap.025
  carried forward for baseline review, not treated as a newly validated threshold.
- Proposed broader gamma bounds: k[.10,1],d[.0005,.008],h[1,10]. Damping widens
  downward; its upper value stays at the previously approved .008 because high
  gammaD has already coincided with roughness. These ranges are proposals, not
  physically validated stable bounds. No new candidate has been queued.

To prevent a feedback filter from gaming the score, the new draft uses raw lengths
and a FIXED OFFLINE analysis differentiator: measured tau=.20s, reference tau=.10s,
same across all compared runs. It does not use changed live velocity/acceleration
plots to claim smoother motion. Those are still retained as control diagnostics.
The roughness proxy remains filtered velocity-error slope, not an independently
measured mechanical jerk. Inspect raw motion/pressure traces with the same
analysis filter, and reassess acceptable lag and2mm tracking on the arm.

Old configurations with no objective/source fields retain their length objective
and logged-velocity scoring. New objectiveMetric and motionMetricSource fields
are immutable within a campaign; do not edit campaign_length_v1 in place or train
the new filtered-controller campaign as though old-controller trials used it.

## Next operator action

1. Save/preserve the current model and keep the existing operator stop/kill method.
2. Paste the entire pacSmoothMeasuredDerivatives.m into the MEASURED
   filteredDerivatives MATLAB Function block (saved model name MATLAB Function5).
   Same input/output ports; no extra wiring. Leave desired lengthDerivatives
   unchanged for this first test. The assistant has not made this edit.
3. Restore original baseline gamma k=.5,d=.002,h=5 and keep lambda18/16.5,
   Ks15/10, activation.2,slew6,bounds/rates/history/switches unchanged. This gives
   a filter-only comparison against existing baseline recordings, not a
   simultaneous new-filter/new-gamma comparison.
4. Full reset/no payload/phi1/T10; run3circles with normal startup, ramp-down and
   final hold. Save out in the parent folder as
   filter_smoke_00g_phi1_T10_3cycles_r01.mat. Send for review before optimization.
5. If the smoke test is acceptable, acquire a12circle filter baseline and assess
   length error, lag, fixed-analysis roughness, pressure/faults and software reset.
   Then finalize the smoothness-first campaign version/bounds. Do not resume the
   old candidate09 suggestion or automatically extend the previous budget.

These are manual-test instructions, not a hardware execution authorization for
the assistant. No guarantee is made for the eventual2s trajectory; filtering and
latency must be reassessed at that speed. If roughness remains, inspect feedback
lambda/Ks and actuator/measurement dynamics instead of increasing gamma alone.
