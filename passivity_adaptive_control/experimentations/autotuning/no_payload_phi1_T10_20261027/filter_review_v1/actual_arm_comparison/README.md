# Actual-arm filter comparison: retain the original filter

Reviewed2026-10-28. User supplied filter_smoke_00g_phi1_T10_12cycles_r01.mat,
a12circle150s run despite the smoke filename. Actual gamma=[.75,.008,7.5],
lambda18/16.5,Ks15/10. This matches candidate08's gains to practical numeric
precision, not the suggested original baseline gamma. Therefore compare against
candidate08, not the original baseline. No live tuning changes, no tracking
faults/pressure limits. Software initial-state checks remain in per-run reports;
physical reset follows standing user attestation.

Exact replay using recorded pollStart_s reproduces new recorded velocity with
pacSmoothMeasuredDerivatives max error0; old filter replay differs by up to
.008204m/s. The new filter was actually installed. No saved model was modified
or hardware executed by the assistant during this review.

Both runs evaluated using the SAME fixed raw-length analysis filter,110-140s:

| Metric | Original filter/candidate08 | New filter |
|---|---:|---:|
| l2 RMSE(mm) |1.17252|1.18936|
| l3 RMSE(mm) |1.30289|1.41650|
| Vector RMSE(mm) |1.75280|1.84961|
| Velocity error RMS(m/s) |.00300670|.00357242|
| Roughness(m/s^2) |.02396971|.03408892|
| Theta-rate SD(rad/s), secondary |.39579|.52262|
| Near-stop fraction, secondary |8.08%|10.10%|
| Phi peak(rad), secondary |1.1673|1.2111|

New roughness42.2% above old and exceeds the provisional.025cap. Length error
also worsens modestly. User perceived no improvement in smoothness; evidence
supports restoring the original measured filter. Added phase lag is a plausible
contributor, not a proven causal explanation from one trial. Offline smoothing
of fixed input data was not a closed-loop performance guarantee.

## Restore instructions

Paste the complete passivity_adaptive_control/pacMeasuredDerivativesOriginal.m
into the measured filteredDerivatives MATLAB Function block, same3inputs/outputs.
It is an exact behavioral copy of the saved original source: velocity tau.20s,
acceleration tau.30s, original frame/timestamp handling. Leave desired derivative
filter at.10/.15s. No need to increase derivative filtering further now. The
assistant prepared source only and has not edited/applied anything to Simulink.

The experimental two-pole source remains for audit but is not selected. The
smoothness-first draft now specifies the ORIGINAL filter. User-approved2mm
coordinate ceilings and roughness-primary objective remain; broader proposed
bounds are still only a draft. Continue fixed raw-length analysis across runs.
Do not ingest this changed-filter experiment into the original-filter optimizer
as if the plant/controller configuration were unchanged. It is recorded here,
outside campaign_length_v1; its ledger remains8observations/2remaining slots.
No new optimizer candidate, campaign, commit or push is created by this review.
