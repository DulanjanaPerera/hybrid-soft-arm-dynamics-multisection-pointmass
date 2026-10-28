# Smoothness-first wider-range campaign: smooth_candidate09

Update: smooth candidates09/10 reviewed; budget exhausted, no next run queued.
See completed_search_review/README.md. Original instructions below are retained
for audit only. Additional trials need budget review.

Prepared2026-10-28 at user request. Old campaign_length_v1 remains unchanged.
Eight prior ORIGINAL-filter recordings were replay-verified and rescored here;
the changed-filter recording is excluded. Historical reuse is training, not fresh
confirmation. Physical reset/calibration metadata follows user attestation.

Minimize fixed-analysis velocity-error-rate roughness, with length RMSE<=2mm
each,velocity<=.01m/s,roughness<=.025m/s^2,zero limiting and95% coverage.
Search bounds: gammaK[.10,1],gammaD[.0005,.008],gammaH[1,10]. Keep the fixed
offline raw-length analysis filter (measured.20s/reference.10s) for fair scoring.

The original baseline repeats are best screened by mean roughness (~.0193m/s^2).
With2mm ceilings they already satisfy the new criteria; this is retrospective
rescoring, not a new confirmation. optimizer.continueSearchAfterFeasible=true
explicitly allows optimizing further in this new campaign. Old defaults remain
confirmation-first. No proposal is declared a confirmed optimum.

## Next manual run

First ensure the ORIGINAL measured filter is restored: velocity tau.20s,
acceleration tau.30s; paste-ready pacMeasuredDerivativesOriginal.m is in the
passivity_adaptive_control folder. Desired filter stays.10/.15s. Change only:

```matlab
k_g = 0.90453500157774203;
d_g = 0.0022428161211506688;
h_g = 2.8192356276278843;
```

Full tuning: candidates/smooth_candidate09/review.md and candidate.mat. Preserve
precision. Lambda18/16.5,Ks15/10,bounds/rates/scales/history,pressure activation.2/
slew6 and switches stay fixed. Full reset,no payload,phi1,T10,12circles:
startup0-20,tracking20-140,ramp140-145,pressurizedhold145-150s,then manualkill.
Keep the existing operator stop method and preserve interrupted runs with notes.

Save out in the PARENT experiment folder as:
smooth_candidate09_00g_phi1_T10_12cycles_r01.mat

Send the recording before another candidate. This is a fitted-surrogate
feasibility_weighted_expected_improvement proposal, not blind initial design.
Predicted normalized mean.72570/std.11699 correspond to roughness about
.01814/.00292m/s^2. These are model-dependent predictions; heuristic feasibility
.58450 is not a physical success probability or stability guarantee.

Budget remains10observations:8historical+2new exploratory slots. Fresh
confirmations need budget review afterward; no automatic extra runs. Do not run
the old unqueued candidate09 suggestion from campaign_length_v1: its objective
and proposal differ. No hardware/model execution, live gain application,
controller/SLX edits, commit or push performed by the assistant.
