# Smooth candidate09 rejected; final budgeted search point is smooth_candidate10

Reviewed2026-10-28. Original filter restored: exact recorded-input velocity replay
matches pacMeasuredDerivativesOriginal with max error0. Gains/software reset and
reference match; no live edits/protocol mismatch, faults or pressure limiting.
Physical reset/calibration per standing user attestation. Source log unchanged.

Last3circles110-140s: l2/l3 RMSE[2.52791,2.75200]mm fail2mm; vector3.73682mm,
velocity.00351522m/s and fixed-analysis roughness.02207634m/s^2 pass their caps.
Scoring coverage99%,pair98%; tracking coverage99%,adaptation98%. This run is
neither more accurate nor smoother than the qualified original baseline
(~.0193roughness with1.55/1.87mm lengths). Do not confirm it. Low hysteresis gamma
also appeared in other slow-learning runs, but simultaneous gain changes and
one run cannot isolate causality. Last-lap lengths[2.3255,2.5333]mm remain outside
the fixed-duration targets even as error declines.

Final current-budget proposal (exploration, not a demonstrated improvement):

```matlab
k_g = 0.25306817650802854;
d_g = 0.00069737816290234929;
h_g = 2.6573766973344188;
```

Full snapshot in candidates/smooth_candidate10. Prediction is normalized objective
mean.81214/std.15115 (~.02030/.00378m/s^2). EI includes uncertainty; the mean alone
does not predict improvement over the baseline. Binary feasibility model does
not supply individual length-error bounds. No physical stability guarantee.

Use original measured.20/.30s and desired.10/.15s filters. Change only gamma
sources. All lambda18/16.5,Ks15/10,bounds/rates/history/pressure activation.2/
slew6/switches remain fixed. Full reset/no payload/phi1/12circles10s,150s total.
Save out in parent folder as smooth_candidate10_00g_phi1_T10_12cycles_r01.mat.
Preserve interrupted runs with notes. Send for review before any further run.

Nine observations counted; one search slot left. After its import the campaign
must stop for budget/search review. Two fresh confirmation repeats would require
a separately reviewed budget; previous screened data are not confirmation.
No live gain application,model/hardware execution,controller edits,commit/push.
