# Candidate06 reviewed; next manual run is candidate07

Reviewed2026-10-28. Candidate06 exact tuning/software initialization match,
no live edits/protocol mismatch. Physical reset per standing user attestation.
Late110-140s length RMSE[2.63375,2.90327]mm/vector3.91990mm; velocityRMS
.00344364m/s; roughness.01935062m/s^2. Both length targets fail; velocity,
smoothness,zero-limiting requirements pass. Coverage98.833%, pair97.667%.
Tracking coverage98.958%,active97.917%,faults0,pressurelimits0. No rate cap
contacts or upper bounds. Length error continues decreasing throughlap12
([2.4442,2.7215]mm), but remains much larger than candidate03 at the same duration.
Smoother than05 but less accurate. Do not attribute the tradeoff to one gamma:
all three changed and physical conditions are not independently held/measured.

Candidate03 remains best observed (vector2.21092mm), not accepted/confirmed.
Next is the last seeded initial-design point:

```matlab
k_g = 0.35438912779204967;
d_g = 0.0030892129738313034;
h_g = 6.8259366991000565;
```

Change only those gamma sources before the same full-reset procedure. Keep all
other gains,bounds,rates,history,switches and pressure activation.2/slew6 fixed.
No payload,phi1,12circles10s,complete150s sequence:20startup,120tracking,
5ramp-down,5pressurizedhold then manualkill. Save out in the parent folder as:
candidate07_00g_phi1_T10_12cycles_r01.mat

Send recording before the next candidate. Six observations counted,four remaining
including confirmations. After this point's observation, remaining new proposals
can use the fitted surrogate acquisition rather than seeded initial design, if
targets remain unmet and budget allows. No gain/performance guarantee. No live
gains,hardware/model execution,controller edits,commit or push by assistant.
Artifacts/candidate06_r01 and status_after_candidate06 hold the audited results.
