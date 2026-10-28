# Candidate07 is best observed; candidate08 is boundary exploration

Reviewed2026-10-28. Candidate07 exact tuning/software initialization match,
no live edits or protocol mismatch; physical reset per standing user attestation.
Late110-140s l2/l3 RMSE[1.23668,1.66263]mm, vector2.07213mm, velocity.00299571m/s,
roughness.02258509m/s^2. Length limits fail; velocity/smoothness/limiting pass.
Scoring coverage98.5%, pair97%; tracking coverage98.625%, active97.042%, faults0,
pressure limits0 and rate caps0. Lap12 length[1.2090,1.5526]mm; no claim of
settling/convergence. Vector score6.3% below03, but not a repeated confirmation.

Seven observations at6 unique gamma triples complete the initial design.
First fitted-GP proposal has reason seek_feasible_region, since none meets all
constraints. Prediction mean2.56631mm/std.60053mm (latent/model-dependent) does
not predict improvement over07. Feasibility frequency.49866 is heuristic, not a
physical success probability. The acquisition favors uncertain, less-sampled
regions when no feasible point exists; all infeasibility labels are binary and
do not distinguish how close a candidate is to individual constraint limits.

Candidate08 is therefore explicitly an upper-bound exploration trial, not a
recommended best setting. High gammaD near this boundary failed smoothness in05;
other gamma values also differ, so that observation does not identify the causal
effect of gammaD alone or guarantee how this combination will behave.

Exact archived proposal (preserve precision for intake):

```matlab
k_g = 0.75;
d_g = 0.0079999999999999967;
h_g = 7.4999999999999991;
```

Change only these gamma sources. Keep all other gains/bounds/rates/history/
switches/pressure settings fixed. Full reset, no payload,phi1,12circles10s,
20startup+120tracking+5ramp-down+5pressurizedhold then manualkill (150s total).
Preserve an interrupted/abnormal run with notes rather than discarding it.
Save out in the parent folder as candidate08_00g_phi1_T10_12cycles_r01.mat.
Send recording before another candidate. Seven observations/three remaining,
including two required confirmations; if the budget is insufficient, stop/review
rather than launch extra trials. No gains applied or hardware/model execution,
controller edits, commit or push by assistant. Reports: artifacts/candidate07_r01,
status_after_candidate07, candidates/candidate08/review.md.
