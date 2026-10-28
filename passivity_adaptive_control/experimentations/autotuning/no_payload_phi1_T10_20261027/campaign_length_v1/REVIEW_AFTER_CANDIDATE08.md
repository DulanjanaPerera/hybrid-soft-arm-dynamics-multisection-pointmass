# Candidate08: best observed, targets not met; proposal review recommended

Reviewed2026-10-28. Imported candidate08_r01 with exact proposed tuning, matching
software initial estimates/history, no live edits/protocol mismatch. Physical
reset/calibration follows standing user attestation. Raw MAT remains unchanged.

Last3circles110-140s: l2/l3 RMSE[1.1725166,1.3028913]mm, vector1.7528037mm,
velocity.00299378m/s, roughness.02391251m/s^2. Both1mm position limits fail;
velocity/smoothness/zero-limiting limits pass under the agreed pooled scoring.
Coverage99.5%, pair99%. Tracking coverage99%, adaptation97.333%, faults0,
pressurelimits0, no tuning changes. Rare early stiffness bound/rate contacts
are reported in the full artifact. Vector score15.4% below07,27.9% below mean
baseline; one observation is not repeatable confirmation or proof of convergence.

Per-lap caveat: roughness lap10 .024674,lap11 .020047,lap12 .026547m/s^2.
The pooled3lap score passes .025, but not every late lap does. Do not change the
acceptance definition retrospectively; retain the caveat for later confirmation.
Last-lap length[1.1433,1.2684]mm also remains above target. Extending the run to
meet limits would be a different protocol, not a pass at the agreed120s duration.

There are8 recorded observations and2 slots remaining, with no accepted candidate
to confirm. Two fresh confirmations cannot fit if another search trial is first
needed. Stop/review rather than automatically extend the budget.

The current unchanged adapter returns seek_feasible_region because every observed
candidate has at least one violated limit. Raw next suggestion is saved ONLY in
suggestion_after_candidate08.mat; it has NOT been queued as candidate09:
gammaK .26075453527454223,gammaD .0010338149762252346,gammaH3.9294515860178212.
GP predicted objective2.69475mm,std.63161mm; heuristic feasibility.49796 is not
a calibrated physical success probability. This proposal does not predict
improvement over08's measured1.75280mm; it favors an uncertain region.

Algorithm limitation: binary feasible/infeasible labels cannot distinguish a
candidate narrowly missing1mm from one with poor length error or roughness. With
no feasible incumbent, the current acquisition uses feasibility score and GP
uncertainty, not objective improvement. Suggest reviewing component constraint
margins before using the final evaluations. A margin-aware acquisition would be
a new algorithm/campaign version, retaining these observations and clearly
documenting the change, not an in-place edit of the immutable campaign config.

No candidate09 proposal event, live gain application, hardware/model execution,
controller edits, commit or push by assistant. Status_after_candidate08 and
artifacts/candidate08_r01 hold complete reports. The campaign's automatic status
may say awaiting_candidate; that is software bookkeeping, not a direction to run
the unscreened suggestion or an authorization for more hardware trials.
