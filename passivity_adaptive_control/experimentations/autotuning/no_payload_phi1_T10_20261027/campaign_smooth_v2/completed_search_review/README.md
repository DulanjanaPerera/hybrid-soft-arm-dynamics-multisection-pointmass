# Completed smoothness-first search: budget exhausted

Reviewed2026-10-28. Eight historical screenings were rescored and two new
original-filter exploratory trials evaluated. No fresh confirmation trials were
run. Campaign status evaluation_budget_exhausted,remaining0,success=false.
No candidate11 was generated. The original campaign is preserved separately;
its unused slots are not permission to bypass this campaign's budget review.

## Final candidate10

smooth_candidate10_00g_phi1_T10_12cycles_r01.mat matches proposed tuning and
software initial state. Original filter verified by recorded-input replay,max
velocity difference0. No live edits/protocol mismatch. Physical reset/calibration
follows standing user attestation. Tracking coverage99.083%,active98.167%,
faults/pressurelimits0. Last3circles110-140s:

- l2 RMSE2.880204mm,l3 RMSE3.166379mm:both fail2mm.
- Vector length RMSE4.280366mm.
- Velocity RMS.003572203m/s:passes.01.
- Fixed-raw roughness.016077616m/s^2:passes.025 and is lowest measured.
- Coverage99.5%,validpairs99%.
- Lastlap length[2.6794,2.9826]mm:still outside targets despite ongoing learning.

Roughness is about16.7% below the baseline-repeat mean, but length errors are
unacceptable under the chosen ceiling. This is a smoothness/accuracy tradeoff,
not an accepted improvement or physical parameter identification result.

## Best qualified observed reference

Original baseline gammaK=.5,gammaD=.002,gammaH=5 with original measured filter
tau.20/.30s and desired.10/.15s remains best screened by roughness under2mm
constraints. Two historical runs have length RMSE[1.5684,1.8887]/[1.5290,1.8569]mm
and fixed-raw roughness.01937096/.01922323m/s^2. These were used in training and
are not held-out confirmations. All lambda18/16.5,Ks15/10,activation.2,slew6,
bounds/rates/history/switches remain as recorded.

Other feasible measured points03,04,07,08 improve some accuracy measures but
are rougher. Candidate08 is best length tracking, not best smoothness. Wider
range points09/10 did not improve the constrained objective over the baseline.
Single runs, simultaneous gamma changes and unmeasured physical drift limit
causal conclusions. No global optimum or stability certificate is established.

## What to review before further experiments

Keep the original filter. Do not expand gamma ranges or weaken criteria merely
to label an infeasible candidate successful. A future search should assess
individual coordinate constraint margins and local experiments around the
feasible baseline, and distinguish slow adaptation from a steady tracking floor.
Inspect actuator/feedback behavior if roughness persists; gamma alone does not
set closed-loop bandwidth. A longer run changes the fixed-duration objective
and must be recorded as a distinct protocol. Any additional search or two fresh
confirmation repeats requires an explicitly reviewed budget/new campaign version.

Artifacts: all_observations.csv,review.mat,accuracy_smoothness_tradeoff.png and
reviewSmoothBudget.m in the parent experiment folder. Source MATs unchanged.
No live gains applied,controller/model edits,hardware execution,commit or push
by the assistant. This closes data evaluation, not the broader control objective.
