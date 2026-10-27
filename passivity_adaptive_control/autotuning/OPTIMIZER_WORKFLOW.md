# Bounded PAC optimizer adapter

Implemented 2026-10-27 (milestone C). This adapter proposes grouped gammaK,
gammaD and gammaH from a manual campaign's observations. It never applies gains,
runs models, acquires experiments or resets the arm. Read MANUAL_WORKFLOW.md
first. Actual arm campaign settings and targets remain unsettled.

## Use with a manual campaign

```matlab
addpath('passivity_adaptive_control','passivity_adaptive_control/autotuning');
p = pacTuneSuggestNext('path/to/campaign');
% Only when p.action is 'propose':
candidate = pacTuneQueueSuggestion('path/to/campaign','candidate02',p);
```

| Returned action | Next step |
|---|---|
| `blocked` | Resolve incomplete campaign details/targets; no gamma is recommended |
| `observe` | Import an operator-run observation for the pending candidate |
| `confirm` | Obtain reserved confirmation repeats for the screened candidate |
| `propose` | Review the proposed gamma and queue its tuning snapshot |
| `stop` | Budget exhausted, confirmed targets met, or finite search pool exhausted |

Queueing archives settings, training counts, predictions, reason and ledger/config
identity in the candidate artifact. It uses the existing writer lock/bound checks
and does not write live gains. If the ledger changes between suggestion and queue,
the write is rejected; recompute rather than editing identity fields. Pending
observations and confirmation take priority over further search. There is no
automatic experiment loop or automatic recovery after a fault.

## Immutable configuration and restart

New manual templates contain an `optimizer` snapshot from
`pacTuneOptimizerConfig`. Set it before campaign creation. The adapter explicitly
maps campaign gamma bounds, grouped baseline and evaluation budget into this
snapshot. Baseline gamma must be equal within each k/d/h group; bias gamma stays
fixed. Older campaigns without an optimizer snapshot need a new campaign version;
do not edit campaign.mat. Incomplete manual campaign details keep proposals blocked.

State is reconstructed from immutable config and append-only events. A local RNG
makes proposals independent of the caller's global RNG. Config plus raw observations
and pending points reproduce a proposal without a live optimizer session. Candidate
artifacts retain algorithm schema `pac-gp-1`; archive source/Git revision, and version
and revalidate future algorithm changes for subsequent campaigns.

## Algorithm and assumptions

1. Transform positive gamma bounds to normalized log coordinates; fixed dimensions
   stay fixed. Zero adaptation is a separate manual baseline/freeze experiment.
2. Propose the exact baseline, then a seeded stratified initial design (six points
   including baseline by default).
3. Fit a small GP with fixed squared-exponential kernel to usable screening scores.
   Center/scale objectives; retain repeated measurements and use their mean/sample
   variance with a nonzero noise floor. Hyperparameters are not optimized.
4. Search a seeded finite pool (1500 points plus corner points). Use expected
   improvement against the best measured feasible objective, weighted by a
   kernel-based feasibility frequency. Before feasibility is observed, prefer
   feasibility score and GP uncertainty; with no usable objective, use distance
   exploration. Numerical EI underflow has an uncertainty/distance fallback.
5. Exclude near-duplicate seen/pending points in normalized log space. Deliberate
   operator-requested screening repeats remain supported via manual imports; the
   optimizer does not automatically retry a fault. Stop at budget or pool exhaustion.

MATLAB bayesopt and a Statistics license were found locally. This adapter uses its
own small proposal-only GP to avoid an objective callback that might launch trials.
It uses base MATLAB numerics/RandStream, not the optimization toolbox. Reading
existing recordings still requires their Simulink data classes.

The feasibility score is a Beta(1,1) kernel-weighted label frequency: a heuristic,
not a calibrated probability or safety certificate. GP predictedStd is latent
uncertainty under fixed kernel/noise assumptions, not a guaranteed tracking interval.
Finite-pool search does not prove a global optimum or stable arm operation.

## Evidence and failure handling

Confirmation observations are reserved and never train either model; manual status
still uses them to accept or invalidate confirmation. Screening rules:

| Evidence | Objective | Constraint evidence |
|---|---|---|
| Comparable valid_feasible | Measured score | Feasible |
| Comparable valid_infeasible | Measured score | Infeasible |
| Comparable logged latched fault | No partial score | Infeasible |
| Early end without verified fault | None | Unknown |
| Missing data, live edits, tuning/reference/reset/metadata mismatch | None | Unknown |

Protocol checks also apply when a fault masks the protocol-deviation classification.
All observations remain in the ledger and count against the actual campaign budget.
Excluded screening attempts remain seen points without invented numeric penalties.
Total-budget checks include confirmations even though training excludes them.
Scores use fixed normalization/windows; periods and payloads are not pooled.

Suggestions expose prediction mean/std, feasibility score, acquisition, config,
usable/excluded counts and repeat statistics. The replicates table includes mean,
sample std and standard error at exact gamma triples. Singleton std=0 is not proof
of zero experimental noise; the GP retains its noise floor. Recorded duration comes
from the logged controller timeline clipped to the planned window, not host or
physical actuation time. Unknown durations remain NaN with a known-duration count.
Final selection uses repeated observations and reserved confirmation, not a lowest
single noisy score.

## Pure offline API and validation

```matlab
cfg = pacTuneOptimizerConfig();
% Supply explicit positive cfg.gammaBounds and cfg.baselineGamma.
observations = pacTuneOptimizerObservations();
p = pacTuneOptimizerSuggest(cfg,observations); % optional Nx3 pending matrix
validatePacTuneOptimizer();
validatePacTuneOptimizerAdapter();
validatePacTuneManualCampaign();
validatePacTuneEvaluator();
pacTuneDemoOptimizer('path/to/new_synthetic_report');
```

Observation rows carry gamma, objective, logical usableObjective/constraintKnown/
feasible flags, candidate ID and duration. Unknown objective uses NaN and
usableObjective=false; no callback is invoked by the core.

Reviewed demo: reports/optimizer_synthetic_v1. It contains 30 synthetic observations,
explicit repeats, seeded noise, a known feasible region, one simulated fault and
one missing-data case. Bounds/gamma/scores are software-test examples, not arm
recommendations. Reports include raw observations, repeat statistics, proposal MAT
state and plots. No prospective improvement on the arm is claimed.

Tests cover bounds/log mapping, exact baseline, pending avoidance, RNG isolation,
uncertainty/repeats, failure separation, deterministic restart, budget handling,
queue provenance, stale snapshots, confirmation isolation and existing regressions.
Temporary test campaigns use explicit test-only metadata; source MAT files remain
unchanged. A historical configuration-incomplete campaign returns blocked.

## Next step

Milestone D requires an agreed reproducible baseline, actual payload/calibration,
gamma bounds, fixed feedback/allocation settings, targets, budget and screening/
confirmation durations. Record resets/preconditioning and review the independent
operator stop procedure. The 2 s command-slew bottleneck and 4 s delayed roughness
remain relevant; optimizing gamma alone may not resolve them. No real hardware
campaign is configured or executed by this milestone.
