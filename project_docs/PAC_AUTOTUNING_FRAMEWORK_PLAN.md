# PAC experiment-based automatic tuning: implementation and handoff plan

Created: 2026-10-26. Updated after subsequent user authorization to develop the framework.
Status: milestones A-C implemented and validated on branch
`codex/pac-autotuning-framework`, starting from 63f679a. See
`passivity_adaptive_control/autotuning/README.md` and `MANUAL_WORKFLOW.md`.
Milestone A is committed as e6fb02c; milestones B/C are included in the subsequent
user-requested framework checkpoint (see git log).
Milestones D-E remain planned;
no hardware automation has been implemented or authorized by this document.

## 1. Purpose and scope

Develop a reproducible framework that proposes PAC tuning values, evaluates recorded experiments, and selects the next candidate. The eventual objective is accurate, smooth phi=1 rad circle tracking down to a 2 s period. Feasibility at that speed has not been established. Begin with an offline evaluator and a manual experiment loop. Keep the current PAC as the inner controller; optimize tuning between experiments.

The user proposed running two circles, collecting data, tuning, and repeating until an error bound is reached. Use two circles as an initial screening concept, not evidence of convergence or sustained stability. Trials need both a circle count and a duration in seconds, and promising candidates need longer confirmation. At 2 s/circle, two circles supply only 4 s of tracking data.

This plan does not implement NMPC, actuator dynamics, a residual model, integral action, or a new adaptation law. Those remain separate proposals. Failure to achieve the target with gain optimization should trigger model/actuator investigation, rather than unlimited gain increases.

## 2. Instructions for the next chat

1. Read `project_docs/PROJECT_STATE.md`, `passivity_adaptive_control/TUNABLE_BLOCK.md`, this document, and the relevant sections of `project_docs/VALIDATION_RESULTS.md`.
2. Inspect Git status and the actual MATLAB sources. Preserve user changes and unrelated files, especially `tmp/paper_draft`. Do not assume a clean checkout or settings from a previous conversation.
3. Inspect recorded inputs before interpreting a trial. Filenames containing `10s` include experiments with actual 2 s and 4 s periods. Current source defaults need not match the live MATLAB Function block.
4. Implement milestone A first. Work offline until the evaluator has passed its acceptance checks. Do not jump directly to an optimizer connected to hardware.
5. Do not edit SLX files unless explicitly requested. Provide paste-ready MATLAB Function source and instructions when integration is needed. Never run hardware models autonomously. Commit/push only when asked.
6. Update PROJECT_STATE and the appropriate specialized log after important findings or code changes. Report assumptions and missing information explicitly.

All paths in this document are relative to the repository root:
`C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl`.

## 3. Verified implementation context

- Mechanical coordinates are q=[l2;l3], with l1=-l2-l3. Three pressure commands are actuator inputs, not three mechanical coordinates.
- Current wrapper: `passivity_adaptive_control/pacAdaptiveTunableBlock.m`. Source tuning function: `pacTuningBlock.m`; validation/mapping: `pacTuningConfig.m`.
- Parameter order: `[k1;k2;k3;d2;d3;h1;h2;h3;b2;b3]`. Lambda and Ks are separate feedback gains. Adaptive damping coefficients d2/d3 multiply the model velocity term; they are not interchangeable with Ks.
- Intended PAC period is 0.05 s. `adaptationOptions(3)=0.075` is the maximum accepted update interval, not the execution period.
- Update before projection: `thetaDot = -(parameterScale.^2).*adaptationGamma.*(Y.'*s)`. Per-parameter rate limits, projection bounds and eligibility gates then apply. Gamma cannot overcome a closed gate or a binding rate/parameter limit.
- `freezeHysteresisWeights=true` freezes h1:h3 only. History states continue evolving; hysteresis force remains active. Keep this switch fixed within each initial tuning campaign.
- Common-pressure adjustment is bounded nominal allocation. Nominal allocation residual is not measured physical force error. Chamber pressure sensing and an actuator-lag model are absent.
- Negative spring estimates are experimental compensation coefficients. Current helper permits minimum lower bounds `[-100;-500;-100]` for k1:k3. The disk tuning source still supplies `[50;-100;0]` for those entries at plan creation; the user's live block may differ.
- Signed box bounds do not guarantee SPD for K=[k1+k3,k3;k3,k2+k3]. Log its minimum eigenvalue, but do not label an indefinite estimate alone as proof of closed-loop instability. There is no new stability guarantee for the signed configuration.

## 4. Evidence that should shape the design

Recordings are under `passivity_adaptive_control/experimentations/adaptCommonP/`. Existing audits in `tmp/` are references to inspect and refactor, not a stable framework API.

| Recording | Relevant finding | Design implication |
|---|---|---|
| `1rad_10s_gainTesting_v8.mat` | Late tracking plateaus without persistent final bound/rate saturation | Nonzero error does not automatically imply insufficient gamma |
| `1rad_10s_gainTesting_v9.mat` | Actual 2 s circles; about 89.25% samples pressure-limited and 5.42% adaptation-active; allocator replay attributed limits to command slew at 2 bar/s | Report adaptation opportunity and actuator-command constraints before proposing higher gamma |
| `1rad_10s_gainTesting_v10.mat` | Actual 4 s circles; roughness grows around 55 s while mean error improves; later gain changes do not reset estimates | Include sustained confirmation and roughness constraints; retain adaptive-state history |
| Same v10 | Tip-invalid grace fault at 104.85 s is distinct from earlier roughness | Separate fault, validity and tracking classifications; do not score unloading as tracking |
| v5/v6/v7 series | Accuracy, smoothness and parameter growth change with feedback and adaptation settings | Optimize a small subset with other settings fixed; do not infer causality from live changes |

Useful starting scripts: `tmp/auditPacGainTestingV8.m`, `tmp/auditPacV8Limits.m`, `tmp/auditPacGainTestingV9.m`, `tmp/auditPacV9Allocation.m`, `tmp/auditPacGainTestingV10.m`, and `tmp/auditPacV10Onset.m`. Check availability and saved reports first. `analyzePacAdaptive()` uses `out`; `analyzePacStaticSweep()` is for static holds, not circles.

## 5. Architecture and boundaries

Data flow:

`campaign definition -> candidate proposal -> manual trial execution -> immutable MAT + manifest -> evaluator -> acceptance classification -> optimizer observation -> next proposal`

Implement five independent components:

1. **Campaign definition:** fixed plant/trial settings, search variables, metric definitions, thresholds, budget and reproducible random seed.
2. **Candidate generator:** initially a supplied candidate table; later bounded Bayesian optimization. It has no hardware access.
3. **Trial intake:** validates manifest, signal provenance, logged tuning and completeness before evaluation.
4. **Evaluator:** deterministic offline metrics, per-cycle trends, plots, validity/constraint diagnostics, and classification.
5. **Campaign ledger:** append-only record of proposals, executions, failures, artifacts and decisions; supports restart without losing observations or rerunning trials accidentally.

No file watcher, model launch, device call or automatic gain write in the first implementation. A future execution adapter must be a separate component with explicit user authorization and verified stop behavior.

## 6. Search variables and fixed settings

First search vector: `x=[gammaK,gammaD,gammaH]`, mapped as:

```matlab
adaptationGamma = [gammaK;gammaK;gammaK;gammaD;gammaD; ...
                   gammaH;gammaH;gammaH;gammaB2;gammaB3];
```

Hold bias gains fixed. Use log-space search for strictly positive gamma values; represent zero/frozen adaptation as an explicit baseline or separate categorical choice, not log(0). Check MATLAB toolbox availability before selecting an optimizer implementation. If Bayesian Optimization support is unavailable, start with a bounded candidate table and retain the same evaluator/ledger interfaces.

Do not specify numerical search bounds until the current baseline and operator-approved range are recorded. Historic gamma values are evidence, not automatically acceptable bounds. Start locally around a reproducible smooth baseline. Do not let the optimizer widen bounds on its own.

Keep these fixed within a campaign: lambda, Ks, scales, parameter rates/bounds/initial values, history settings, bias gamma, enable/freeze switches, pressure rails/activation/slew/common-pressure mode, trim/calibration, supervisor settings, trajectory, payload and controller period. A material change creates a new campaign version. If later necessary, tune feedback gains in a separate stage; do not start with every PAC input as a free variable.

## 7. Repeatable trial protocol

Before the first trial, obtain a complete baseline configuration from actual inputs. Confirm payload mass and attachment timing relative to trim/calibration, model-parameter changes, direction, phi, period, and the user's acceptable operating envelope.

For each candidate:

1. Record trial/campaign/candidate IDs and configuration hash. Confirm the candidate is within the approved search box.
2. Apply settings before a fresh run. Reinitialize adaptive estimates, history and supervisor state using the same procedure. Verify from logs; changing gamma alone does not reset learned estimates.
3. Use a consistent physical starting condition and record available settling/preconditioning observations. A software reset does not reset pneumatic/material history or temperature.
4. Retain the existing startup/trim/feedback blend sequence. Use explicit phase labels when available; otherwise derive windows from verified logs and mark that method.
5. Hold candidate tuning constant during the evaluated run. Log all tuning outputs and switches. Classify live edits as exploratory segments, not clean candidate observations.
6. Save the original MAT and manifest without overwriting earlier runs. Evaluate only after acquisition ends in the initial manual workflow.
7. Produce a concise result and candidate decision. The user performs any next hardware trial.

Use screening and confirmation tiers. Exact durations require agreement before running. Confirmation must extend beyond the previously observed roughly 55 s deterioration region, with a suitable margin selected by the operator; two circles cannot substitute for this check. Compare baseline/candidate repeats in an interleaved order to expose drift or changing physical conditions. Use identical initialization and duration when comparing candidates.

## 8. Data contract and manifest

Each trial should carry a JSON/MAT manifest with:

- Schema version, campaign/candidate/trial IDs, creation time, source MAT path and file hash.
- Git revision plus dirty status and an archived tuning snapshot; revision alone does not describe uncommitted/live block settings.
- Candidate values and all fixed tuning vectors, lambda/Ks, switches and controller/sample periods.
- Desired trajectory, direction, phi, cycle period/count, startup/learning/scoring/confirmation intervals.
- Payload, attachment/calibration timing, calibration/trim identifiers, nominal model version and operator notes. Missing fields remain unknown, not inferred.
- Initialization/reset procedure, physical preconditioning and any operator intervention.
- Metric schema, normalization constants, thresholds and campaign budget version.

Intake must inventory Dataset elements by index, name, block path, dimensions, units and time base. Repeated names exist; choose by provenance and store the resolution. Never silently assume the first match is correct. Handle timeseries `IsTimeFirst=false` and matrix channels explicitly. Keep native timestamps and document alignment; do not interpolate through invalid gaps or extrapolate nonconstant signals.

Required signal groups, subject to verifying exact names/provenance:

- Reference/measured q, desired/measured velocity, desired theta/phi and measured theta/phi.
- Validity/freshness flags, controller status/fault, adaptation activity, limiting flags, nominal allocation residual and final pressure commands.
- Parameter estimates, history state, all tuning vectors, feedback gains and enable/freeze flags.
- Desired/pre-slew pressure and force diagnostics when available; optional compensated XYZ with a documented coordinate convention.

Common existing names include `sampled des_length`, `meas_length`, `filtered_d_length`, `d_length`, `theta`, `phi`, `des_theta`, `des_phi`, `parameterEstimate`, `historyState`, `lambda_gains`, and `damping_gains`. Some recordings use `hystoryOptions` (typo). Resolve semantics from actual sources. Do not invent missing measurements.

## 9. Evaluator definitions

Version metric definitions so campaign results remain comparable. Use fixed normalizers for the whole campaign, not candidate-specific min/max normalization.

**Windows and masks:** report startup/blend, learning, steady scoring and confirmation separately. Tracking metrics require valid/fresh measurements and the appropriate normal tracking state (verify status definitions; currently 3 normal, 4 invalid hold, 5 latched fault). Report invalid/held/fault exposure separately. Do not allow a candidate with fewer valid samples to win by excluding difficult portions: require minimum coverage and complete planned duration. After a fault, mark the trial failed and preserve the valid pre-fault diagnostics.

**Tracking:** per-coordinate q RMSE in mm; vector q RMSE `sqrt(mean(e2.^2+e3.^2))`; wrapped theta error `atan2(sin(thetaDes-theta),cos(thetaDes-theta))`; phi MAE/RMSE; per-lap and full-window values. Theta requires observability masking. Also report maximum/upper-quantile errors and coverage. Optional XYZ RMSE uses the same compensated sensor frame and geometry for desired and measured coordinates; no fitted alignment that conceals error.

**Smoothness:** desired-minus-measured velocity RMSE; cycle harmonic residual or another explicitly defined repeatable roughness metric; theta-rate variability and near-stop fraction as secondary diagnostics only where observable and timestamps valid. Derive angles across valid consecutive samples with unwrap/gap handling; differentiate with a documented fixed filter. Do not interpret noisy second differences as physical acceleration. Compare against the desired motion so faster trajectories are not penalized merely for moving faster.

**Adaptation:** active fraction, eligible fraction where reconstructible, gate reasons, per-parameter bound occupancy, actual and reconstructed unclipped update rates, rate-cap occupancy, per-cycle estimate changes, and force-component RMS. Respect output-before-update versus post-update flag timing. If Y or gates cannot be reconstructed exactly, label rate/gate attribution unavailable or approximate.

**Constraints/data quality:** commanded pressure range, command slew, limiting occupancy, nominal residual distribution, invalid runs, measurement age, timing gaps and faults. Distinguish pressure rails, slew and allocation infeasibility using available evidence. Diagnostic forces may return zero on invalid samples; mask these before interpreting parameter resets or physical force loss. No measured chamber pressure or physical force error is available.

**Model diagnostics:** estimated K minimum eigenvalue/determinant, damping and hysteresis weight growth, history bounds, bias evolution. Record these as diagnostics or agreed empirical constraints; do not equate physical parameter identification with adaptive compensation.

## 10. Objective, constraints and candidate classification

Initial objective proposal: minimize normalized q tracking RMSE over a fixed scoring window. Record learning-curve area and time to a sustained target as secondary metrics. This keeps the objective interpretable while the user's desired error thresholds remain undecided.

Use separate acceptance constraints for theta/phi error, smoothness, coverage, pressure limiting and sustained confirmation. Specify all numeric thresholds and units in campaign configuration before optimization; unset required thresholds must prevent a campaign from claiming success. Do not invent an acceptable theta error or roughness limit from this plan.

Classification should distinguish:

- `valid_feasible`: complete comparable trial satisfying configured constraints.
- `valid_infeasible`: usable observation with an operating/performance constraint violation.
- `invalid_data`: missing/ambiguous required signals or inadequate coverage; not a numeric success.
- `protocol_deviation`: live gains, changed payload/settings, incomplete/reset mismatch; exploratory only.
- `aborted_or_faulted`: retain event/censoring time and preceding metrics; cannot win.

When integrating the optimizer, use supported feasibility constraints or an explicitly documented failure model. Do not encode missing data as zero error, ignore failures, or replace every failure with an arbitrary finite score without recording the reason. A bounded optimizer and empirical rejection rules do not establish closed-loop stability.

## 11. Proposed MATLAB files and artifact layout

Proposed location: `passivity_adaptive_control/autotuning/` (not created by this plan).

| File/component | Responsibility |
|---|---|
| `pacTuneCampaignConfig.m` | Versioned fixed settings, variable mapping, bounds, metrics, constraints, budget |
| `pacTuneReadTrial.m` | Load original MAT, resolve provenance, normalize shapes and timestamps |
| `pacTuneValidateTrial.m` | Compare logged settings to manifest, verify completeness and phase windows |
| `pacTuneEvaluateTrial.m` | Deterministic metrics, classifications and per-lap table |
| `pacTuneReportTrial.m` | Human-readable report and scientific plots |
| `pacTuneProposeCandidate.m` | Manual-table/baseline/local-search interface, later Bayesian proposals |
| `pacTuneRecordObservation.m` | Append result and artifact references; guard duplicate IDs/hash |
| `pacTuneCampaignStatus.m` | Best confirmed candidate, pending work, budget and stop reason |
| `README.md` | Exact offline/manual workflow and example calls |

Suggested campaign artifacts: immutable campaign config, `candidates.csv`, append-only observation ledger, per-trial manifests and reports, per-cycle metrics, optimizer state/seed and a final recommendation report. Reference original MAT files rather than duplicating large recordings unnecessarily. Write outputs outside source MAT locations or into dedicated analysis directories; never overwrite raw data.

## 12. Milestones and acceptance checks

### A. Offline evaluator and historical replay

Implement intake, masks, metrics, classifications and reporting. Refactor useful audit logic into reusable functions while preserving original scripts. Replay v8, v9 and v10; optionally v5/v7 for accuracy/smoothness tradeoffs. Historical live gain changes must be detected and segmented, not treated as one stationary candidate.

Acceptance: reproduce existing reported metrics on explicitly matched windows within declared numerical tolerances; recognize actual 2 s/4 s periods; expose v9 slew limitation and low adaptation activity; expose v10 pre-fault roughness and later separate fault. Show coverage next to every score. These checks establish evaluator consistency, not causal proof or hardware readiness.

### B. Manual candidate loop

Create configuration/manifest templates, candidate table, ledger and baseline report. Generate a human-reviewable proposed tuning vector without writing to a model. Import a manually acquired trial, validate actual settings, and produce a next-step report. Resume a campaign from disk and prevent duplicate ingestion.

Acceptance: a complete offline dry run from candidate creation through import/evaluation/status works on existing data or clearly labelled synthetic fixtures. Missing fields and deviations produce explicit classifications. No hardware/model launch occurs.

### C. Bounded optimization

After agreeing baseline, bounds, thresholds and budget, add the optimizer adapter. Treat scores as noisy; retain replicates and use a reserved confirmation set rather than selecting solely on the best single trial. Seed reproducibly; include the baseline. Start with the three grouped gamma variables. Report uncertainty and the number/duration of evaluations.

Acceptance: test optimizer bookkeeping on a known noisy synthetic objective with feasible/infeasible regions, validate bound mapping/log transforms and failure handling, and replay historical observations without claiming prospective optimization success. Confirm toolbox dependencies and a restartable state.

### D. Operator-run experimental campaign

Start at an agreed reproducible period with adequate command headroom, then progress toward 4 s and 2 s in separate campaign versions. No automatic extrapolation of performance across periods. Repeat baseline/candidate trials, inspect confirmation roughness and constraints, and keep the best confirmed configuration available.

Acceptance: repeated confirmation meets all configured targets for the specified payload, direction, duration and period. Report unsuccessful outcomes honestly. Do not declare success from falling MAE alone or from a censored short run.

### E. Optional execution automation

Separate future request. Before implementation, agree how the user authorizes runs, where commands are applied, independent operator stop behavior, logging/heartbeat failures, reset verification and recovery. Offline constraints cannot replace the existing controller/hardware supervision. Automatic restarts after a fault require a separately agreed policy; do not add them implicitly.

## 13. Tests and review

Use meaningful fixtures for duplicate names/provenance, time-last data, matrix estimates, angle wrap, unobservable theta, invalid gaps, fault censoring, changed tuning, constant Inf stop time, bound/update timing and incomplete trials. Compare Inf-containing configuration values with appropriate equality rather than `diff` alone (`Inf-Inf` becomes NaN).

Test that a faulted, rough or poorly covered trial cannot beat a valid candidate simply by masking its bad samples. Test consistent score normalization, unit conversion and window selection. Use real recording regressions for metric definitions and synthetic tests for parser edge cases/optimizer bookkeeping. Run only offline MATLAB tests; no Simulink/hardware execution is necessary for milestones A-C.

## 14. Decisions needed before the first new tuning campaign

1. Verified baseline: all live tuning inputs, feedback gains, signed stiffness bounds, payload/trim and initial estimates.
2. Initial period/direction/phi and progression criteria toward 2 s.
3. Screening and confirmation durations in seconds and circles; initialization/preconditioning procedure and repeat count.
4. Acceptable q/theta/phi errors, roughness, coverage and pressure-limiting thresholds; operator abort criteria.
5. Approved gamma bounds, evaluation/time budget, and available MATLAB optimization toolboxes.
6. Whether the task prioritizes learning rapidly from a fresh reset or final tracking after a long learning phase. Use one protocol per campaign; report both learning and final performance.

These answers do not block building the offline evaluator. They do block claiming a hardware campaign is configured or successful. Do not ask the user to repeat information that can be recovered from the actual logs and sources.

## 15. Immediate next action

Milestone A is now implemented: native-signal inventory/provenance, bounded alignment,
explicit constant-switch lifetime policy, metrics and classifications, adaptation
gate/rate diagnostics, nominal allocation replay, reports and regression tests.
Reviewed replay reproduces v8 per-lap errors, v9's 1071 slew mismatches and 65/1200
active samples, and v10's pre-fault roughness increase plus the 104.85 s fault.
Read the implementation README for current limitations (including constant-speed
circle assumptions and no harmonic/XYZ score). No hardware execution.

Milestone B is implemented (2026-10-27): immutable campaign snapshots, bounded
operator-supplied gamma proposals, full tuning review artifacts, manifests,
offline trial intake, append-only event files, immutable-config hash checks,
duplicate ID/file checks, evaluation budgets, restart/status and repeat-based
screening/confirmation selection. Longer confirmation windows are supported;
failed confirmation prevents or withdraws confirmation. See MANUAL_WORKFLOW.md.
The saved v9 demo remains unconfigured/ineligible, not a tuning result.

Milestone C is implemented (2026-10-27): bounded proposal-only GP search in
normalized log gamma space; seeded baseline/initial design; repeated/noisy scores;
explicit usable-objective and known-constraint flags; held-out confirmation;
ledger reconstruction and stale-suggestion rejection. See OPTIMIZER_WORKFLOW.md.
The 30-observation demo is synthetic only. Known comparable faults contribute
constraint evidence without partial scores; protocol/data failures stay unknown.
Fixed-kernel GP uncertainty and kernel-weighted feasibility are model-dependent
diagnostics, not physical safety guarantees or evidence of arm improvement.

Next: milestone D, a configured operator-run experimental campaign. Thresholds,
search bounds, physical reset/preconditioning and confirmation duration remain
unsettled. Do not invent real-arm settings, treat exploratory historical trials as
clean candidates, or claim successful physical tuning. Hardware execution remains
manual. Do not run hardware or commit/push without the relevant user request.
