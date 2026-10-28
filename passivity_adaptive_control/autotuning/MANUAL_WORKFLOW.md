# Manual PAC tuning campaign workflow

Implemented: 2026-10-27, milestone B. This workflow archives proposals and evaluates
saved trials. It never applies gains, runs a model, starts a device, or chooses
an optimized next candidate through its manual API. Milestone C now adds a
proposal-only optimizer adapter; see OPTIMIZER_WORKFLOW.md.

## Files and entry points

- `pacTuneManualConfig`: campaign template with unknown settings left unset.
- `pacTuneBaselineFromTrial`: extract all fixed inputs from a recorded snapshot.
  This is evidence of what was logged, not a recommended or verified baseline.
- `pacTuneManualCampaign`: actions `create`, `propose`, `import`, `status`.
- `pacTuneTrialManifest`: operator-completed metadata for one trial.
- `pacTuneReportCampaign`: export a fresh human-readable status and CSV tables.
- `validatePacTuneManualCampaign`: synthetic offline end-to-end checks.
- `pacTuneDemoManualCampaign`: existing-v9 replay demonstration; no new trial.

Use these from MATLAB with both `passivity_adaptive_control` and its `autotuning`
subdirectory on the path. Read the evaluator README for signal and metric semantics.

## Configure an immutable campaign

```matlab
addpath('passivity_adaptive_control','passivity_adaptive_control/autotuning');
C = pacTuneManualConfig();
C.id = 'my_campaign_v1';
```

Complete these fields before a real manually executed campaign:

| Field | Meaning |
|---|---|
| `purpose` | `manual` for real operator-run observations; `draft`/`offline_demo` cannot qualify candidates |
| `baselineTuning` | All fifteen fields returned by `pacTuneBaselineFromTrial`, reviewed against actual live inputs |
| `gammaBounds` | 3x2 positive finite limits for gammaK/gammaD/gammaH; bounds are never widened automatically |
| `evaluationBudget` | Positive integer number of imported observations, including failures |
| `confirmationRepeats` | At least two independent accepted confirmation observations |
| `minimumConfirmationDuration_s` | Required duration of each confirmation, chosen to expose delayed deterioration |
| `evaluator` | Fixed windows, period, normalization, coverage and required length/smoothness/limiting targets; angular limits optional in pac-tune-2 |
| `confirmationEvaluator` | Optional evaluator copy with longer tracking/scoring windows only; all other settings must match |
| `conditions` | Payload mass/attachment timing, calibration/model IDs, reset/preconditioning procedures, phi/direction and code revision |

For example, extract an existing snapshot for inspection:

```matlab
trial = pacTuneReadTrial('path/to/baseline.mat',C.evaluator);
C.baselineTuning = pacTuneBaselineFromTrial(trial,20);
```

The snapshot includes gamma, scales, bounds, rates, initial estimates, history,
pressure, supervisor, feedback and switch inputs. Only the grouped gamma entries
1:8 change in a candidate. Bias gamma remains fixed. Compare snapshots to live
block code: a saved source function and a live block can differ.

Targets default to NaN, and experimental details default to unknown. Do not fill
optional theta/phi/vector-RMSE limits just to make the campaign ready: only the
required limits in pacTuneTargetPolicy need values. Old pac-tune-1 snapshots retain
their original policy; adopt the new length-based policy in a new campaign.
Do not invent values for required limits or experimental conditions. Drafts can be archived, and
bounded proposals can be reviewed when baseline/bounds/budget are present, but
unconfigured campaigns cannot produce accepted or confirmed candidates.

Create a new campaign directory after completing the intended snapshot:

```matlab
state = pacTuneManualCampaign('create','path/to/my_campaign_v1',C);
```

`campaign.mat` is immutable and its SHA-256 is checked against the creation event
on every subsequent operation. If settings need revision, create a new campaign
directory/version. Do not edit this MAT in place. A configuration-ready flag means
the software fields are populated; it is not a stability or hardware safety claim.

## Propose a candidate for manual review

```matlab
candidate = pacTuneManualCampaign('propose','path/to/my_campaign_v1', ...
    'candidate01',[gammaK;gammaD;gammaH]);
```

Gamma values are operator inputs, not automatically generated recommendations.
The returned full tuning snapshot and `candidates/candidate01/review.md` contain
paste-ready vectors. `candidate.mat` retains the exact snapshot. Review files map
the logged `hystoryOptions` spelling to the live `historyOptions` input, and
`lambda_gains`/`damping_gains` to lambda/Ks. No gains are applied by this command.

IDs use 1-64 letters, digits, hyphens or underscores; paths/traversal are rejected.
Duplicate candidate IDs and out-of-bounds proposals are rejected. A fixed finite
evaluation budget is required for proposals/imports. Trial budget counts raw
imports, not only successful observations.

## Import an operator-run recording

After the user has executed and saved a trial:

```matlab
m = pacTuneTrialManifest();
m.kind = 'screening'; % later use 'confirmation'
m.conditions = C.conditions; % only when actual conditions match
m.resetVerified = true;      % only after the operator verifies the procedure
m.operator = 'operator identifier';
m.notes = 'Actual observations/interventions';
observation = pacTuneManualCampaign('import','path/to/my_campaign_v1', ...
    'trial01','candidate01','path/to/recording.mat',m);
```

The importer checks logged tuning against the full proposed snapshot, initial
parameter estimates/history, reference phi/direction/period, actual conditions,
operator/reset attestation, coverage, faults and target limits. It hashes the
original MAT and rejects an already imported byte-identical file even under a
new filename or trial ID. It cannot detect every semantically duplicated re-save,
prove physical reset, measure payload attachment, or establish trial independence.
Those remain operator responsibilities, recorded in the manifest.

Each trial has original source path/hash, manifest/candidate MAT, evaluator MAT,
report and diagnostic CSV/PNG files. The original recording is not modified or
duplicated. The report distinguishes evaluator criteria from final campaign
eligibility; a low numerical error does not override a metadata mismatch.

Malformed schemas/missing required manifest fields are rejected before ingestion.
Represent genuinely unknown conditions explicitly and retain protocol-deviation
observations; do not mark reset verified merely to bypass rejection.

## Screening, confirmation and next action

```matlab
state = pacTuneManualCampaign('status','path/to/my_campaign_v1');
pacTuneReportCampaign('path/to/my_campaign_v1','path/to/new_status_report');
```

Status returns candidate and observation tables, remaining budget, configuration
issues, best screened/confirmed candidate, stop reason and next action. Reports
require a new output directory so earlier snapshots are retained.

A candidate needs accepted screening and the configured number of accepted
confirmation observations, each long enough for the configured minimum. Screening
and confirmation may use different windows, but share metric definitions,
normalizers, fixed tuning and conditions. Ranking uses mean objective within each
tier; it does not compare screening-window scores against confirmation scores.

Every screening/confirmation attempt remains in the ledger. A failed screening
prevents that candidate being selected; a failed/short confirmation prevents or
withdraws confirmation. Passing trials cannot be cherry-picked around failures.
Review such a candidate and use a new candidate/campaign version when justified.
No recommendation is made to rerun automatically after a fault.

Campaign success means the configured empirical criteria were met on the recorded
screening/repeats. It is not a formal stability claim or an automatic hardware
authorization. Hardware-independent operator review remains necessary.

## Restart and recovery

The authoritative append-only ledger is `events/000001.mat`, `000002.mat`, etc.
Each event is written to a temporary file then renamed. Status reconstructs from
these files after a MATLAB restart; CSVs are exported snapshots, not mutable state.
A `.writer.lock` prevents simultaneous writers. Read-only status does not modify
the campaign. Configuration hashing detects configuration edits; the event store
is an audit convention, not a cryptographically signed tamper-proof database.

A process failure may leave a lock, `.pending` event, or a trial/candidate artifact
without its event. The importer refuses to overwrite orphan artifacts. Inspect the
campaign and preserve evidence before manually recovering; do not blindly delete
locks or regenerate events. Automatic crash repair is not implemented.

## Offline demo and tests

```matlab
validatePacTuneManualCampaign();
validatePacTuneEvaluator();
pacTuneDemoManualCampaign('path/to/new_offline_demo');
```

The persistent reviewed demo is `reports/manual_demo_v1`. Its v9 recording is
historical, its bounds are illustrative software-test values, payload/reset data
remain unverified, and it cannot yield success. Synthetic test campaigns use
temporary files and permissive test-only targets; they are not experiments.

Milestone C now supplies bounded GP suggestions through pacTuneSuggestNext and
pacTuneQueueSuggestion. Actual bounds, baseline, thresholds, budget and confirmation
durations still need agreement before a real campaign. Hardware execution remains
manual; no hardware adapter has been implemented.
