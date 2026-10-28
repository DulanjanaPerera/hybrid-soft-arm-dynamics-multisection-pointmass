# PAC offline tuning evaluator

For the implemented continuous Simulink tuning blocks, separate-process BO
worker and new online model copy, see `ONLINE_SIMULINK_WORKFLOW.md`. Only
synthetic runtime/library/transport tests have been run; hardware integration
verification is pending. The original offline/manual workflow remains available.

The earlier host-only prototype (`pacTuneOnlineStep`, `ONLINE_WORKFLOW.md`)
is retained for history; the new runtime integration does not use its host timer.
The descriptions below cover the offline evaluator/proposal APIs.

Milestones A-C of `project_docs/PAC_AUTOTUNING_FRAMEWORK_PLAN.md` are implemented.
For the manual candidate/manifest/ledger loop, read `MANUAL_WORKFLOW.md`.
For bounded proposal-only optimization, read `OPTIMIZER_WORKFLOW.md`.
The evaluator reads saved MAT data; the optimizer only proposes candidates.
Neither opens/runs models, operates devices or changes live gains. No optimization
toolbox is required.
MATLAB and the Simulink data classes are required to load these recorded logs.

## Run the historical replay

From the repository root in MATLAB:

```matlab
addpath('passivity_adaptive_control/autotuning');
validatePacTuneEvaluator();
results = pacTuneReplayHistorical( ...
    'passivity_adaptive_control/autotuning/reports/my_review');
```

Choose a new output directory; the reporter refuses to overwrite existing reports.
Replay includes v8 (10 s), v9 (2 s), and v10 (4 s) despite their filenames. v10's
planned tracking window ends at 140 s, but the recording ends early and faults.
It is classified `aborted_or_faulted`, even when the selected pre-fault scoring
window has good coverage. v8/v9 both contain live gamma edits. v8 is classified
`protocol_deviation`; v9 is `invalid_data` at the default 95% coverage requirement
(its tracking coverage is about 94.83%). Both remain exploratory observations.

## Evaluate another recording

```matlab
cfg = pacTuneCampaignConfig();
cfg.trackingWindow_s = [20 140];
cfg.scoringWindow_s = [94 102];
cfg.expectedPeriod_s = 4;
trial = pacTuneReadTrial('path/to/recording.mat',cfg);
result = pacTuneEvaluateTrial(trial,cfg);
pacTuneReportTrial(trial,result,'path/to/new_report_directory');
```

Set windows from the actual intended protocol, not the filename. A recording
ending before the planned window remains incomplete. Missing required signals
produce `invalid_data`; duplicate signal names require an explicit selector.
Inspect the Dataset directly to choose a block path/index, then set, for example:

```matlab
cfg.selectors.parameterEstimate = 'model/PAC'; % exact block path
% or cfg.selectors.parameterEstimate = 84;     % exact Dataset index
```

The chosen index/path, native dimensions, units and time spans are saved in the
inventory and in `trial.signals`. The original MAT receives a SHA-256 hash and
is never modified. The saved evaluation contains the configuration and metrics;
native raw arrays remain in the original MAT to avoid needless duplication.

## Alignment, coverage and metric definitions

- The parameter-estimate timeline anchors evaluation. Other signals use the
  previous native sample within 0.075 s, without extrapolation or gap filling.
  Native NaNs remain invalid. Time-last arrays are reshaped explicitly.
- Three declared constant switch inputs may be held beyond their final logged
  timestamp only when every native value is identical. This policy is explicit
  in `cfg.constantInputFields` and the report. Varying signals never extend.
- Tracking requires aligned finite signals, model/orientation/configuration/base
  validity, measurement freshness, and status 3. Theta also requires observability.
  Diagnostic zero forces on invalid samples do not enter force RMS.
- Coverage uses the planned window duration divided by the median controller
  sample interval. This version targets the existing nominally uniform PAC logs;
  irregular/event-based controller timelines need a time-weighted extension.
  No interpolation through missing measurements is used to improve coverage.
- q2/q3 RMSE and vector RMSE are in mm. The objective is scoring-window vector
  RMSE divided by the fixed `lengthNormalizer_mm`. Theta errors are wrapped;
  phi errors are ordinary differences. Velocity RMSE is the vector norm in m/s.
- Theta-rate SD uses consecutive valid angle increments. `nearStopFraction`
  uses absolute rate below 20% of the reference rate. `slowOrReverseFraction`
  also counts reverse motion and matches the older signed slow-rate diagnostic
  on these positive-direction circles. Neither is an independent speed sensor.
- Per-cycle windows use measured reference angular speed and the configured
  tracking start. Windows are half-open with numerical boundary tolerance.
  This assumes constant-speed circles, not arbitrary trajectories or reversals.

## Adaptation and allocation diagnostics

The evaluator reconstructs first-failing update gates using logged inputs and the
current wrapper's timestamp logic. It reports active/gate disagreement counts.
Observed parameter rates account for pre-update output timing. Bounds refer to
stored estimates; projected post-update flags can differ by one sample. Optional
history/s signals allow raw gradient-rate reconstruction using the current Y
structure; rate statistics explicitly distinguish caps from raw rates.

Allocation replay is performed only for valid normal-control samples with common
pressure enabled and the required forces available. It uses the existing nominal
allocator, logged pressure settings and logged trim force. A target/command
mismatch is attributed to slew only when replay agrees with the logged limiting
flag. Check `loggedFlagMismatch`; future allocator changes require revalidation.
This is command/nominal-model evidence, not measured pressure or force delivery.

## Acceptance and limitations

New configs use `pac-tune-2`. Required finite nonnegative limits are:
`q2Rmse_mm`, `q3Rmse_mm`, `velocityRmse_mPerS`,
`velocityErrorRateRmse_mPerS2`, and `limitedFraction` (0..1).
Vector length RMSE remains the optimization objective. Its acceptance limit and
theta/phi/theta-rate limits are optional: leave NaN to disable them. Numeric
values are intentionally unset; baseline measurements are not agreed targets.
Unset required limits produce `unconfigured_targets`
on otherwise valid fixed-tuning data. Faults, incomplete runs, low coverage and
live tuning edits take precedence. Angular observability/coverage gates acceptance
only when the corresponding angular limit is enabled; angles remain diagnostics.
Required length/velocity valid-pair coverage still gates length-based acceptance.
`velocityErrorRateRmse_mPerS2` is the vector RMS of consecutive velocity-error
differences divided by sample dt on valid pairs with no gaps. It is a filtered
length-based roughness proxy, not independently measured physical acceleration.
Keep filtering/sample period fixed when selecting/comparing that limit.
Existing `pac-tune-1` campaigns keep their original mandatory angular policy;
create a new immutable campaign to adopt version 2.
`cfg.expectedTuning` can supply exact expected
vectors for mismatch checks. The manual workflow archives candidates/manifests
and adds campaign-level metadata and repeated-confirmation checks.

`result.success` means this one offline observation met the supplied evaluator
criteria. It does not establish repeated confirmation, physical stability, a
verified reset, or an approved hardware campaign. Historical data here never
claim success. Smoothness currently includes velocity error and angular-rate
diagnostics; cycle harmonic residual, time-to-target, uncertainty estimates,
XYZ scoring remains future work. Bounded GP candidate generation is available
through the proposal-only optimizer; manual proposals, manifests and the ledger
remain the execution/intake interface.
There is no automatic hardware execution component.

## Regression evidence

`validatePacTuneEvaluator` checks duplicate provenance, time-last matrix signals,
single-point/gapped alignment, unset targets, masking attacks, constant Inf
settings, matched v8 lap errors, v9's exact 1071 slew mismatches and 65 active
samples among 1200, and v10's falling MAE alongside rising pre-fault rate SD plus
the separate 104.85 s fault. Synthetic fixtures are created only in a temporary
directory and removed after testing. Raw experiment files are not rewritten.

The reviewed evaluator reports are under `reports/historical_reviewed_v1/`.
The manual replay demo is under `reports/manual_demo_v1/`; the synthetic optimizer
demo is under `reports/optimizer_synthetic_v1/`. Next is milestone D, an explicitly
configured operator-run campaign; real bounds/targets/baseline remain unset.
Keep optimizer integration and hardware execution separate.

2026-10-28 smoothness-first extension: `pacTuneSmoothCampaignDraft` prepares a
new, non-active configuration with user-approved2mm coordinate ceilings and
roughness as the objective. `pacTuneObjectivePolicy` keeps old length objectives
unchanged when the new field is absent. Optional fixedRawLength motion scoring
uses the same offline differentiator across runs, so changing live filters cannot
artificially improve the optimization score. See the first unloaded experiment
folder's filter_review_v1/README.md for filter replay and manual smoke instructions.
