# PAC offline tuning evaluator

Milestone A of `project_docs/PAC_AUTOTUNING_FRAMEWORK_PLAN.md` is implemented.
These functions read saved MAT data; they do not open/run models, operate devices,
change gains or select the next candidate. No optimization toolbox is required.
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

Targets are intentionally NaN. Enter agreed finite nonnegative limits for every
target before evaluating feasibility. Unset targets produce `unconfigured_targets`
on otherwise valid fixed-tuning data. Faults, incomplete runs, low coverage and
live tuning edits take precedence. `cfg.expectedTuning` can supply exact expected
vectors for mismatch checks; archived candidates/manifests are not implemented yet.

`result.success` means this one offline observation met the supplied evaluator
criteria. It does not establish repeated confirmation, physical stability, a
verified reset, or an approved hardware campaign. Historical data here never
claim success. Smoothness currently includes velocity error and angular-rate
diagnostics; cycle harmonic residual, time-to-target, uncertainty estimates,
XYZ scoring, candidate generation, ledger, manifests and Bayesian optimization
remain future work. There is no automatic hardware execution component.

## Regression evidence

`validatePacTuneEvaluator` checks duplicate provenance, time-last matrix signals,
single-point/gapped alignment, unset targets, masking attacks, constant Inf
settings, matched v8 lap errors, v9's exact 1071 slew mismatches and 65 active
samples among 1200, and v10's falling MAE alongside rising pre-fault rate SD plus
the separate 104.85 s fault. Synthetic fixtures are created only in a temporary
directory and removed after testing. Raw experiment files are not rewritten.

The reviewed report set is under `reports/historical_reviewed_v1/`. Next work is
milestone B: a manual candidate/manifest/ledger loop, using agreed campaign
settings. Keep optimizer integration and hardware execution separate.
