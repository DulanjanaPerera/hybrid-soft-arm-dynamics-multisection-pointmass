# GRU capacity comparison: completed offline experiment

2026-11-03. Nine architecture/seed cases: 1x64, 1x128, 2x64; seeds113/211/307.
Training gradients use October19/27 TRAIN only. November3 whole-trial holdouts participate in checkpoint and architecture selection.
November2 remains excluded from fitting/selection, but is previously inspected evaluation, not blind confirmation.

## Matched protocol

- Stage1: original streaming protocol, up to1200 Adam updates at3e-4; checkpoint by old holdouts.
- Stage2:1200 extra old-only updates at1e-4 with fresh Adam from stage1 selected weights.
- Same batch32, clip1, .5s measured-state TBPTT chunks/.5s auxiliary free path, .5s burn-in, losses, physics and integration.
- Unchanged old TRAIN normalizer, old/Nov3 splits and208 internal anchors. No interpolation or new training data.
- Stage2 checkpoint score=.5old+.5Nov3; each cohort=.7one-step+.3half-second path, equal recording means.
- Same10% old-regression guard relative to EACH case's stage1 checkpoint; starting checkpoint eligible.
- Architecture selected by three-seed mean composite; seed113 fixed in advance as representative, not luckiest seed.
- Exact existing64-unit stage1 fits and seed113 continuation reused by hash: six new stage1 and eight new continuation fits.
- Same update budget is not proof all architectures converged. Parameter counts differ; no parameter-matched depth claim.

## Internal architecture comparison

| Architecture | Parameters | Mean score (mm) | Sample SD | Seed113 | Seed211 | Seed307 |
|---|---:|---:|---:|---:|---:|---:|
| 1x64 | 14162 | 0.58531 | 0.02495 | 0.58716 | 0.55949 | 0.60928 |
| 1x128 | 52882 | 0.48684 | 0.02226 | 0.46154 | 0.49553 | 0.50344 |
| 2x64 | 38930 | 0.53814 | 0.05209 | 0.50814 | 0.50799 | 0.59828 |

Predeclared winner: **gru128**. Frozen representative candidate: candidate_force_model.mat.
SHA256 `8f2e0bb55ea9399b62abd5abe0b8f6dd15ab0470a400ad3e593e1edcd8f9f545`.
The existing fine-tuned GRU64 research release remains unchanged. This candidate is an offline result, not deployed or controller-qualified.

## All-seed short-path evaluation

Scores below are per-seed recording-macro scores, then means/SD across three seeds. Short paths use the same common anchors.

| Cohort | Architecture | One-step mean +/- SD (mm) | Half-second path mean +/- SD (mm) | Mean recording force-rate P95 (N/s) |
|---|---|---:|---:|---:|
| old_internal | 1x64 | 0.45492 +/- 0.01844 | 0.36129 +/- 0.02766 | 12.803 |
| old_internal | 1x128 | 0.34287 +/- 0.01913 | 0.26554 +/- 0.00473 | 10.825 |
| old_internal | 2x64 | 0.41233 +/- 0.04017 | 0.29697 +/- 0.02785 | 10.464 |
| nov3_internal | 1x64 | 0.75112 +/- 0.05558 | 0.72669 +/- 0.02898 | 18.570 |
| nov3_internal | 1x128 | 0.69187 +/- 0.03718 | 0.56564 +/- 0.03530 | 18.663 |
| nov3_internal | 2x64 | 0.74199 +/- 0.07239 | 0.59721 +/- 0.09731 | 15.260 |
| nov2_external | 1x64 | 0.33904 +/- 0.03756 | 0.27698 +/- 0.01405 | 10.628 |
| nov2_external | 1x128 | 0.26730 +/- 0.01547 | 0.20983 +/- 0.02488 | 9.883 |
| nov2_external | 2x64 | 0.34491 +/- 0.04681 | 0.23219 +/- 0.02527 | 10.373 |

## Whole-record one-step evaluation: fixed seed113

All eligible anchors, equal means across five recordings; length-vector RMSE, not single-coordinate error.

| Cohort | 1x64 | 1x128 | 2x64 |
|---|---:|---:|---:|
| nov2_external | 0.36143 | 0.28691 | 0.31188 |
| nov3_selection | 0.75485 | 0.64815 | 0.72412 |
| nov3_train | 1.04670 | 0.93190 | 1.04314 |

Nov3 TRAIN is diagnostic only; these weights were not optimized on Nov3 TRAIN, but it is not the declared whole-trial validation subset.

## Long recorded-input free predictions: fixed seed113

One observed initialization plus causal observed recurrent prefix; thereafter predicted q/v only. Future recorded pressure/time known.
Missing observed truth stays NaN/unscored. All Nov3 comparison windows are the existing held-out whole trials, including helix trial2.

| Day | Motion | Channel | Trial | 1x64 | 1x128 | 2x64 |
|---|---|---:|---:|---:|---:|---:|
| Nov2 | circleSweep | 0 | 3 | 0.62291 | 0.92961 | 0.56961 |
| Nov2 | SineWave4 | 1 | 1 | 1.43739 | 0.67309 | 1.06268 |
| Nov2 | SineWave4 | 2 | 3 | 1.30422 | 1.99627 | 1.54568 |
| Nov2 | SineWave4 | 3 | 5 | 1.23154 | 1.67860 | 1.47629 |
| Nov2 | helix | 0 | 1 | 0.68036 | 0.84065 | 0.80918 |
| Nov3 | circleSweep | 0 | 3 | 1.00465 | 1.36564 | 1.16264 |
| Nov3 | SineWave | 1 | 1 | 1.20376 | 0.82382 | 0.91425 |
| Nov3 | SineWave | 2 | 3 | 1.43673 | 1.35716 | 1.23764 |
| Nov3 | SineWave | 3 | 5 | 1.77377 | 2.07234 | 1.97720 |
| Nov3 | helix | 0 | 2 | 4.13058 | 2.99437 | 3.34507 |

**Interpretation:** GRU128 improves 4/10 fixed long examples relative to GRU64, while others worsen; these correlated examples are not independent statistical trials. Circle free RMSE worsens on Nov2 0.623->0.930mm and Nov3 1.005->1.366mm. It improves several sine/helix examples. The short-horizon winner is therefore a candidate for measured-state model augmentation, with mixed long-horizon accuracy. Retain the current GRU64 reference; do not infer improved physical tracking from short prediction scores.

## Bounded inverse self-consistency: fixed seed113

Each model generates its own target force at known logged pressure, then solves from a displaced initial pressure guess.
This checks numerical invertibility/state handling, not physical force accuracy, pressure correction quality or closed-loop stability.
Offline command bounds0..3bar; slew6bar/s at assumed.05s; fixed mean logged pressure. Candidate evaluations keep immutable pre-current state; selected state committed once.

| Architecture | Feasible cases | Max force-equation error (N) | Median (ms) | P95 (ms) | >50ms |
|---|---:|---:|---:|---:|---:|
| 1x64 | 208/208 | 3.57238e-05 | 2.181 | 5.410 | 2 |
| 1x128 | 208/208 | 9.16772e-05 | 2.349 | 5.198 | 0 |
| 2x64 | 208/208 | 5.35283e-05 | 2.647 | 6.223 | 0 |

Timing is this workstation only, excludes acquisition/control/communication, includes any cold/JIT outliers and is not a real-time deadline guarantee.
Force-rate measurements quantify predicted residual variation, not measured mechanical jerk or physical pressure smoothness.

## Retained artifacts and limitations

capacity_overview and six circle/sine/helix figures have editable FIG, vector PDF and300dpiPNG copies.
Predictions, CSV metrics, trial-window metadata, checkpoints, source snapshots, hashes and numerical verification retained.
See VERIFICATION.md. Existing models, raw data, normalizer, controllers, Simulink files and user edits preserved.
Three optimization seeds do not represent three independent hardware sessions. Architecture ranking is budget/data/protocol-specific.
A fresh unseen day/condition and matched PAC/PID FB controller tests remain necessary; no hardware, deployment, migration, commit or push performed.
