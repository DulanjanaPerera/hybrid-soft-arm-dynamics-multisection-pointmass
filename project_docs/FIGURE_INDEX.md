## 2026-11-03 - Three-seed GRU capacity comparison completed

Seven capacity comparison figure bundles at residual_dynamics_learning/force_model/capacity_comparison/campaign_20261103_v1:
capacity_overview plus motion_comparison/nov[2,3]_[circle,sine,helix]_before_after. Each FIG/PDF/300dpiPNG.
Architecture legends1x64/1x128/2x64; fixed seed113. Long predictions are known-command free runs, not controller tracking.
All Nov3 windows held-out, missing truth unfilled. MAT/CSV and source snapshots retained; all negative results preserved.

## Retention convention

Keep figures in versioned analysis folders; never overwrite earlier runs.
For new figures retain an editable MATLAB .fig, vector .pdf, 300 dpi .png,
plotting source, underlying numerical data and provenance. Keep raw experiment
files at their existing locations and reference acquisition/run/split IDs.
When only a legacy raster exists, preserve it and identify that limitation;
do not present a raster wrapped in PDF as a vector original.

Each bundle should record captions, units, scoring definition, selections,
parameters, filters/timebase, fitting/evaluation split, and limits on claims.
Retain per-run/per-trial data rather than just rounded aggregate bars.
Screenshots are useful context but are not substitutes for quantitative data.
Local retention is separate from Git commits: commit/push only when requested.

## 2026-11-02 model diagnostics v1

Folder: residual_dynamics_learning/figures/20261102_model_diagnostics_v1.
Plot source: residual_dynamics_learning/archiveResidualFigures.m.
Bundle includes a source snapshot, source reports, copied CSVs, plotted-value
CSVs, source-path manifest and SHA256 manifest. Reproduction reads saved results
only; it does not refit models or run simulations/hardware. Details and captions
are in the bundle README.

| ID | Content | Interpretation / paper-use limit |
|---|---|---|
| F01 | Nonlinear prediction error versus physical horizon, with matched baselines | 30 selected October19 validation windows; prediction, not controller tracking; overlapping windows are not independent repeats |
| F02 | Static equilibrium error for all four validation holds | Illustrates remaining equilibrium mismatch; small selected set, no identified physical stiffness claim |
| F03 | Training-only damping profile with other parameters fixed | Diagnostic of parameter compensation; not a damping measurement or a reason to expand controller gains |
| F04 | Preliminary forward-screening validation/external errors | One native sensor step, log y-axis; cross-day timing/coverage differs; kinematic ablation is included |
| F05 | Original example nonlinear length rollouts, one per waveform | Preserved original PNG and source reference; examples are not all windows and are not evidence of control improvement |

F01-F04 have .fig/.pdf/.png. F05 is an existing raster only; its original
source and underlying windows remain in the calibration output. Do not compare
F04 one-step values directly with F01 multi-step values: models/windows/horizons
differ. None establishes residual-model success or readiness for arm testing.

## Other retained sources

- 2026-11-03 force-GRU fine-tuning comparison:
  residual_dynamics_learning/force_model/fine_tuning/campaign_20261103_v1.
  fine_tune_overview plus motion_comparison/nov2_circle/sine/helix_before_after
  and nov3_circle/sine/helix_before_after. Seven FIG/PDF/300dpiPNG bundles;
  full trajectories/forces/truth/masks/pressure/time, score/windows CSVs,
  checkpoints/optimizer states, source/model hashes retained. Original, old-data
  continuation and old+Nov3 replay, including regressions. Nov3 motion windows
  all heldout; helix uses trial2, unlike prior TRAIN-trial1 illustration. Frozen
  weights, measured initialization then no state resets; no physical tracking claim.

- 2026-11-03 November3 frozen force-model motion plots:
  residual_dynamics_learning/figures/motion_predictions_20261103_v1/frozen_model_predictions_v1.
  circle/sine/helix_free_prediction and matching *_one_prediction, plus nominal
  separately. Seven FIG/PDF/300dpiPNG sets, full predictions/pressures/physical
  times/masks MAT, score/window/coverage CSVs and source hashes. Chronological
  first examples, same weights and plotting conventions as November2; faster
  cycles and larger helix. No training, actual-gap filling or new controller run.

- 2026-11-03 requested November2 frozen force-model motion plots:
  residual_dynamics_learning/figures/motion_predictions_20261102_v1.
  circle/sine/helix_free_prediction and circle/sine/helix_one_prediction;
  nominal_free_prediction separately. Seven FIG/PDF/300dpiPNG bundles with
  predictions MAT, score/selection CSVs and source hashes. Exact frozen GRU64
  one-step versus long free propagation; inferred-length/angle/PCC response,
  not new controller tracking or independent sensor XYZ. Missing labels blank,
  free predictions carry recorded pressure/time with no measured-state resets.

- 2026-11-03 bounded NDI interpolation feasibility:
  residual_dynamics_learning/figures/motion_predictions_20261103_v1/short_gap_audit.
  short_gap_audit FIG/PDF/300dpiPNG: actual gap counts/spans, artificial-mask
  errors and explicit one-frame candidate. Full raw/preview flags/values and
  source/error CSVs retained. Actual gap truth unknown; no training, validation
  or runtime filled measurements. Read REPORT.md before using for a paper.

- 2026-11-03 streaming force-GRU sequence:
  residual_dynamics_learning/force_model/streaming/evaluation_20261103_v1.
  streaming_comparison, circle_streaming_example, circle_short_free_example
  each FIG/vectorPDF/300dpiPNG. Six-model old-only screen/learning curves,
  repeated-seed results, frozen Nov2 comparison versus nominal/CV; explicitly
  separate measured-state one-step from free integrated prediction. All source,
  predictions, force spikes, pressure/time, CSV/MAT and SHA256 retained.
  Prediction/model-inverse self-consistency only, no controller/physical-force claim.

- 2026-11-03 first generalized-force GRU:
  residual_dynamics_learning/force_model/evaluation_20261103_v1.
  force_model_evaluation FIG/vectorPDF/300dpiPNG: older-only selection curve,
  November2 matched one-step/2s/5s inferred-length prediction comparisons.
  Source, all per-window trajectories/forces/pressure/time, CSV/MAT/hashes and
  model release retained. New outputs are forces integrated through corrected
  mechanics, not old position corrections. Show weak short-horizon results as
  well as favorable5s scores. No force-truth, controller or stability claim.

- 2026-11-03 requested sine prediction grouping:
  residual_dynamics_learning/figures/sine_predictors_20261103_v1. Four sets:
  sine_free_running_cv_gru, sine_free_running_others, sine_one_step_cv_gru,
  sine_one_step_others (FIG/vectorPDF/300dpiPNG). Rows P1/P2/P3, columns l2,l3,
  phi; measured black. Free predictions use one continuous mature ~31.4s
  window/channel, initialized once; full-record one-step uses measured history
  every sample. Different start phases/trials, no hidden gap bridging. Spikes
  retained. Grouping does not imply GRU wins P2 (ridges better there). Model
  predictions, not controller tracking or force-residual results. CSV/MAT,
  window provenance, source/hash and verification/report retained.

- 2026-11-03 requested phi1 circle comparison of five existing predictors:
  residual_dynamics_learning/figures/circle_phi1_predictors_20261103_v2.
  full_circle_xy and full_circle_angles: first chronological20s circle,
  observed initial history then no measured-state resets. Requested red circle,
  measured black response, colored predictions; XY is reconstructed PCC, not
  independent trackerXYZ. one_step_five_circles_angles_labeled: measured history
  supplied eachsample across100s; keep separate from autonomous propagation.
  FIG/vectorPDF/300dpiPNG, CSV/MAT/source/hashes/report retained. Requestedphi1,
  actualphi.052-.230; GRU closest length/XY on this single circle, no control or
  architecture-winner claim. Original legend export and incompletev1 preserved.

- 2026-11-02 frozen predictors on five fresh experiments:
  residual_dynamics_learning/validations/fresh_20261102_v2. Two FIG/vectorPDF/
  300dpiPNG sets: fresh_prediction_summary (per-run one-step error and .5/2/5 s
  path error) and fresh_rollout_examples (first eligible 5 s window per run,
  including unfavorable near-straight drift, no error-based example selection).
  Forward/inverse trial/phase scores, raw-export equality, hashes, source,
  selected-window/row manifests and q2/q3 traces retained. Errors are inferred
  length prediction, not controller tracking. Locally frozen nominal transition,
  known pressures/query times, no measured-state resets. Physics-prior benefit
  remains mixed; no full-circle, force-residual or controller-success claim.

- 2026-11-02 SineWave4 same-reply status review:
  residual_dynamics_learning/audits/validity_sine4_20261102_v1. Two FIG/vectorPDF/
  300dpiPNG sets, event/native/packet-size/summary CSVs and source hashes. Main
  plot: channel3 pressure/measuredphi, raw tipstatus2 versus orientation/base
  validity, last-valid XZ poses before losses. Missing poses unknown; no
  interpolation.10recorded MISSING events, no unique physical cause established.
  Acquisition-quality evidence, not tracking or learned-model improvement.

- 2026-11-02 SineWave3 recurring losses:
  residual_dynamics_learning/audits/validity_sine3_20261102_v1. SineWave3_validity
  and dropout_pose_comparison as FIG/vectorPDF/300dpiPNG, full native/event/neighbor
  data and source hashes. Pose plot is XZ projection: grey every tenth valid
  sample, all markers actual last-valid poses before missing frames, colored
  by phi. Missing poses are unknown; no interpolation. Individual runs are not
  registered to one world frame. Several loss regions and intermittent valid
  returns do not identify a hardware/field cause. Initial3D attempt retained
  separately after vector-export failure. Quality diagnostics only.

- 2026-11-02 corrected-bound counterfactual:
  residual_dynamics_learning/validations/actuator_bounds_20261102_v1.
  nominal_bound_counterfactual FIG/vectorPDF/300dpiPNG, full time-series/summary
  CSV, source hashes and verification. Caption: six device-free nominal ODE
  simulations, individual0->3bar ramp/hold, unchanged mechanics except legacy
  versus corrected actuator-bound force. Equilibrium asymmetry disappears.
  No experimental measurements/fitting or real-arm accuracy/stability claim.

- 2026-11-02 nominal bound asymmetry:
  residual_dynamics_learning/audits/bound_asymmetry_20261102_v1. Two per-run
  FIG/vectorPDF/300dpiPNG sets with numerical plot/peak tables and source hashes.
  Rows show first sine trial per channel; left predicted versus orientation-
  measured phi, right all three simulated lengths with+-35mm reference lines.
  The implemented soft bound acts on l2/l3 only, illustrating a directional
  modeling artifact, not identified hardware stiffness. Soft restoring force
  starts before the35mm line. No corrected-model counterfactual shown.

- 2026-11-02 fresh sine validity audit:
  residual_dynamics_learning/audits/validity_20261102_v1. Two FIG/PDF/300dpiPNG
  plots (SineWave_validity, SineWave2_validity), all event/native-interval CSVs,
  full exported input tables, plotting source, hashes and report. Captions: upper
  panel configuration flag with all invalid samples marked; middle P3 command;
  lower first channel3 tip/base/configuration flag comparison. Isolated tip
  dropouts explain every configuration rejection; detailed physical cause not
  identified. Quality diagnostics only, no predictive/control-performance claim.

- Calibration: residual_dynamics_learning/calibration/greybox_20261102_v2,
  including original example_validation_rollouts.png, reports, selected windows,
  fit stages, sources and per-window numerical results.
- Earlier residual screening: residual_dynamics_learning/benchmarks/screen_20261101_v2,
  including trained weights, protocols, per-run scores, rollouts and ablations.
- Existing PAC presentation plots and experiment-analysis outputs remain in their
  original experiment folders. This index does not claim that those older
  artifacts have all been inventoried or audited for publication.
