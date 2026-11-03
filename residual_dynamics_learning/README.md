# Residual dynamics learning

## Current status: November3,2026

The generalized-force residual model is implemented, trained and evaluated.
Primary next workstream: source-only PDFB/PAC residual-model integration, then
operator-run matched model-off/on comparisons. HIL batch model updates follow
that first comparison; no automatic hardware or SLX execution is included.

- Primary short-horizon candidate: force_model/capacity_comparison/
  campaign_20261103_v1/candidate_force_model.mat (GRU1x128).
- Retained reference: force_model/fine_tuning/campaign_20261103_v1/
  release/force_model.mat (GRU1x64).128improves short prediction, but long
  circle free runs worsen. Do not infer better physical tracking yet.
- Canonical latest reports: capacity_comparison/campaign_20261103_v1/REPORT.md
  and fine_tuning/campaign_20261103_v1/REPORT.md under force_model.
- Runtime uses streaming/sgCell and model.cfg.widths; older cell schemas differ.
- Handoff: project_docs/RESIDUAL_CONTROLLER_HANDOFF.md.
- Versioned artifacts/local-only data policy: VERSION_CONTROL.md.

The sections below retain the earlier investigation history.

Fresh five-experiment evaluation completed (2026-11-02): canonical report
validations/fresh_20261102_v2/REPORT.md. Frozen weights/normalizer unchanged,
50,078 guarded anchors; supplied exports exactly match raw MAT re-exports.
Ridge improves 5 s path prediction over constant velocity (1.620 vs 3.408 mm),
but kinematic ablation remains competitive and no force residual is trained.
No controller or hardware qualification. This cohort is now evaluated.

Reproduce into a NEW output folder from the repository root in MATLAB:

    addpath(pwd,fullfile(pwd,'residual_dynamics_learning'));
    folder=fullfile(pwd,'residual_dynamics_learning','validations','CHOOSE_NEW_VERSION');
    validateFreshResidualModels(pwd,folder);
    assessFreshResidualValidation(pwd,folder);

These are offline data/predictor operations only. Models, raw recordings,
exports and normalizers are located explicitly in the source; neither helper
trains models or calls Simulink/serial/DAQ. V1 is incomplete and marked DO_NOT_USE.

Current nominal bound correction (2026-11-02): see
RL_data_collection/ACTUATOR_BOUND_CORRECTION.md and
validations/actuator_bounds_20261102_v1/REPORT.md. Shared current force model
includes all three actuator lengths at unchanged+-35mm. Historical residual/
greybox helpers retain explicit legacy mode0 by default to preserve existing
weights/reports; a corrected baseline uses mode1 and a new model/data version.
Numerical tests passed; standalone Coder toolchain/build check did not pass.
No SLX/hardware execution. Manual collection inverse script installation remains.

Historical November2 discussion: resume missing-dynamics residual learning with
a declared frozen nominal baseline; further nominal calibration is not a prerequisite.
The force-residual model was subsequently implemented and trained November3;
see the current status above. It differs from the earlier next-position screen. See the research plan and
latest specialized log. Retained paper-candidate figures are indexed in
project_docs/FIGURE_INDEX.md, with reproducible versioned exports under figures/.

Latest offline diagnosis (2026-11-02):
calibration/greybox_20261102_v2/REPORT.md. Constrained mechanical calibration
improves on the original nonlinear model; shared effective actuator memory adds
little and both remain worse than matched constant-velocity predictions.
Fitted damping reaches upper bounds and delay/T are confounded. No controller
changes or hardware qualification. Report contains reproduction commands,
parameter bounds, static/dynamic results, source snapshots and next diagnosis.

Step3 preliminary benchmark complete: see
benchmarks/screen_20261101_v2/REPORT.md. Fitted nominal-residual ridge/MLP/GRU
and diagnostic baselines are saved offline. Shortlist winner on validation is
ridge, but cross-day and shortcut ablations do not yet qualify a learned
missing-dynamics model for hardware. No deployment/control model changed.

Step2 complete: use datasets/two_day_v2/README.md for the prepared-data contract,
counts, quality limits and scripts. two_day_v1 is an invalid incomplete attempt
and must not be trained on. Trained screening models are described above.

Read MODEL_SPECIFICATION.md for completed step1: forward/inverse interfaces,
pressure allocation, target semantics, model shortlist and two-day split policy.
The six-step research plan is logged in
project_docs/RESIDUAL_DYNAMICS_LEARNING_PLAN.md; specialized findings are in
project_docs/RESIDUAL_DYNAMICS_RESULTS.md.

Reproduce interface audit from repository root in MATLAB:

    addpath('residual_dynamics_learning');
    auditResidualInterfaces(pwd,fullfile(pwd,'residual_dynamics_learning', ...
        'audits','new_interface_audit'));

Output directory must be new. Audit reads raw MATs, never merges/changes them,
and never opens/runs Simulink models. Primary results are under
audits/interfaces_20261101. SHA256 provenance and source extraction were added
in the accompanying offline inspection; raw data stays in its original paths.

To reproduce step2 in a NEW directory:

    addpath('residual_dynamics_learning','RL_data_collection');
    validateResidualNominal();
    prepareResidualData(pwd,fullfile(pwd,'residual_dynamics_learning','datasets','new_version'));
    % Then validateResidualData(newFolder), summarizeResidualData(newFolder),
    % verifyResidualRecordedMaps(newFolder).

Reproduce a new screening output with benchmarkResidualModels(dataset,newFolder),
then benchmarkResidualAblations and validateResidualScreening. Output folders
must be fresh. Defaults use CPU and fixed small training budgets. October27
has now been evaluated under frozen protocols; collect fresh confirmation for
later claims rather than treating it as perpetually untouched test data.
