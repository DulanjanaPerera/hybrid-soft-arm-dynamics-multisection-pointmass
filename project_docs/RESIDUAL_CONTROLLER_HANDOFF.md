# Residual controller integration handoff

Prepared2026-11-03 for branch `codex/residual-controller-integration` after the
residual-learning checkpoint. No controller implementation or hardware run is
included in this handoff. Existing operator SLX files/rename and paper drafts
remain outside these source/results commits.

## Immediate sequence

1. Trace the active PDFB outer feedback and inverse path. Repository PID/PDI/PD
   names are not sufficient to infer an integral term. Start with a source-only,
   paste-ready PDFB augmented inverse, preserving the existing supervisor,
   pressure allocation, gains, filters and measurement-validity policy.
2. Derive consistency of the fixed-nominal learned residual with PAC's adaptive
   K/D/hysteresis/bias before adding it to PAC. Do not double-count existing
   adaptive compensation or claim the old passivity proof is preserved.
3. Perform offline controller/inverse/state/limit/runtime replay. Then prepare
   operator-run model-off/on comparisons with frozen weights and matched settings.
4. Evaluate actual length/velocity tracking, smoothness and pressure limiting.
   Use collected controller-domain data to design subsequent model improvement.
5. Defer the proposed HIL collect-fit-validate automation until this first
   implementation and matched comparison work. HIL note is a proposal only.

Do not run hardware models autonomously or edit SLX unless explicitly requested.
Provide paste-ready MATLAB Function code and log important findings/changes.

## Frozen model choices

Primary short-horizon candidate:
`residual_dynamics_learning/force_model/capacity_comparison/campaign_20261103_v1/candidate_force_model.mat`
SHA256 `8f2e0bb55ea9399b62abd5abe0b8f6dd15ab0470a400ad3e593e1edcd8f9f545`.
GRU1x128, fixed representative seed113 selected after a three-seed comparison.

Retained reference:
`residual_dynamics_learning/force_model/fine_tuning/campaign_20261103_v1/release/force_model.mat`
SHA256 `1c38edc36ba6ff4a41da8f974a88ca1c5ea775655426da91124a35b255a9c870`.
GRU1x64, old-data continuation, representative seed113.

Mean internal weighted scores across three seeds:64units .58531mm,
128units .48684mm,2x64 .53814mm. Representative Nov2 one-step vector RMSE
.36143->.28691mm; Nov3 holdout .75485->.64815mm. Long free circle errors worsen
with128: Nov2 .62291->.92961mm; Nov3 1.00465->1.36564mm. Retain both models and
these negative results; physical controller improvement remains untested.

## Model/runtime contract

Use `sgCell`, `sgZero` and `model.cfg.widths` from force_model/streaming. The
legacy `hiddenUnits` field/older cell are not the streaming model contract.
Load the exact numeric parameters and original eight-feature normalizer.
Features: q2,q3,v2,v3,p1,p2,p3,previous physical dt. Hidden memory uses a causal
observed-state prefix/reset; it evolves during inference with frozen weights.

The output is TWO generalized forces in N on the RHS:

    M(q)ddq+(C(q,dq)+D0)dq+K0q+bound(q)+G(q) = B*p+r(q,dq,p,h_previous).

Nominal: K0=650*[2 1;1 2], D0=40I; corrected all-three-actuator bounds,
l1=-l2-l3; pressure map B*p has zero activation. Align these definitions with
the active controller's pressure/activation map and adaptive decomposition.
It is not a network output of three pressures or an identified physical force.

Since the residual uses candidate CURRENT pressure, inversion is implicit.
Evaluate all candidate pressures with the SAME pre-current hidden state; commit
one state update per real sample using the chosen applied command. Respect
actual common-pressure, trim, pressure and slew settings, not copied offline
test limits. Refresh measured q/v on valid samples; specify invalid/recovery
state behavior explicitly. Offline interpolation previews are not live inputs.

Intended control period .05s.128-unit offline inverse self-consistency208/208,
max equation error9.17e-5N, median2.349ms/P955.198ms. These model-generated
targets and workstation times exclude the full hardware/control pipeline.

## Evidence and repository checkpoint

Read `RESIDUAL_CONTROLLER_COMPARISON_PLAN.md` for the full matched experiment
design, `PROJECT_STATE.md` and `RESIDUAL_DYNAMICS_RESULTS.md` for the latest
findings. Capacity/fine-tuning reports and verification live next to the models.
Current frozen numerical results and original models remain unchanged.

See `residual_dynamics_learning/VERSION_CONTROL.md` and its inventory for what
is versioned versus retained locally. Full raw recordings, generated training
datasets and intermediate optimizer checkpoints stay local pending the user's
deferred artifact migration. No remote push is requested or performed.
