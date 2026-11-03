# Equal actuator-length bounds: 2026-11-02 source correction

The limit was already **+/-0.035 m = +/-35 mm**. Do not increase it to fix the
directional discrepancy. The correction includes the previously omitted
actuator length l1=-l2-l3 in the existing soft restoring-force law.

Standing user reminder: apply this consistently when working on other PAC and
PDI/FB/PDFB controllers using this same reduced-coordinate model. Audit their
embedded and external forward/inverse/compensation code; changing shared source
does not migrate every saved model. See project_docs/DECISIONS.md D-12 and
project_docs/OPEN_ISSUES.md OI-10 for the project-wide handoff.

## Corrected force and scope

For q=[l2;l3], reconstruct l=[-l2-l3;l2;l3]. Apply the same scalar force F(l)
to each actuator, then map by virtual work:

    f_bound_q = [F(l2)-F(l1); F(l3)-F(l1)].

armS_single_boundForce.m is the shared implementation. Bounds, strength1e6
and smoothing mu2000 are unchanged. This is a soft restoring force, not a
hard guarantee that all lengths stay within35mm. Its returned coefficient
matrix is a force coefficient, not the derivative/tangent stiffness.

Current forward RHS, fixed/adaptive PAC force helpers and common-pressure
inverse source functions now use this shared force. No gains, adaptive update
laws, pressure limits, sensor processing or geometry/mass terms changed.
No SLX was edited and no hardware/Simulink model was run by the assistant.

## Collection model: operator installation

Do this with the model stopped, for a fresh run:

1. Keep the project root and RL_data_collection on the MATLAB path.
2. The existing forwardModel calls armS_single_entry -> armS_single_dynamics.
   On a fresh model rebuild it picks up the corrected external RHS/helper.
3. Replace the complete script in the collection **inverseModel MATLAB Function**
   with RL_data_collection/rlInverseModelBlock.m. It keeps the same three inputs
   and five outputs, invalid-sample hold and pressure-allocation policy.
4. Rebuild/start manually using the established operator procedure. A currently
   running or previously compiled model may still contain the old implementation.
   Record the code version and save a new recording; keep the earlier files.

Updating an external .m file does not replace an embedded inverse script.
The saved SLX still contains the old inverse until step3. Any embedded copies
of common-pressure inverse scripts also require their updated source to be
pasted if that particular model is used. PAC wrappers referencing external
force helpers pick up the correction on rebuild. This source correction is
not a claim that a hardware model has already been updated/qualified.

## Model versions and residual learning

- **actuator-bounds-v2 / boundMode1:** corrected default in armS_single_dynamics
  and armS_single_boundForce; used by current inverse/PAC force source.
- **coordinate-bounds-v1 / boundMode0:** historical l2/l3-only formulation.
  Set params.boundMode=0 to reproduce its RHS.

Old datasets, models, reports and source snapshots remain untouched. Existing
residualNominal and residualFrozenStep default to mode0 for the frozen two-day
screening contract; pass a fifth argument1 explicitly for a new corrected
baseline. greyboxSimulate similarly retains mode0 by default, with optional
fourth argument1. The historical calibration/nonlinear benchmark explicitly
select mode0. prepareResidualData remains a legacy two-cohort reproducer and
records cfg.boundMode=0; its old artifacts must not be silently relabeled.
The legacy screening reader rejects explicitly mode1 datasets.

A residual trained around the old nominal context must not simply be combined
with the corrected baseline. New corrected labels/context and any retraining
need a separately versioned dataset/model/evaluation. This is a concrete force-
mapping correction, not a renewed search for fitted nominal parameters.

## Offline verification

Results: residual_dynamics_learning/validations/actuator_bounds_20261102_v1.
Checks passed for independent three-actuator force mapping, potential gradient,
cyclic virtual-work symmetry, legacy-RHS agreement, forward/inverse/PAC static
force agreement, legacy/corrected nominal roundtrips, greybox causality, fixed/
adaptive PAC regressions and inverse dropout/fault handling.

Device-free same3bar ramp/hold counterfactual: legacy equilibria
3.0051/2.5268/2.5268rad; corrected equilibria2.5268rad for all three channels.
Transient responses need not be identical because other model terms remain.
This verifies removal of this directional artifact, not real-arm accuracy.

Numerical MATLAB checks passed. Standalone MATLAB Coder attempts hit a toolchain
initialization/build error after C source emission, including a shorter-path
retry. Full C build and actual Simulink integration remain unverified; see
verification.txt and codegen_status.csv. No pressure/hardware run was attempted.
