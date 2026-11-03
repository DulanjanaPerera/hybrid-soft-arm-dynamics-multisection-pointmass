# Three-actuator bound correction: offline verification

2026-11-02. User authorized correcting the directional bound issue. The common
limit remains[-0.035,+0.035]m, strength1e6,mu2000. This is not a bound increase.
Shared helper maps three scalar restoring forces by J'F, J=[-1 -1;1 0;0 1].
See RL_data_collection/ACTUATOR_BOUND_CORRECTION.md for exact installation,
scope, version defaults and the paste-ready collection inverse block.

## Numerical checks passed

- Independent scalar actuator forces versus generalized-force mapping.
- Cyclic actuator permutation/virtual-work identity: maximum relative error
  6.5344e-17 over zero/interior/near/outside-bound poses.
- Force is the gradient of independently integrated scalar potential:
  maximum relative finite-difference error5.1067e-10.
- Positive semidefinite force coefficient and nonnegative q'*f checks.
  Coefficient is not the tangent stiffness.
- Isolated l1-bound case q=[-18;-18]mm now produces the missing restoring
  contribution; legacy mode remains approximately zero there.
- Legacy RHS versus independent former diagonal-K formula: maximum relative
  error3.2839e-16.
- New fixed/adaptive PAC static force and paste-ready inverse agree with
  the corrected forward-model force at all tested poses.
- residualNominal both modes matches RHS, including implied-l1-bound poses;
  feasible inverse roundtrip/allocation checks pass.
- Greybox original-RHS/future-state independence/delayed-input tests pass in
  its declared legacy mode.
- Fixed PAC inertia/skew/derivative/force/energy/invalid/base/gravity tests and
  matched-model ODE pass; maximum balance error1.4211e-14N, final length error
  2.9021e-11m. This is a synthetic matched model, not measured tracking accuracy.
- Adaptive PAC force/regressor/history/projection/ideal-energy/freshness/fault
  tests pass; common-pressure inverse single-frame/sustained-fault tests pass.
- Against archived pre-correction source helpers, legacy frozen transition
  error0 and nominal helper maximum relative difference1.44e-16.

## Device-free counterfactual

Six offline ode15s runs start from rest; one channel at a time ramps0->3bar in
one second then holds through12s. Mechanical parameters/geometry/mass/gravity,
pressure map and damping are identical; only boundMode differs.

| Mode | Channel | Final phi (rad) | Maximum phi (rad) |
|---|---:|---:|---:|
| Legacy0 | 1 | 3.005137 | 3.038081 |
| Legacy0 | 2 | 2.526837 | 2.526852 |
| Legacy0 | 3 | 2.526837 | 2.526852 |
| Corrected1 | 1 | 2.526837 | 2.528654 |
| Corrected1 | 2 | 2.526837 | 2.526852 |
| Corrected1 | 3 | 2.526837 | 2.526852 |

Corrected equilibrium spread4.885e-15rad. These simulations test the identified
nominal artifact; no actual-arm measurements are used as targets and no gains/
parameters are fitted. Other model terms can still produce transient direction
dependence. Remaining hardware mismatch is not resolved by this result.

counterfactual_summary.csv and counterfactual_trajectories.csv retain every
value. nominal_bound_counterfactual is retained as FIG/vectorPDF/300dpiPNG.
Caption: same nominal arm under individual3bar ramp/hold inputs; original
coordinate-only soft bounds versus equal three-actuator soft bounds. No
experimental-control, learned-model or hardware-stability claim.

## Build status and limitations

Numerical MATLAB verification completed. Additional standalone Coder checks
failed during the selected Visual C++ toolchain initialization/build step
with 'Index exceeds array bounds', after emitting bound-helper C source.
A shorter workspace path did not resolve the failure. Further per-entry
outcomes are recorded in codegen_status.csv; emitted source is not a passed
library build or Simulink compilation. Failed attempts are retained rather
than overwritten. No installed compiler settings were changed.

No SLX/hardware model was compiled/run. The collection inline inverse still
requires the operator to paste rlInverseModelBlock.m before a paired updated
model test. Numerical/standalone checks cannot qualify hardware behavior.
Raw recordings, old datasets/weights/results and historical source snapshots
remain intact. Current baseline default is mode1; historical reproducer paths
explicitly retain mode0. source_snapshot/source_manifest.csv retain exact
corrected implementations/dependencies and hashes.

Reproduce numerical checks from repository root:

    addpath('residual_dynamics_learning');
    validateArmSBounds(newOutputFolder);

The destination must be new. Additional regression commands and all numerical/
Coder output are retained in verification.txt. No commits or push performed.
