# Residual controller implementation: offline verification

2026-11-03. Source-only PDFB/PAC wrappers implemented on
codex/residual-controller-integration. Read-only source audit of the user's
provisional FBL_commonP_circle_BRot.slx confirms an outer Ki*ie term and the
saved inverse's old bound. PAC saved adaptive-tunable interface matches source.
Both SLX inputs exactly match their HEAD bytes; no original controller helper,
feedback, gains, filters, trajectories or weights were edited.

## Implementation

- pdfbResidualInverseBlock: preserved first6 input/first8 output names/shapes,
  appended timing/enable/common/config inputs and status/diagnostic outputs.
  Virtual acceleration remains external, including existing integral feedback.
- pacResidualAdaptiveTunableBlock: existing23 inputs/28 outputs retained;
  append enableResidual/options and typed diagnostic bus. Model-off exactly
  delegates existing pressure supervision; full first28 outputs are identical.
- Fixed numeric64/128 GRUs, common128 padding, constant MAT packet loader,
  analytic dr/dp, bounded local inverse, activation/gravity re-expression and
  consistent augmented trim anchor. Reuse existing fixed/variable-common
  allocation and evaluate final post-slew candidates. No fmincon/runtime Python.
- PAC uses theta-theta0 remaining corrections; existing gamma/projection/rate/
  history/freeze logic preserved. MODEL_INTEGRATION.md states algebra/limits;
  no new passivity/stability theorem or physical force accuracy claim.
- Explicit warmup/gap/duplicate handling, one selected-pressure state commit,
  invalid hold/fault/slew and active model/network-change fault policy.
- setupResidualControllers prepares weights/path/typed Bus without loading a
  model. Installation inputs, diagnostics and configuration review in README.

## Canonical MATLAB checks: validation_20261103_v4

Frozen GRU cell mismatch0; analytic pressure Jacobian relative FD error
9.7972e-11. Activation-map identity for0/.2/.8 and PAC Y*(theta-theta0)
identity passed. All28 model-off PAC outputs identical in both common modes,
including adaptive updates/dropouts. Matched zero-activation zero-residual
pressure-path reduction passed.48 reachable numerical inverse cases:
max error8.8177e-6N; impossible requests remain bounded/limited. Cold-inclusive
inverse max12.4143ms/P951.4477ms on this host.

804 full PDFB/PAC synthetic-kinematic component ticks, both common modes,
passed pressure/slew/single-state-update/short-gap recovery/sustained-fault
checks. Wrapper P95.5787ms/max5.9104ms. Timing excludes acquisition, Simulink,
NI/communications and other system load. Synthetic imposed states are not
physical closed-loop tracking or force validation. Setup/bus definition test
also passed without opening a model. Earlier incomplete test attempts retained
with STATUS notes, not silently deleted or used as qualification.

## C emission, environment failure and native replay

MATLAB Coder emitted27/31 C files for PDFB/PAC, then returned the project's
existing generic toolchain error (index exceeds array bounds / compiler build
error). Generated reports/status retained in local
`tmp/residual_controller_codegen_20261103_v1`; the MATLAB build was NOT called
successful. No global MATLAB compiler configuration was changed.

The emitted source was compiled independently with local GCC13.1, explicit
MATLAB SDK includes and OpenMP. Both standalone DLL builds passed. Component
replay against MATLAB covered16 independent initialization cases /2576 ticks:
both controllers,64/128 units, residual-off/on, fixed/variable common, gravity
variation and short/sustained invalid samples. Max numeric discrepancy
5.00222085975e-12; pressure discrepancy 2.88657986403e-14bar.
11 explicitly enumerated logging-only allAboveActivation flag
comparisons at the exact.8bar threshold differ from roundoff (~3e-16bar).
Raw traces and exception counts retained; pressure/status/adaptation behavior
is numerically consistent. Native replay summary/source hashes:
`c_replay_20261103_v1/native_results.json` and generated_source_hashes.csv.
This is standalone component evidence, not a successful user SLX/Desktop
Real-Time build or deadline qualification.

## Comparison contract and next step

PDFB's NEW off/on pair shares the corrected three-actuator bound and PAC
supervisor/allocator. The legacy embedded PDFB path differs in activation,
blending/allocation, timestamp/enable/clock-reversal handling. Run a fresh
new-wrapper nominal baseline; do not credit old-versus-new allocator changes
to the learned model. Original first8 force outputs preserve shapes; available
force is now calibrated model-equivalent force, not a measured pressure force.

PAC first28 model-off outputs remain exactly compatible; model-on forceBias,
equivalent force and allocation residual use the augmented nominal definition.
Short gap memory hold with actual dt is an explicit deployment policy, distinct
from strict invalid-gap segmentation used in model fitting/evaluation.
Gravity re-expression preserves the trained model; it is not generalization
or adaptation to a new physical base orientation/payload.

Install sources manually in duplicates using README, confirm actual live
settings/timestamp origin/downstream pressure permit, compile with the operator's
model environment, then perform a model-off wiring/smoke check before enabling
the residual in a fresh trial. Physical gains, reference and filters stay
matched. HIL training and gain optimization remain deferred. No hardware,
Simulink edit/run, retraining, remote push or Git commit this implementation.
