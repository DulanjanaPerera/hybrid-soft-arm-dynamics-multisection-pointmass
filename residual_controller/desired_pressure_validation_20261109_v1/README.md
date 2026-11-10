# Desired-motion pressure block: offline verification

2026-11-09. New paste-ready rlResidualDesiredPressureBlock.m uses6x1 desired
q/v/a for its nominal force target; both branches share that target. A single
GRU call stays conditioned on actual measured q/v, actual command, physical dt
and immutable actual history. Historical observation blocks remain unchanged.
Manual wiring/semantics: ../DESIRED_PRESSURE_BLOCK.md.

## Numeric evidence

validateResidualDesiredPressure.m passed for both64/128-unit frozen packets:
desired force/sign construction, independent original collection inverse,
desiredq/v sensitivity with unchanged observed residual, repeatability/context
immutability, ignored compatibility inputs, zero nominal demand/no preload,
maximum-pressure scaling, actual-command slew bounds and invalid handling.
Cold context, bad state/pressure/time and all6 NaN reference positions return
NaNs/validfalse without evaluating the GRU. Bus schema unchanged.

Every saved test4 row was replayed:5433rows,5118valid and315invalid.
No dropout interpolation, gate bypass or recurrent-state commit occurs.
The source MAT saves only out; the explicit pressure settings in summary.json
were independently verified against its existing single-pass behavior in the
preceding target_mismatch_review_20261109_v1 audit.

| Check | Maximum difference |
|---|---:|
| Desired force vs recorded tauDesired |3.56e-15N|
| Actual-state residual vs recorded residual |1.07e-14N|
| Independently mapped pre-slew pressures |1.34e-15bar|
| Independently clipped post-slew pressures |1.12e-15bar|
| Nominal pre/post vs collection pressureRequested |4.45e-16bar|
| Allocation + slew error identity |1.01e-14N|

First70s corrected peaks before AND after slew:[.5228,.4682,.4835]bar,
versus recorded nominal[.2423,.2403,.2423]bar. No corrected slew limiting
in that interval. The included PNG/PDF/FIG shows the expected larger demand
and the new nominal comparator coincident with the original nominal command.

## Limits of this replay

Across the entire recording, corrected maximum-pressure scaling occurs126
times and slew clipping3027/5118valid rows; nominal slew clipping0. Shadow
slew remains anchored to the recorded actual command at each observation.
It cannot progressively ramp an imaginary command while actual pressures
remain unchanged. Consequently many larger desired demands are infeasible
within that per-observation box in this counterfactual replay.

Frozen-residual force-balance error norm median/P95/max:[7.107,14.143,19.958]N
after limits. Reevaluating the pressure-dependent model at the new limited
candidate gives[7.615,14.541,20.370]N. These are model diagnostics, not
physical measured force errors. They are retained rather than reporting
pressure as though it achieved the requested motion. First70s peaks are
illustrative inverse approximations, not validated hardware requirements.

## Generated C / independent native verification

checkResidualSinglePassShadowCodegen now accepts the new entry with a fixed
6x1 reference while preserving existing2x1 entries. MATLAB Coder emitted
rlResidualDesiredPressureBlock.c and fixed-size6-input header, but its local
standalone build failed with the existing 'Index exceeds array bounds' /
'C compiler produced errors' toolchain failure. It is NOT a successful
MATLAB standalone build or user-SLX compilation claim.

Independent MinGW GCC compilation and native replay passed764 fixtures,
including both widths, both compatibility switch values, synthetic states,
recorded histories/reference motions, and invalid/cold contexts. Maximum
absolute difference across outputs/diagnostics:2.203e-13. Fixtures, expected/
native CSV, build log and native_result.json are in
../desired_pressure_native_20261109_v1. The generic native checker now accepts
the6-element reference for this entry, preserving2-element legacy layouts.

No Simulink model, hardware, acquisition, pressure or training was run.
No SLX changed. No commit/push performed. Operator model compilation and
real-time timing/physical validation remain to be performed manually.

## Artifacts and reproduction

- summary.json,replay.csv,replay.mat.
- desired_pressure_first70.png/.pdf/.fig.
- source_snapshot and source/packet/MAT/SLX hashes.
- codegen_result.json and generated C/header snapshot.
- See artifacts_sha256.json for retained file hashes.

From the project root, offline only, using fresh output directories:

```matlab
addpath('residual_controller');
validateResidualDesiredPressure(pwd,freshValidationFolder);
checkResidualSinglePassShadowCodegen(pwd,freshCodegenFolder,'rlResidualDesiredPressureBlock');
writeResidualDesiredPressureFixtures(pwd,freshNativeFixtureFolder);
```

Then run checkNativeResidualSinglePass.py with generated-folder,fixtures-folder,
and rlResidualDesiredPressureBlock arguments. No model build/run is involved.
