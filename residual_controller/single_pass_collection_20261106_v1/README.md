# Collection-policy single-pass block verification — 2026-11-06

User authorized the updated observation code after diagnosing the0.3bar
startup artifact. New source: rlResidualSinglePassCollectionShadowBlock.m.
Installation: SINGLE_PASS_COLLECTION_SHADOW.md. The old PAC-policy observation
source is preserved, not silently redefined. Five inputs, three outputs and
ResidualSinglePassDiagnosticsBus are unchanged. No SLX, controller, NDI,
reference, actual pressure command, frozen network, physical Kill or acquisition
wiring is changed. No hardware/model run, training, commit or push.

## Implemented change

Replace PAC/common-pressure allocation with the collection inverse's minimum
nonnegative pressure rule, followed by scaling at the configured maximum.
No common preload or trims are added. Same measuredq/v+supplied desiredddq
target, once-evaluated residual at actual pressure, immutable pre-current
history and subtraction sign. Pressure rails and actual-command-anchored slew
are retained. The unused common/trim/s-a settings and variableCommon input
are compatibility slots; this policy does not enforce the PAC correction caps.
The diagnostic bus exposes both pre-slew values and limited outputs.

The model remains pressure dependent, so this is a frozen-residual approximate
inverse, not a converged nonlinear solve, physical pressure validation or
complete tracking controller. The existing iterative inverse still uses PAC
allocation; it is not a matched solver comparator after this policy change.

## Numerical evidence

validateResidualSinglePassCollection passed:

- Both64/128 frozen widths and unchanged bus dimensions.
- Constructed feasible force target recovering[0,.02,.04]bar, verifying the
  subtraction sign and minimum-pressure convention.
- Identical outputs for either variableCommon input and changed unused PAC
  common/trim/correction-cap entries, with unchanged input contexts.
- Zero nominal force gives zero pre-slew pressure. Large requests scale to
  the configured maximum; rail and slew limits hold. From actual[3,3,3], a
  zero pre-slew nominal demand produces limited[2.7,2.7,2.7], confirming that
  the slew protection remains active. Missing/cold/nonfinite data returnNaN.
- All5266valid contexts from helix_FWD_INV_singlePass_test.mat replayed using
  the actual logged dd_length input, not a reconstructed target. Targets match
  the former block exactly; neural-force difference<=1.43e-14N. The new
  pre-slew allocation matches the separately calculated collection-policy
  counterfactual<=1.34e-15bar. Error decomposition difference<=1.01e-14N.

At firstvalid1.60s, actual commands arezero:

| Output,bar | P1 | P2 | P3 |
|---|---:|---:|---:|
| Nominal, before and after slew | .009561 | .038671 | 0 |
| Corrected, before and after slew | .022919 | .011821 | 0 |

All853zero-command contexts avoid the all-three-.3bar startup pattern, versus
853with the former PAC policy. Nominal/corrected pressures differ at all5266
valid samples before and after slew. Corrected requests scaled30times and
slew-limited861times (16.35%, formerly5099/5266=96.83%). Nominal slew limiting
remains3336times; this method does not make all demanded forces feasible.

Frozen residual equation mismatch median/P95/max=[4.72e-16,2.6530,7.0613]N.
Additional GRU evaluation at the NEW predicted pressure, offline only and with
no state commit, gives[.1010,3.2008,7.6572]N. Frozen and reevaluated mismatch
are different diagnostics; neither is measured physical force or tracking
error. The lower mismatch relative to the old policy is from an allocation
change, not retraining or evidence of improved closed-loop control.

## Generated/native component checks

The MATLAB Coder command emits C/header sources, but its standalone compiler
stage reports the existing Index exceeds array bounds error. Command result:
tmp/sp_collection_codegen_v1/codegen_result.json. This is not a successful
MathWorks standalone build or actual SLX build.

Emitted C was independently compiled with MinGW GCC. All764native-reference
rows pass, maximum absolute output/diagnostic difference1.99e-13. Cases include
both widths, both retained common-port values, invalid contexts and360real
circle/helix contexts. Evidence: single_pass_collection_native_20261106_v1
inputs/expected/native CSVs,compiler log,native_result.json. This qualifies the
component math only, not operator SLX compilation/full-loop timing or hardware.

Offline utility functions checkResidualSinglePassShadowCodegen,
writeResidualSinglePassCReplay and checkNativeResidualSinglePass.py now accept
an optional entry name for this same-ABI new variant; their default remains
the original entry. No runtime/controller code is changed by this extension.

## Retained artifacts / reproduce

summary.json,replay.csv/MAT,startupPNG/PDF/FIG,source/input SHA256,source and
generated-C snapshots. The source log/counterfactual is the independently
audited single_pass_helix_review_20261106_v2. Original raw files, packets and
old source are retained.

```matlab
addpath(fullfile(pwd,'residual_controller'));
validateResidualSinglePassCollection(pwd,fullfile(pwd,'tmp','new_collection_check'));
entry='rlResidualSinglePassCollectionShadowBlock';
checkResidualSinglePassShadowCodegen(pwd,fullfile(pwd,'tmp','new_collection_c'),entry);
writeResidualSinglePassCReplay(pwd,fullfile(pwd,'tmp','new_collection_fixtures'),entry);
```

```
python residual_controller/checkNativeResidualSinglePass.py tmp/new_collection_c/generated tmp/new_collection_fixtures rlResidualSinglePassCollectionShadowBlock
```

Use fresh output folders. The user installs/pastes/compiles the block manually;
it remains disconnected from valve commands for observation.
