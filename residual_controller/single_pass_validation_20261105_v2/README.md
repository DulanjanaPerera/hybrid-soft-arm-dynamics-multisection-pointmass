# Single-pass shadow inverse verification — 2026-11-05

User authorized a separate block implementing one residual evaluation followed
by subtraction and pressure allocation. Source and manual wiring are in
`residual_controller/rlResidualSinglePassShadowBlock.m` and
`residual_controller/SINGLE_PASS_SHADOW.md`. The existing iterative inverse,
shared controller allocator, forward observer, NDI source, packets and SLX
files were not changed. No model/device/pressure execution, training, commit
or push occurred. Documentation and offline helpers/artifacts were added.

## Contract

The nominal force demand remains the existing shadow demand: measured q/v
and supplied desired acceleration. The block computes exactly one GRU force
at the actual current pressure and immutable pre-current history, subtracts
that force, and calls the existing bounded nominal pressure allocator. It
then clips pressure to the same actual-command-reanchored slew box. There is
no pressure iteration, dynamics integration, neural state commit or persistent
state in this block. Forward readiness/missing-sample policy is unchanged.

This is a **frozen-residual approximate inverse**: the GRU depends on pressure,
so r(pActual) need not equal r(pPrediction). The bus reports frozen equation
balance, not nonlinear convergence or measured force accuracy. No actuator
connection or full tracking controller is provided. Pre-slew diagnostics are
hypothetical allocations, still subject to rails and correction caps.

## Numeric checks

MATLAB R2025a `validateResidualSinglePassShadow` passes:

- Bus element count, both64/128 network widths, a constructed exactly feasible
  force target recovering known pressures (including subtraction sign).
- Identical repeated calls, unchanged input contexts and no hypothetical
  accumulation; large-demand rails/slew and invalid/cold/nonfinite handling.
- Every nonlimited pre-slew allocation has force mismatch below1e-8N.
- Full circle/helix saved-context replay, including unavailable rows. Validity
  matches the existing inverse; 6908+5056=11964 valid rows. Nominal target and
  nominal post-slew pressures match logged results exactly, GRU force within
  1.07e-14N, error decomposition within1.31e-14N.

The saved recordings did not retain complete original option vectors.
Source-default variable-common settings previously reproduced the existing
inverse; the new replay uses those same settings and validates the nominal
comparator/force target against the logs. This does not uniquely identify
every hidden historical option or establish physical inverse accuracy.

## Saved-context results

| Quantity | Circle | Helix |
|---|---:|---:|
| Valid rows | 6908 | 5056 |
| Equal corrected/nominal before slew | 0 | 0 |
| Equal corrected/nominal after slew | 5440 (78.75%) | 3491 (69.05%) |
| Allocation limited | 236 | 79 |
| Corrected slew limited | 6646 | 4891 |
| Maximum pressure difference before slew | 1.1135 bar | 0.9767 bar |
| Maximum pressure difference after slew | 0.5799 bar | 0.6000 bar |
| Frozen balance <=0.02N, after slew | 110 | 201 |
| Candidate-pressure model balance <=0.02N | 23 | 1 |

The correction is visible before slew at every valid sample. Retaining the
slew constraint still masks much of it after clipping. These results do not
justify removal of the real pressure protection.

The offline analyzer makes an **additional GRU call outside the real-time
block** at pPrediction to quantify the pressure-dependence approximation:

| Norm of model equation mismatch | Circle median | Helix median |
|---|---:|---:|
| New block, residual frozen at actual pressure | 0.4197N | 0.4427N |
| New pressures, residual reevaluated at predicted pressure | 0.4083N | 0.5815N |
| Existing iterative inverse from saved log | 0.3772N | 0.4841N |
| Nominal inverse from saved log | 9.4967N | 7.7555N |

Candidate-pressure mismatch P95 is4.3583/4.8579N. Neither the new single-pass
approximation nor the existing constrained local solver is an exact inverse
on these recordings. The single-pass method is simpler and numerically
comparable in median mismatch here, not demonstrated superior. All values
are model algebra errors, not measured physical force or tracking error.
Pressure commands are not independently measured chamber pressures.

Measured MATLAB numeric P95 time is0.144/0.132ms for this block on these
saved contexts. This excludes acquisition, scheduler, forward observer,
plotting/logging and the offline extra evaluation; it is not full-loop
real-time qualification or a worst-case bound.

## Generated C / independent native verification

The fixed-size MATLAB Coder command emitted C/header sources but reported
the existing standalone compiler/toolchain `Index exceeds array bounds`
failure. The raw command result is retained at
`tmp/single_pass_codegen_20261105_v1/codegen_result.json`. This is **not** a
successful MathWorks standalone build or an actual Simulink-model build.

The emitted source was separately compiled with MinGW GCC. All764 reference
rows pass independent native replay, maximum absolute difference
1.92e-13. Fixtures cover both widths, fixed/variable common allocation,
invalid contexts and360 distributed real circle/helix rows. Every output
and diagnostic field, including NaN masks/Boolean flags, was compared.
Evidence: `residual_controller/single_pass_native_20261105_v1`, with
inputs/expected/native CSVs, native_result.json and compiler log. Source and
generated-C snapshots/hashes are retained with this validation folder.

The operator still needs to paste/type/wire and compile the new block in
their SLX. No pressure output should actuate hardware in this observation test.

## Retained artifacts / reproduction

- `validation.json/.mat`, `summary.csv` and per-case `samples.csv/replay.mat`.
- Per-case `single_pass_pressure_comparison.png/.pdf/.fig`: pre/post-slew
  pressure overlays with actual missing/cold intervals preserved as gaps.
- Source snapshots, frozen packet and raw input SHA256 references.
- v1 numerical replay was an earlier pass before moving the evaluation-count
  diagnostic assignment immediately after the GRU call. Valid-row numeric
  outputs are unchanged; v2 is canonical and adds explicit nonlimited-force
  assertions. v1 is preserved, not overwritten or presented as final provenance.

```matlab
validateResidualSinglePassShadow(pwd,fullfile(pwd,'tmp','single_pass_new_numeric'));
checkResidualSinglePassShadowCodegen(pwd,fullfile(pwd,'tmp','single_pass_new_codegen'));
writeResidualSinglePassCReplay(pwd,fullfile(pwd,'tmp','single_pass_new_fixtures'));
```

Then independently compile/replay emitted C (offline):

```
python residual_controller/checkNativeResidualSinglePass.py tmp/single_pass_new_codegen/generated tmp/single_pass_new_fixtures
```

Use fresh output folders. No raw data, production settings, neural weights,
sensor-validity rules or protections need changing to reproduce this review.
