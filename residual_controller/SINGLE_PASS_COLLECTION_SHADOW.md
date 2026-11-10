# Updated single-pass observation block: collection pressure convention

For the new desired-trajectory pressure use case, replace this historical
observed-state source with rlResidualDesiredPressureBlock.m and follow
[DESIRED_PRESSURE_BLOCK.md](DESIRED_PRESSURE_BLOCK.md). Its input2 takes the
full6x1 desired motion, not acceleration alone. This source is retained to
reproduce earlier recordings.

2026-11-09 clarification: this block reconstructs pressure around measured
q/v with supplied desired acceleration. It is NOT a desired-trajectory inverse
and has no desired q/v inputs. Its nominal pressure should not be expected to
equal the collection desired-state nominal command. Circle test4 confirms the
two different targets exactly; see
Experiments/Model_Val/target_mismatch_review_20261109_v1/README.md.
For the intended desired inverse comparison, a separate shared nominal force
demand is needed; keep actual-state/history conditioning of the GRU.

2026-11-06 operator verification: both new *_test_2 recordings (model1.53)
reproduce this source over all12218valid contexts with saved settings.
The0.3bar startup artifact is absent; scaled/slew behavior and validity pass.
See Experiments/Model_Val/collection_runs_review_20261106_v2/README.md.
This establishes operator execution of the observation block, not residual-on
controller performance or full-loop timing qualification.

Use **rlResidualSinglePassCollectionShadowBlock.m**. Copy the **entire file,
including its local minimumPressure function**, into the existing single-pass
MATLAB Function block. This is a separate source variant; the original
rlResidualSinglePassShadowBlock.m is preserved to reproduce earlier logs.
No SLX modification, model/hardware run or valve connection is performed by
these files. Outputs remain for Scopes/logging only.

## What changes

The network is evaluated once at actual pressure and measured q/v, with the
forward observer's immutable pre-current history. The force demand is still
the same measuredq/v plus supplied desired acceleration. Subtraction is still
correctedForce=nominalForce-residualForce. No desired length/velocity or outer
feedback input is added, and no inverse iteration or GRU state commit occurs.

Force-to-pressure conversion now matches the collection inverse convention:

```matlab
d = force / forcePerBar;
p1 = max([0; -d(1); -d(2)]);
p = max([p1; p1+d(1); p1+d(2)],0);
if max(p) > maximumPressure
    p = p*(maximumPressure/max(p));
end
p(p<=1e-8) = 0;
```

It adds no common-pressure preload. The nominal comparator uses the same
conversion without the residual. Pressure rails and the existing slew cap
are retained. The bus still exposes predictions before slew and limited
outputs separately. A genuine excessive change can still reach the0.3bar
step allowance; that number is not a learned minimum pressure or deadzone.

Only pressure conversion matches the collection code: this shadow target still
uses measuredq/v+desiredddq, whereas the collection inverse takes its supplied
desired q/dq/ddq. Consequently its nominal output need not match an existing
collection pressureRequested signal. Freezing r at actual pressure remains
an approximation because the network also depends on candidate pressure.

## Same ports, sizes and bus

Keep the block **Discrete0.05s** and its existing connections:

| Input | Size | Connect |
|---|---|---|
| inverseContext | Bus: ResidualShadowContextBus, scalar | Same forward observer output5, without delay |
| desiredAcceleration_mPerS2 | double2x1 | Same desired[ddl2;ddl3] line |
| shadowPressureOptions | double10x1 | Same existing Constant |
| variableCommon | Boolean scalar | Keep existing connection; ignored by the new minimum-pressure policy |
| singlePassOptions | double2x1 | Same[network width;frozen balance tolerance] Constant |

The three outputs and their names/shapes are unchanged:

1. singlePassPressurePrediction_bar,double3x1, **after slew**;
2. nominalPressurePrediction_bar,double3x1, **after slew**;
3. singlePassDiagnostic,scalar **Bus: ResidualSinglePassDiagnosticsBus**.

There are no new ports or bus fields. If the current base workspace already
contains the bus/options, no setup change is needed. To recreate them while
the model is stopped, from the project root with the existing shadow setup:

```matlab
addpath(fullfile(pwd,'residual_controller'));
[singlePassOptions,ResidualSinglePassDiagnosticsBus] = ...
    setupResidualSinglePassShadow(pwd);
singlePassOptions(1)=shadowOptions(1);
```

The new code uses these existing pressure-option entries:

| Entry | Requirement/default | Meaning |
|---|---|---|
| 4 | 0 | Minimum pressure, matching collection convention |
| 5 | 3bar | Maximum pressure; excessive requests are scaled |
| 8 | 6bar/s | Retained pressure-slew allowance |
| 9 | 0 | As-trained zero-activation nominal pressure map |
| 10 | pi*.0065^2*1e5 N/bar | Must match frozen packet |

Entries1:3(common/trim) and6:7(PAC correction caps) are **not used** in this
collection convention; there is no need to change their existing Constant
values. Port4variableCommon likewise no longer selects a common-pressure
policy. Keep it connected for port compatibility. The0..maximum rail and
slew bounds still apply; do not describe the PAC s/a caps as retained here.

## What to view

Use the current Bus Selector to display/log:

- correctedPressureBeforeSlew_bar and nominalPressureBeforeSlew_bar;
- residualForce_N and correctedForceRequest_N;
- allocationLimited (now peak-pressure scaling),slewLimited and valid;
- frozenBalanceError_N and frozenBalanceSatisfied.

Compare output1 andoutput2 separately for the limited outputs. Each slew box
is still anchored to the actual pressure in the context at that observation,
not the prior hypothetical prediction. Missing/cold/invalid contexts produce
NaN predictions, valid=false; do not replace them with held/zero measurements.
The forward observer retains ownership of history and readiness.

The old iterative block still uses the PAC policy. Its pressures therefore
are **not a matched solver comparison** with this new variant. Use the new
block's own nominal output as its comparator. Matching an iterative variant's
allocation policy would be separate work; neither PAC nor PDFB is changed.

## Offline evidence and manual next run

Canonical checks: single_pass_collection_20261106_v1/README.md. Both frozen
widths, sign/force construction, zero-force/no-preload, scaling, slew, invalid
handling, unchanged context and compatibility-port invariance pass. All5266
valid saved helix contexts match the independently calculated collection
mapping. At firstvalid1.60s:

- Nominal:[.009561,.038671,0]bar.
- Corrected:[.022919,.011821,0]bar.

There are0all-three-.3bar outputs among853zero-command contexts, versus853
with the former PAC policy. This is offline replay of the saved recording,
not a new hardware test or proof of physical pressure/tracking accuracy.

After your manually started observation test, save to a **new filename**:

```matlab
save('helix_FWD_INV_singlePass_collection_r01.mat','out','shadowOptions', ...
     'shadowPressureOptions','singlePassOptions','-v7.3');
```

Keep the actual commands, reference/filter/Kill/acquisition wiring intact.
The operator still needs to compile the pasted block in the actual SLX;
numeric and native component verification do not establish full-loop timing.
