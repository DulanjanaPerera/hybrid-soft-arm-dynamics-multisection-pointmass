# Updated single-pass operator runs — 2026-11-06

The two new recordings confirm the **collection-policy single-pass block is
installed and operating as implemented**. The0.3bar startup artifact is gone,
the residual pressure correction is visible, and all checked pressure rails,
source-command and slew identities pass. Remaining limitations include
pressure limiting, the frozen-residual approximation and tracker dropouts.

Offline audit only: no hardware, acquisition, pressure or Simulink model was
started; no SLX/controller/NN/settings/raw-data changes, training, commit or
push. Raw filenames found:

- circlesweep_FWD_INV_test_2.mat (the supplied .at extension was a typo).
- helix_FWD_INV_singlePass_test_2.mat (the supplied shorter helix filename
  does not exist; this is the new single-pass recording).

Both are model residual_FWD_INV_validation1.53, ReachedStopTime with
no saved model errors. Unlike previous logs, both contain the actual saved
shadowOptions,shadowPressureOptions,singlePassOptions,shadowInverseOptions.
The new block uses GRU128,0.05s horizon,0..3bar,6bar/s slew,zero activation,
and0.02N frozen-balance diagnostic tolerance. The stored common1.6bar/trim/PAC
correction-cap entries remain compatibility slots, unused by this variant.

## Main pressure results

Counts/percentages use the valid inverse context rows, not the whole record.

| Quantity | Circle sweep | Helix |
|---|---:|---:|
| Duration / total rows | 375s / 7501 | 290s / 5801 |
| Valid inverse rows | 6961 | 5257 |
| Corrected maximum-pressure scaling | 64 (0.92%) | 40 (0.76%) |
| Corrected slew limiting | 774 (11.12%) | 859 (16.34%) |
| Nominal slew limiting | 4254 (61.11%) | 3334 (63.42%) |
| Corrected=nominal before slew | 0 | 0 |
| Corrected=nominal after slew | 3 | 0 |
| Max absolute correction before slew | 1.62788bar | 1.62110bar |
| Max absolute correction after slew | .60000bar | .60000bar |
| Zero actual-command valid rows | 979 | 852 |
| All-three-.3bar corrected outputs during zero command | 0 | 0 |

All corrected outputs remain within0..3bar and within6bar/s times the logged
horizon of the actual command at that tick. At least one pre-slew pressure is
zero on every valid sample, as required by minimum-nonnegative allocation.
Slew clipping may make the final vector have all nonzero entries when actual
commands are high; this does not reintroduce common-preload allocation.

First valid sample is1.65s in both recordings, actualcommands[0,0,0]bar:

| Prediction | Circle,[P1,P2,P3]bar | Helix,[P1,P2,P3]bar |
|---|---|---|
| Nominal | [0,.0189,.0114] | [0,.0043,.0043] |
| Corrected | [.0200,.0168,0] | [.0254,.0093,0] |

Exact values and component correction ranges are in case summary.json.
See each case pressure_startup.png and pressure_full.png, also PDF/FIG.
Left panels show actual commands/limited predictions; right panels show the
two allocations before slew. These pressures are commands, not measured chamber
pressures. Pressure agreement with the actual command is not independent
inverse accuracy: the actual command is also a network input and slew anchor.

## Source and wiring checks

Full replay of all6961+5257valid contexts using SAVED settings reproduces
pressure outputs within1.34e-15bar and target/residual/pre-slew diagnostic values
within1.07e-14. There are0valid/allocation/slew flag mismatches. The separate
algebraic minimum-pressure calculation and clip identities match<=4.45e-16bar.
Source-command age is0 on valid rows, sampled NI commands equal context pressure
exactly, and the inverse residual equals the forward observer residual exactly.
The recorded networkEvaluations is1 on every valid row.

The older PAC-policy single-pass source does not reproduce these outputs
(distributed-context max discrepancies~1.60bar), confirming the mapping change.
The separate old iterative inverse still matches the variable-common PAC
policy (~1.1e-8bar circle/~8e-12bar helix). It remains a DIFFERENT allocation
policy, so cross-block differences are not a matched iterative-versus-single-
pass solver performance test.

On valid rows, actual NI commands exactly equal the original collection
pressureRequested signal. This is consistent with a shadow observation test;
these recordings do not demonstrate physical residual-on PAC/PDFB tracking.

## Model equation consistency and limiting

Extra GRU evaluations at proposed pressure were done offline with immutable
history, not committed by either running observation block. All quantities
below are generalized-force MODEL equation mismatch, not measured force.

| Error norm,median/P95/max (N) | Circle | Helix |
|---|---|---|
| Frozen residual, final limited pressure | [~0,3.4459,7.7025] | [~0,2.6862,7.6932] |
| Residual reevaluated at final proposed pressure | [.06110,4.0416,8.4019] | [.10035,3.2014,8.1804] |
| Residual reevaluated at pre-slew allocation | [.06110,.9976,1.7398] | [.10035,.9381,1.7466] |
| Nominal final pressure, without residual | [6.8407,19.8236,32.3292] | [5.2662,18.4926,31.1834] |

Of735circle samples with final model mismatch>1N,733coincide with scaling or
slew limiting; the helix counts are750/751. Without either corrected-pressure
limit, model mismatch median is.04957/.06838N. Limiting dominates the larger
error tail here. The residual still changes with pressure, so zero frozen
balance does not mean exact nonlinear inversion: only1154/6961circle and
1321/5257helix rows meet the0.02N criterion after candidate-pressure reevaluation.
This distinction matters for any later PAC adaptation gate.

The low median errors support the numerical single-pass approximation on
much of these datasets, not a uniform bound, pressure-model identification,
passivity guarantee or closed-loop qualification. The target remains measured
q/v with supplied desired acceleration; it is not the full desired-state
collection inverse demand or a feedback controller.

## NDI and readiness

| Quantity | Circle | Helix |
|---|---:|---:|
| Orientation-invalid / tip-invalid | 16 | 18 |
| Those with raw tip MISSING status | 16 | 18 |
| Base-invalid | 0 | 0 |
| Valid orientation but theta unobservable | 0 | 1119 |
| Logged acquisition overruns | 3 | 1 |
| Forward numerical-guard status rows | 0 | 0 |

Every invalid-orientation sample is rejected by the observer. All1119helix
theta-unobservable samples have valid orientation and are accepted; their
inferred bends range.0000529–.0086535rad, below0.5degree. These are near-straight
direction ambiguities, not additional tracker losses or neural resets.
Raw MISSING identifies the rejection stage, not a unique physical cause.
Actual losses retain the established segment reset/warmup behavior. Do not
interpret held measured displays as fresh data through these losses.

## Forward length prediction check

The raw NDI/observer q identity, logged-error recalculation and strictly prior
issued-to-matched prediction identity all match EXACTLY. Valid scored coverage
is6927/7501=92.35% and5236/5801=90.26%, with gaps retained.

| 50ms model prediction RMSE,[l2,l3]mm | Circle | Helix |
|---|---|---|
| Nominal | [4.4070,4.8237] | [4.5256,4.6504] |
| Residual+nominal | [.8482,1.2829] | [1.0855,1.2813] |
| Constant velocity | [.0997,.1138] | [.1242,.1370] |

Residual augmentation still improves the nominal predictor substantially,
while CV remains better at50ms. This observation matches earlier cohorts;
different physical runs are not a causal test of the allocation change.
One-step predictions are measurement-initialized, not free whole-motion
propagation, physical pressure validation or feedback tracking error.

## Interpretation / next work

The user-installed collection-policy block has now been exercised in an
actual model run and independently reproduced. There is no implementation
fix indicated by these checks. The remaining issues to carry into the PAC
math are bounded pressure delivery, nonzero pressure-dependent residual
change, adaptation gating and explicit recovery from true missing samples.
Do not equate the frozen-balance flag with true/model force delivery or
assume a matched comparison against the still-PAC-policy iterative block.

Retained percase CSV/MAT/JSON, full/startup PNG/PDF/FIG, root summary and
forward metrics, saved settings, source/raw/model hashes and independent CSV
checks. Reproduce with analyzeCollectionShadowRuns(root,freshFolder), no model
run required. v1 is an incomplete analyzer attempt at ambiguous overrun names;
v2 explicitly selects the NDI MATLAB System diagnostic and is canonical.
