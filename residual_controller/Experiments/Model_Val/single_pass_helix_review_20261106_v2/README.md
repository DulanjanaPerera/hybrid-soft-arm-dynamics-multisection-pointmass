# Operator single-pass helix pressure audit — 2026-11-06

Input: helix_FWD_INV_singlePass_test.mat, out-only, model
residual_FWD_INV_validation1.52,290s/5801rows, ReachedStopTime, no saved model
error. There are5266 jointly valid context/iterative/single-pass samples.
Offline inspection/numerical replay only. No SLX/runtime/NN/protection edits,
hardware/model/acquisition run, training, commit or push.

## The0.3bar startup is a shadow-output artifact

The first valid prediction is at1.60s. The actual NI commands are[0,0,0]bar.

| Quantity | P1 | P2 | P3 | Units |
|---|---:|---:|---:|---|
| Actual command | 0 | 0 | 0 | bar |
| Nominal allocation before slew | 1.59348 | 1.62259 | 1.58392 | bar |
| Single-pass allocation before slew | 1.61134 | 1.60024 | 1.58842 | bar |
| Nominal, iterative and single-pass outputs, each | .3 | .3 | .3 | bar |

The reused PAC policy chooses common pressure near1.6bar. Both diagnostic
inverse wrappers then clip around the **actual sampled NI command**, with
allowance6bar/s*.05s=.30bar. The stateless shadow output cannot accumulate
its own hypothetical ramp while actual command remainszero. All853valid
zero-command samples produce[.3,.3,.3] in allthree inverse outputs. This is
neither a learned deadzone nor a0.3bar physical pressure minimum, and it does
not show the actual NI commands starting at0.3bar. See pressure_startup.png.

## Why the two inverse curves look similar

- Single-pass and iterative pressure vectors coincide3960/5266=75.20% of
  valid samples (max-axis equality tolerance1e-8bar).
- Single-pass slew limiting5099/5266=96.83%; nominal slew limiting5266/5266.
- Single-pass equals nominal3570/5266=67.79%; iterative equals nominal3647/5266.
- Single-pass differs from nominal at **every** valid pre-slew sample;
  maxaxis pre-slew correction reaches0.97574bar.
- Single-pass versus iterative maximum-axis difference: median0bar,
  P95.28666bar,max.60000bar. They are not identical everywhere.

Both wrappers use exactly the same nominal force target and same actual
command context. The single-pass GRU force exactly matches the forward
observer's force. Nominal outputs from the two blocks match exactly, and
the logged clipping formula reproduces single-pass and nominal outputs exactly.
No evidence of a disabled residual or wrong single-pass wiring was found in
these checks. The comparison is strongly masked by the shared pressure policy.

The older inverse does an **iterative nonlinear pressure search**, not ODE
integration. Forward prediction separately integrates the equations of motion.

## Training versus pressure allocation

Training used recorded commanded pressures as forward-model inputs and learned
two generalized-force residuals through length-prediction loss. It did not
contain this PAC inverse allocation, its1.6bar nominal common pressure, or the
shadow0.30bar step restriction. These were inverse-side implementation choices.

Converting two generalized forces to three pressures always requires a
convention. The saved collection inverse chart currently uses the minimum
nonnegative solution, then scales to3bar if necessary, with zero deadzone:

```
d = requestedForce / forcePerBar
p1 = max(0,-d1,-d2)
p = [p1,p1+d1,p1+d2]
if max(p)>3: p = p*3/max(p)
```

This is an implicit allocation policy even though it is not a standalone
allocator block. The saved SLX was inspected as ZIP/XML only, never loaded or
executed; unsaved operator chart changes cannot be inferred from this snapshot.
The chart source and SLX hash are retained. Its original demand uses supplied
desired q/dq/ddq, whereas the shadow demand uses measuredq/v plus desiredddq;
the counterfactual below holds that shadow demand fixed to isolate allocation.

The supplied shadow wrappers reused a DIFFERENT PAC convention: common
pressure near1.6bar plus correction/rail limits and slew. This choice should
have been kept explicit for the first comparison with the collection inverse;
it explains the artificial preload seen in the observation output. It is not
a property the GRU learned or something required by residual-force training.

Offline copying the collection minimum-pressure convention on the **same
first valid force targets** gives:

- Nominal:[.009561,.038671,0]bar.
- Single-pass:[.022919,.011821,0]bar.

No1.6bar common offset is added. On5186samples where both mappings preserve
the requested frozen force, their B*p difference is<=7.11e-15N: common pressure
accounts for the remaining difference in pressure vectors. The alternative
scales30corrected targets to3bar. It is retained in
offline_collection_allocation.csv/json/mat as a calculation ONLY, not applied
to a live block. The residual itself depends on absolute candidate pressures,
so equal B*p does not prove equal augmented dynamics or physical response.

## Single-pass approximation remains measurable

Additional GRU calls at the new predicted pressure were performed OFFLINE with
immutable history, not committed by the running block. Generalized-force
equation-error norm statistics (not measured physical force) are:

| Model equation mismatch | Median | P95 | Max | Units |
|---|---:|---:|---:|---|
| Single-pass, residual frozen at actual pressure | .46364 | 4.94911 | 9.53669 | N |
| Single-pass, residual reevaluated at predicted pressure | .60208 | 5.28348 | 9.90478 | N |
| Iterative inverse | .49672 | 4.06277 | 8.07370 | N |
| Nominal inverse | 8.03121 | 19.55838 | 24.09891 | N |

Changing from actual to predicted pressure changes residual force norm by
median.25187N,P95.66319N,max1.19395N. This is why frozen-residual subtraction is
an approximation. Iterative converged at241/5266samples under0.02N criterion;
the existing constrained local method also has previously demonstrated
projection plateaus. The current post-limit curves do not isolate algorithm
quality, physical inverse accuracy or controller benefit.

## Recorded inputs and replay

Only out was saved; complete option vectors were not saved. Source-default
variable-common replay over120distributed valid contexts matches all logged
single-pass pre/post outputs within8.89e-16bar and iterative outputs within
7.27e-9bar. Fixed-common settings do not reproduce the run (maxerrors.225/.374bar).
This supports the default-policy explanation without uniquely identifying
every unlogged option. Targets match exactly in replay. No timing/gain changes
were inferred from filenames alone.

Native NI command signals use a faster grid than the0.05s observer. The audit
aligns each NI signal to t-commandAge using bounded previous sampling(.075s),
then checks the context. Recorded commandAge iszero on valid rows here, and
both current and timestamped NI commands exactly match context pressure.
There is no evidence of a source-pressure mismatch in this recording.

## Next recommendation

For the next **passive observation** comparison, use an explicit allocation
choice consistent with the collection inverse and expose nominal/corrected
values before rate clipping. Keep actual commands, hypothetical allocated
pressures and protected output separate. Do not interpret the clipped0.3bar
trace as learned compensation. No live pressure protection needs to be removed
to inspect the already logged pre-slew values. Matching allocation convention
does not retrain the network or establish closed-loop readiness.

No change to either runtime inverse, PAC/PDFB, source settings or SLX is applied
in this review. The proposed policy change remains next work.

## Artifacts

summary.json,first_valid_samples.csv,pressure_comparison.csv,analysis.mat,
signal_inventory.csv; full/startup PNG/PDF/FIG; offline collection-policy
counterfactual CSV/JSON/MAT; input/source/artifact SHA256 and source snapshots.
analyzeSinglePassHelixLog and compareCollectionShadowAllocation reproduce the
analysis without loading a model. The emptyv1attempt is marked incomplete:
an analyzer assumption of identical NI/observer grids failed; v2 handles the
native rates correctly. This was not a runtime block failure.
