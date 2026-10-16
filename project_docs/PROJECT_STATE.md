Last Updated: 2026-10-16
Current Milestone: Tunable adaptive PAC with variable common-pressure allocation implemented; initial physical circle data reviewed; directional lag remains.

# Current project state

Read this first, then the relevant specialized document. Use recorded inputs
for experiment settings; do not infer them from filenames or current SLX defaults.

## Objective and mechanics

Control a backbone-constrained pneumatic continuum arm with paired NDI poses.
Mechanical coordinates q=[l2;l3] (m), l1=-l2-l3. Three pressure commands are
actuator inputs, not three mechanical coordinates. NDI lengths are inferred
from orientation and constant-curvature geometry, not independent strain sensors.
Nominal active length .17411 m, radius .013 m, mass .1 kg; actuator area radius
.0065 m. Model/pressure-map values are assumptions, not uniquely identified truth.
Gravity is transformed using the base orientation. Nominal sensor period .05 s.

## Code and integration

Branch: codex/passivity-adaptive-controller. Last user-requested commit:433750a;
subsequent common-pressure/tuning changes and experiments are uncommitted.
Do not commit/push or run hardware without the user's request.
Do not edit SLX models unless requested; the user installs pasted block code.

- passivity_adaptive_control/pacAdaptiveTunableBlock.m: current tunable wrapper.
  First11 inputs/26 outputs preserve the previous adaptive wrapper interface;
  adds common-pressure Boolean, ten tuning vectors, mean-pressure diagnostics.
- pacTuningBlock.m: one no-input MATLAB Function supplying all ten vectors.
- TUNABLE_BLOCK.md: exact ports, defaults, dimensions and integration steps.
- pacAdaptiveCommonPressureBlock.m: retained older fixed-mean wrapper.
- Parameter order: k1 k2 k3 d2 d3 h1 h2 h3 b2 b3. K is structured SPD under
  nonnegative spring coefficients and positive k1/k2; D diagonal. History
  weights/bias adapt; internal normalized history law is separately configured.
- Lambda and Ks are externally selected feedback gains, not adaptive estimates.
  Model damping d2/d3 is distinct from Ks. Mass/gravity remain nominal.
- Explicit .05 s PAC execution is the intended setup. Valve outputs may run .01 s.
  Acquisition slowdown once cleared after restarting MATLAB; profiling preflight
  failed before PAC execution, so it did not establish PAC computational cost.

## Common-pressure allocation

OFF: retained activation-aware fixed mean1.6bar allocation with bounded fallback.
ON: all-active nominal force projection onto a feasible pressure-difference
polygon; choose feasible mean nearest nominal. Requests outside feasible set
receive closest Euclidean nominal force. Rails/correction caps retained.
Defaults: pressure .2-3bar, activation .8bar, s/a correction caps2.5bar,
slew2bar/s. ON commands remain at least activation (or configured minimum).
Common-pressure selection is allocation, not a gradient-estimated model parameter.
Trim/anchor, freshness, invalid-frame hold and fault supervisor remain.
Structural tuning/calibration/initial values should be set before a fresh run.
Exact nominal force balance is not measured physical force delivery. Changing
common pressure may change stiffness/preload outside the ideal nullspace model.

## Latest validated evidence

- Offline allocator, adaptive algebra/supervisor, variable-common projection,
  switching, custom initial-value and MATLAB Coder checks passed. These do not
  establish physical sampled-loop stability or passivity.
- Ks40/lambda50 produced recurrent zero-return oscillation. Ks20 and reduced
  lambda suppressed it. Lambda40/Ks20 was a useful fixed-mean static baseline.
  k3 hitting zero was NOT the onset: it stayed zero during quiet and oscillatory
  intervals while K remained SPD. Negative k3 is not authorized or implemented.
- Fixed-mean circle allocation left force residuals that gated adaptation;
  activation-consistent fixed-mean update did not improve tracking because many
  requested forces were infeasible at mean1.6bar. Variable mean addresses this.
- Latest reviewed file:
  passivity_adaptive_control/experimentations/adaptCommonP/1rad_20s_tr1.mat
  Ten20s circles,phi1rad. ACTUAL logged gains lambda28/Ks5 (not40/20).
  Mean pressure1.597-1.911bar,commands .8-2.803bar. Most nominal residuals
  roundoff; brief lap7 transient. Adaptation active97-99.5% per lap.
  Theta MAE lap2 .436 -> lap10 .288rad; phiMAE .072 -> .064;
  lengthRMSE3.52/4.03 ->2.54/2.72mm. Dominant positive theta error indicates lag.
  Changed speed/gains confound direct fixed-vs-variable comparison.
- User reports later better gains; those recordings/settings have NOT yet been
  supplied/analyzed. Parameters cycle with direction and some hit bounds;
  do not interpret as identified physical parameters or proven convergence.

## Analysis workflow

analyzePacAdaptive() and analyzePacStaticSweep() use base-workspace out; explicit
out and MAT paths also supported. Workspace calls do not automatically save.
Static-sweep analysis is for holds, not circle tracking. Static lambda inference
is invalid for circles with nonzero desired rates. Circle metrics so far were
computed offline per lap; no dedicated current PAC-circle report script added.
plotCommonPressureTracking() supports PAC/FBL signal names and duplicates use
first matching name. Theta plots/errors wrap to[-pi,pi]; exclude near-straight
unobservable theta, especially final return. Missing essential signals still fail.

## Next decisions / open issues

1. Review later better-gain data or matched ON/OFF circle test at same speed,
   gains, initialization and reference. Do not increase gains/speed blindly.
2. Quantify signed phase error, actuator response and allocation constraints by
   direction. Nominal residual alone cannot identify physical force uncertainty.
3. Uncertainty-model discussion: favor simple actuator dynamics first, then
   bounded residual compensation if justified. No observer, lag model, robust
   term, neural model or revised stability proof has been implemented.
4. Common-pressure dependence of stiffness/hysteresis and pressure-to-force
   mismatch remain unmodeled. No chamber pressure sensors are available.

## Documentation maintenance

PROJECT_LOG.md is append-only history. Update this state in place and relevant
MODEL_ARCHITECTURE, VALIDATION_RESULTS, OPEN_ISSUES, DECISIONS or TECHNICAL_NOTES
when code or important conclusions change. Root PROJECT_LOG.md is a pointer.
Use local READMEs and TUNABLE_BLOCK.md for integration, not historical settings.

