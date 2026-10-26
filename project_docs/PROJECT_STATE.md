Last Updated: 2026-10-23
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

Branch: codex/passivity-adaptive-controller. Verified latest commit: 4c9f930
(100 g payload comparison documentation); no remote push was performed.
Do not commit/push or run hardware without the user's request.
Do not edit SLX models unless requested; the user installs pasted block code.

- passivity_adaptive_control/pacAdaptiveTunableBlock.m: current tunable wrapper.
  First11 inputs/26 outputs preserve the previous adaptive wrapper interface;
  adds common-pressure Boolean, ten tuning vectors, mean-pressure diagnostics.
- pacTuningBlock.m: one no-input MATLAB Function supplying all ten vectors.
- Tunable wrapper appends input23 freezeHysteresisWeights: true freezes only
  h1/h2/h3 estimates; history states/forces keep evolving, other adaptation
  keeps its existing gates. False resumes; default false for old MATLAB calls.
  Baselinefalse vs steptrue at120s tests freezing after10x10s circles.
  Offline freeze/resume/global-disable and existing validation checks passed;
  no SLX edits or hardware tests. See TUNABLE_BLOCK.md for paste/wiring.
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

- Trial3 freeze recording verified: new switchtrue121.11s,PAC121.15s;h held
  exactly1.53026/1.35607/1.32159 while z/otheradaptationcontinue. Historyforce
  amplitude and mean stiffness trends plateau. Final-lap thetaMAE.520,phi.150,
  lengths4.281/4.900mm vs trial2continuousadapt .416/.102/3.520/4.026.
  Supports useful continuing history compensation,not physical identification;
  prefreeze differences/unlogged tuning/preconditioning limit causal claims.
  58invalid-tip single-sample events,no faults/limiting. See VALIDATION_RESULTS.

- Trial2 constant-gain30x10s reviewed: lambda18/16.5,Ks23/10. History weights
  keep increasing while z stays periodic; mean spring coefficients decline.
  Final-lap thetaMAE.416,phi.102,lengthRMSE3.520/4.026mm.48 single-sample
  invalid-tip events (47 in circles) zero force-evaluation outputs,hold pressure,
  and freeze stored estimates; no faults/limiting. adaptiveForce diagnostic is
  misleading on invalid samples due to unconditional nominal subtraction;
  mask by modelValid pending a deliberate fix. Tuning vectors not found in MAT
  despite intended logging; observed k2/k3 rate250 differs from source cap200.
  See VALIDATION_RESULTS for mechanism,limits and exact evidence.

- Later 1rad_10s_tr1.mat reviewed: thirty10s circles,phi1,live gain changes
  (lambda initially21/15 finally21/16; Ks initially25/20 finally25/30).
  Final-lap thetaMAE .326,phi .089,lengthRMSE2.693/3.474mm; no faults,
  almost no limiting. 55 zero diagnostic samples are invalid-data holds, not
  estimate resets. Live tuning confounds attribution to adaptation. See
  VALIDATION_RESULTS for exact changes; gamma/rate-cap inputs not logged.

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

plotPacPresentation(filePath,timeWindow,outputDir,prefix) exports presentation
PNG/PDF/SVG/FIG from saved logs only. 100g10s trial2 figures include full0-330s
and final3cycles290-320s; sensor-coordinate XYZ reference includes base/tip
offsets. See recording's presentation_1rad_10s_tr2_100g/README.md for conventions
and actual live gain changes. No hardware/model execution for plotting.

10 s circle hysteresis audit: history states remain near +/- .666 while adaptive
history weights grow to ~3.1-3.4 N; force growth is weight drift, not history
state explosion. No convergence established; missing dynamics/parameter tradeoffs
are hypotheses. See VALIDATION_RESULTS before increasing history gamma.

Read-only 2026-10-16 audit of unloaded 20 s circle: 31 zero damping-output
samples all match invalid model/status 4; stored damping estimates do not reset.
Spring/bias observed rates reach current source caps; gamma/rate-cap inputs
were not logged. See VALIDATION_RESULTS for limits and measured peak rates.

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
5. Integral PAC action discussed as a proposal only: existing adaptive biases
   already integrate filtered error subject to gates/bounds. Any extra integral
   action needs anti-windup and coordination with bias adaptation plus revised
   stability analysis. Separate static-offset correction from directional lag;
   actuator dynamics remain the first proposed lag investigation. See
   TECHNICAL_NOTES for options and adaptation-speed tuning distinctions.

## Documentation maintenance

PROJECT_LOG.md is append-only history. Update this state in place and relevant
MODEL_ARCHITECTURE, VALIDATION_RESULTS, OPEN_ISSUES, DECISIONS or TECHNICAL_NOTES
when code or important conclusions change. Root PROJECT_LOG.md is a pointer.
Use local READMEs and TUNABLE_BLOCK.md for integration, not historical settings.

100g20scircle reviewed: same28/5 gains,no faults/minimallimits,lap10thetaMAE.311 vsunloaded.288,phi.0636 vs.0642. Payloadtiming/modelchangeunverified; see validation. No payloadmass estimator.

## 2026-10-19 - Current additional task
Read-only RL/open-loop acquisition audit completed for Dynamic_model_Val_FWDnINV_v2.slx. See RL_DATA_LOGGING_INSPECTION.md. No model edits; settle signal semantics/duplicates, raw orientation retention and trial metadata before bulk data. NDI20Hz, pressure outputs100Hz, commands not measured pressure.

First RL trial forward_test_RL_30s.mat validated offline: P1=1.5*(1-cos(.2*t))bar,600/601NDIvalid,XYZcorrectionexact. Hostelapsed31.48s vs30sim,10overruns; preserveactualtimestamps. See validation. One partial sinecycle only; no embedded trialmanifest.

Open-loop pressure sources now in RL_data_collection/: Steps, Sine, Combined, with trial labels and complete endings. Start from previously tested single-P1 sine; see README for wiring/stop times. Raw commands for complete dynamics learning, no source deadzone compensation. No model edits.

RL_data_collection adds desired-angle circle sweep(.2,.8,1,1.2,1.5;5x20s each) and4loading/unloadinghelixtrials. These requireopenloopfeedforwardpressuremapping; notdirectvalveinputs. Sequentialpressuredefaults nowall3channels. See README.

RL CSV/report exporter exportRlTrainingData.m nowavailable;actualSineWave.mat exportedto RL_data_collection/exports/SineWave_reviewed and ZIP. DATASET_HANDOFF.md explainsrealactions/observations vsmodelpredictions. 12710/12721valid,hostgapsflagged;metadata stillrequirescompletion. No modelchanges.

Latest RL_dataCollection.slx feedforwardaudit:angles->desiredlength/rates->inverse->valves wired,deadzone0confirmed. CriticalinversevalidityinputFrom3/A1 is measuredtheta;reconnectpropervalidity. Stop342sneedscircle625. Helixunconnected. See inspectionlog; no model edits.
