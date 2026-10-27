Last Updated: 2026-10-27
Current Milestone: Tunable adaptive PAC with variable common-pressure allocation implemented; initial physical circle data reviewed; directional lag remains.

# Current project state

2026-10-26 framework development: branch `codex/pac-autotuning-framework`.
Milestone A offline evaluator is implemented in `passivity_adaptive_control/autotuning`;
see its README and `PAC_AUTOTUNING_FRAMEWORK_PLAN.md`. Offline tests reproduce
v8 lap errors, v9 slew/gate counts and v10 pre-fault roughness/fault separation.
Reviewed reports are in autotuning/reports/historical_reviewed_v1. No targets
have been agreed and no hardware automation exists.

2026-10-27: milestone A committed as e6fb02c at user request. Milestone B is now
implemented: immutable campaigns, bounded manual gamma proposals, reviewed tuning
snapshots, operator manifests, offline imports, append-only event ledger, restart,
budget/duplicate checks and separate screening/confirmation status. See
autotuning/MANUAL_WORKFLOW.md and reports/manual_demo_v1 (historical replay only,
not configured for success). Offline tests pass, including failed-confirmation
withdrawal and longer confirmation windows.

Milestone C is now implemented (2026-10-27): proposal-only bounded GP optimizer
in normalized log gamma space, with seeded initial design, noisy/repeated
observations, explicit failure/missing-data separation, held-out confirmation,
restartable ledger-derived proposals and stale-suggestion checks when queueing.
See autotuning/OPTIMIZER_WORKFLOW.md and reports/optimizer_synthetic_v1.
Core/adapter/manual/evaluator offline tests pass. Local bayesopt/Statistics are
available; this adapter uses base MATLAB numerics and invokes no trial callback.
Milestones B/C are included in the user-requested framework checkpoint; use git
log for its revision. Next is milestone D, an operator-run campaign;
real baseline, bounds, targets, conditions and confirmation protocol remain
unsettled. No hardware run, SLX edit, automatic gain write or push performed.

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

2026-10-26 user-requested experiment: pacTuningBlock lower(2)=-100 N/m;
pacTuningConfig now permits lower(2)>=-500 after subsequent user experiments,
while the source block default remains -100. The helper also accepts k1/k3
lower bounds >=-100; source defaults remain [50;-100;0] for k1:k3.
Legacy defaults stay positive. Negative k2 is signed compensation, not identified
physical stiffness. Relaxed bounds no longer guarantee SPD K; current update
does not enforce coupled SPD and modelValid does not check it. Monitor
k1*k2+k1*k3+k2*k3>0 with k1+k3>0. No SLX edits/hardware runs; see validation.

Branch: codex/pac-autotuning-framework (created from 63f679a). Use git log for the latest checkpoint;
the checkpoint preceding this review was bd83708. No remote push was performed.
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

- v10 actual4scircles,slew6,caps850/20/5/.2,initialgamma_k.5,d.02,h5.
  Growingoscillationnear55precedesgammachanges78.7/93.5;estimatedd~86/90by55,
  ~104/109by78.7;loweringgammadoesnotresetestimates. Kindefinitecontinuously
  from55.2(notproofcause);k2boundfirst62s. Meanerrorfallswhileroughnessgrows.
  Finalshutdown104.85isseparatetip-invalidgracefault;excludeafterfaultfrom
  trackingmetrics. Seev10onsetaudit;nohardware/controller/modelchanges.

- v9 verifiesactual2scircles20-80s(phi1),despite10sfilename. Trackingpoor
  (finallaptheta1.067,phi.275,length8.07/9.76mm). Adaptactive5.42%;89.25%
  sampleslimited. Exactnominalallocatorreplay:all1071limitsareslew mismatch,
  no allocationinfeasibility;2bar/scapdominates,beforepressureupperrails.
  k3/h1/h3spendmosttimeatlowerbounds,butfeweligibleupdatesdominateproblem.
  Needactuator/commandslewheadroomassessmentratherthangammaincrease;
  donotdisablegates. No faults;estimatedKSPDthisrun. Seev9allocationaudit.

- v8 saturationaudit:fixedfeedback18/16.5,15/10;hgamma5.8,k.7->.5at221.65,
  d.002->.02at271.45. Theta<.2bylap8;lateplateau~.1rad/phi.062. k2lower
  24samples,k3lower93,k3ratecap25early;nofinal3cyclevalue/ratecaps. Active96.15%,
  nopressurelimiting. 2sgoalrequires5xreferencevelocity/25xacceleration;simple
  compressionofloggedpressurewaveformneeds~7-10bar/svs2limit(not2sprediction).
  Bandwidth/lag/SPDissuesneedattention,fasteradaptationaloneinsufficient.
  Detailedartifacts:v8_analysis/bound_and_rate_audit.txt;seeVALIDATION_RESULTS.

- v7 livek/h/dgammaaudit:feedback18/16.5,15/10fixed;hgamma5.8coincideswith
  fasterlearning,finallaptheta.064,phi.051,length.787/1.013mm. Butgamma_d
  .02->2at254.95 drivesestimatedd~50->107/118by272.6;late~112/124.
  ThetaSD/nearstop/phiovershoot/limitingworsen;loweringgammadoesnotresetd.
  Do not usefinalMAEaloneasgainselection. Proposedfreshfixedcandidate
  gamma_k.5,h5.8,d.002,bias1;notvalidatedoptimum. Springcaps650,hupper20;
  signedKindefiniteperiodsstilloccur. SeeVALIDATION_RESULTS/v7analysis.

- v6 verifies springgamma.8,fixedlambda18/16.5,Ks15/10,historygamma1.5,
  hupper20. No cleartrackinggainvsv5;springratecapsnowhit10-16%updates,
  estimatedKindefinite12.9%circlesamples(v5.7%). Nofaults/limiting,but
  moreaggressivesignedstiffnessadaptationnotfavoredbydata. Recommendfixed
  feedbackbaselinewithspringgamma.4forcleancomparison;nocodechangesmade.
  SeeVALIDATION_RESULTS/v6auditforconfoundsandmetrics.

- v5 livefeedback endsat88.2s withlambda18/16.5,Ks15/10;historygamma1.5
  unchanged. Latevsv4:thetaRatevariability37%lower,nearstopfraction3.2v14.6%,
  l3peakvelocity11.8v15.7mm/s,nopressurelimiting. Morelag:errorlast3theta.196v
  .116,phi.064v.053,length1.927/2.263v1.211/1.453mm;continuesimprovingafter
  finalgainchange. Smoothbaselinepromising,butnotgammacomparison. Nextfixed
  finalgainrunproposedbeforemoregammaadjustment. SeeVALIDATION_RESULTS/v5audit.

- v4 secondKs.1->1 verified,otherloggedtuningmatchesv3:smalllateerrorgain
  (last3meanlaptheta.116v.122,phi~.053both),l3velocityRMSE~5%lower,but
  theta near-stopfractionstill14.6%;roughnesspersists. Adaptactive96%;history
  weightsnear9.26of10;no faults. No clear smoothness cure fromthisKsstep.
  Nextinspectfilters/actuatorresponseonfixedbaselinebeforemoregainsweeps.
  SignedKindefiniteearlycycles/return;seeVALIDATION_RESULTSandv4analysis.

- v3 historygamma1->1.5,otherloggedtuningfixedversusv2:finalthetaMAE.119v.168,
  phi.053v.054,length1.200/1.543v1.643/1.881mm;theta<.2bylap17v24.
  Adaptactive96.02%. Historyweightsapproach9.25of10bound. Within-cycle l3
  velocitybursts/thetaslowdownsremain;latephiovershootmax1.220v1.168.
  Fasterlearningnotacureforroughness;nextseparatedamping/actuatordynamics
  comparisonproposed,nocodechange. SignedKindefiniteearlycycles/return.
  SeeVALIDATION_RESULTSandv3analysisforcomparisonandlimitations.

- Fixed-gain1rad_10s_gainTesting_v2:lambda28/25.5,Ks25/.1 unchanged;thetaMAE
  lap2.435->lap30.168,phi.139->.054,length3.989/4.191->1.643/1.881mm.
  Adaptactive96.18%;off105limited+69invalid+55dt,zero deadbandcases.
  History/dampingupdatesnotratecapped;furthercapincreaseunlikelytoacceleratethem.
  Historyweightsstillgrow(~7Nlate),slowergradientlearningplausible;notproof
  of convergenceorphysicalidentification. SignedK indefiniteearlycycles despite
  nofaults. SeeVALIDATION_RESULTS forsettings,limitsandofflineartifacts.

- 2026-10-26 rad_10s_gainTesting: all tuning inputs logged. Circleadaptation
  active92.23%;off466samples=358limiting+67invalid+41dt,zero deadband cases;
  gate reconstruction exact. Kssecond100 nearly blocks adaptation (0.8%active)
  and worsens tracking. Ratecaps rarelyreached;gamma/rate increases not main
  fix indicated. Finallaptheta.188,phi.060,length1.967/1.904mm,live gains
  confoundcausal comparisons. SignedK indefinite23.95-25.35s and finalreturn,
  despite nofaults;currentmodelValiddoesnotcheckSPD. SeeVALIDATION_RESULTS
  andadaptCommonP/rad_10s_gainTesting_analysis/README.md foractualsettings.

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
computed offline per lap; the autotuning evaluator now provides reusable circle
reports and v8/v9/v10 replay checks.
plotCommonPressureTracking() supports PAC/FBL signal names and duplicates use
first matching name. Theta plots/errors wrap to[-pi,pi]; exclude near-straight
unobservable theta, especially final return. Missing essential signals still fail.

## Next decisions / open issues

2026-10-27 tuning framework: evaluator, manual ledger and bounded grouped-gamma
proposal adapter are implemented offline. Feedback gains remain fixed campaign
settings. The user starts hardware runs. Score tracking plus roughness/limiting,
use defined adaptive initial states, and confirm over repeated cycles. Two circles may screen
candidates but do not establish learned performance. Address the 2s command-slew
bottleneck and signed-stiffness constraint before autonomous exploration.
No automatic hardware runner or unattended physical tuning is implemented or
authorized here. Optimizer proposals alone do not establish stable operation.

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

2026-10-26: pacTuningConfig accepts inclusive experimental stiffness lower bounds [-100;-500;-100] for k1/k2/k3. Only the assertion changed; live tuning vectors must supply the desired bounds. No Simulink edits or hardware runs. Signed coefficients do not guarantee positive-definite estimated stiffness. Offline signed-bound projection/recovery validation passed for all three coefficients.
