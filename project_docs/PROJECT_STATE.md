Last Updated: 2026-10-28
Current Milestone: Tunable adaptive PAC with variable common-pressure allocation implemented; initial physical circle data reviewed; directional lag remains.

# Current project state

2026-10-28 R033: three5s circles completed29.6-44.6, no tracking limiting or
latched PAC fault; constant gamma. Last2circle exploratory lengths[4.6818,
6.7801]mm, velocity.012102m/s, roughness.030820m/s^2. Coverage.884422<.95,
zero accepted GP observations despite one completed-motion batch. Score200
samples:192rawvalid,184controlled,176/199pairs; eight tip-validity dropouts,
base valid/overrun0. Raw-only pairs183/199 still<95%; do not loosen posthoc.
Need inspect tip tracking/frame quality before gain search. Rejection-ID bug
fixed: reply1 while runtime waited2 caused30s hold thenreason11. Now errors
reply2, don't waitstationary, retain invalid attempt/checkpoint separately.
Recorded regression and runtime/validity suites pass. Actual commandszero86.9,
stop87.9, allDAQzero. Detailed R03/review_r033_v1. No gains/policy/sensor changes
or hardware runs during review; existing session logs preserved.

2026-10-28 R03 reviewed: no circle commanded. Phase0 until startup timeout
50.05s reason9; all non-heartbeat gates ready20.55s, actual phi/theta0,
PAC status3, velocities <.0005m/s. Zero commands at shutdown, end56.30s.
No accepted worker heartbeat. Current netstat showed model32151 and worker32152
on same MATLAB backend, while model cache PID differed: stale bridge in worker
instance is supported conflict explanation. Worker now clears its model bridge
with models stopped; UDP bind failure explicit reason13/warning. Added
pacOnlineConnectionCheck and r03/check_connection.m pressure-free preflight;
must PASS in model MATLAB before retry. Alternate-port regression seeded stale
socket, tested cleanup/probes11s apart and synthetic completion. R03 recording
preserved; same fixed-gain smoke can retry after worker restart and preflight,
using new recording filename. No gain tuning or hardware run by assistant.

2026-10-28 r02 physical smoke reviewed: startup/worker/baseline acknowledgement
worked, constant-speed run29.6s, abort34.0s reason4 after0.88circles; zero
completed batches. No tracking limiting or latched PAC faults. Shutdown worked:
phase10 at39.0s, zero commands39.3s, finalDAQzero, stop40.3s. False-positive
dropout grace joined raw invalid33.90/valid33.95/invalid34.00 via delayed PAC
diagnostics, with strict floating-point >.10 boundary. Corrected raw-current
validity grace and 1e-9 tolerance, preserving PAC fault/limit/heartbeat handling;
invalid hold frames wait, fresh movement still aborts. Recorded-prefix and
synthetic sustained-invalid tests pass; not a counterfactual completion claim.
Review artifacts under smoke_phi1_T5_scale5_r02/review_v1. Fresh r03 fixed-gain
three-circle model/session prepared; no gain search or extra hardware run by
assistant. Repeat monitored smoke with same gains before optimization.

2026-10-28 smoke build-path correction: r02 compilation hit Windows 260-character
limit because cwd was the nested session folder. Added pacOnlineConfigureBuild
and called it from all generated open_model scripts: model-process cache/code
folders now go to short, PID-isolated local TEMP/pac_bld directories. No model
execution, project relocation or cache/log deletion. Existing listening r02
worker can remain; retry requires setting this in model MATLAB, not worker.

2026-10-28 smoke startup correction: user reported Desktop Real-Time rejects
simulation pacing, then worker zero packets/batches and overwrite refusal.
Installer now sets EnablePacing=off; generated initial/r01 model files corrected.
Fresh PAC_autotuning_online_smoke_r02.slx and smoke_phi1_T5_scale5_r02 session
prepared; use matching r02 worker, not an old loaded model. Worker initial
connection wait increased from five to thirty minutes; live heartbeat watchdog
unchanged. Old attempt/log preserved. No hardware run performed by assistant.

2026-10-28 next hardware step: prepared separate fixed-gain smoke session
autotuning/online_sessions/smoke_phi1_T5_scale5_r01 and new
PAC_autotuning_online_smoke.slx. Three 5 s circles, score last two, gamma fixed
at [.0515;.0005999;1.49999] by equal lower/upper bounds, budget one, no gain
search. Operator monitored only; save full out log and worker audit. Expected
normal terminal reason23 is budget exhausted, not target confirmation. No
hardware execution performed; inspect smoke logs before multi-candidate runs.

2026-10-28 continuous tuning implementation: created PAC_autotuning_online.slx
as a NEW wired copy of saved PAC_autotuning_manual_v2.slx, standalone block
library and session MAT under autotuning/online_sessions/phi1_T5_scale5_initial.
Starting gamma [.0515;.0005999;1.49999], 5 s circles, ten circles/batch, score
last three, budget ten completed candidates, two identical-gain operational
passes to finish. Fixed hysteresis scale/rate/upper = 5/10/50; lambda/Ks remain
[18;16.5]/[15;10], verified in saved model. Initial global search bounds half
to twice baseline; local multiplicative change <=1.5. These are engineering
limits, not proven stable bounds. Preserve adaptive estimates/history, gate
mechanical updates during straight holds, independently retain common-pressure
allocation. Worker runs in a separate MATLAB process over localhost UDP.
Runtime controls C2 reference ramps, actual batch timing, independent fixed-raw
scores, candidate acknowledgement, settings/gamma checks and worker heartbeat.
Terminal sequence withdraws reference, disables PAC pressure with its existing
slew, then stops after command pressures <.01 bar for one second. New copy's
dashboard pressure switch is rebound to the actual enable constant and defaults
Off; saved source dashboard referred to missing Kill. Original models preserved.
No arm model compiled/run. Synthetic library/runtime/transport and device-free
clone wiring tests performed; real hardware integration remains unverified.
See ONLINE_SIMULINK_WORKFLOW.md for launch, wiring, reasons and limitations.
Normal desktop simulation only; watchdog cannot protect against Simulink/PC
stall. No commit/push, hardware execution or automatic old-campaign reopening.

2026-10-28: user requested unattended between-batch tuning. Added an offline
host supervisor (pacTuneOnlineStep/Config) with hold/apply acknowledgement,
settle/run/score sequencing, local bounded BO, repeated operational targets,
budget and terminal failure handling. See autotuning/ONLINE_WORKFLOW.md.
Not arm-ready: model/transport/watchdog/evaluator adapters and hardware checks
remain; no SLX edited or hardware run. Gamma values and new scale-5 search
bounds are unconfigured. Preserved-state BO remains a history-dependent heuristic.

2026-10-28 proposal: continuous operation with constrained Bayesian gain tuning
between circle batches, preserving adaptive estimates for operational improvement,
then potentially scheduling validated gains by speed. User reports stable but
inaccurate 5 s circles; this new observation has not been independently evaluated.
The detailed extension in PAC_AUTOTUNING_FRAMEWORK_PLAN.md records the sequence,
history-dependent comparison problem, objectives and staged integration.
Proposal only: no controller/SLX changes, hardware execution or new campaign budget.

Latest2026-10-28: smooth_candidate10 evaluated; current campaign_smooth_v2
budget exhausted(10/10),no candidate11 generated. Original filter verified by
exact replay,settings match,no tracking faults/limits. Late l2/l3 RMSE[2.88020,
3.16638]mm fail2mm; roughness.0160776m/s^2 is smoothest measured but infeasible.
Best qualified smoothness reference remains original gamma.5/.002/5: length
~1.55/1.87mm,roughness~.0193 with old filters. Historical screenings are not fresh
confirmations. Review objective/constraint-margin modeling and budget before
more experiments. See campaign_smooth_v2/completed_search_review/README.md and
tradeoff plot/CSV. Old campaign history and raw data remain unchanged.

Latest2026-10-28: user approved broader-range performance search. Prepared NEW
campaign_smooth_v2 in unloaded T10 experiment folder: old measured.20/.30s filter,
smoothness-primary fixed-raw scoring,2mm ceilings, gamma k[.1,1],d[.0005,.008],
h[1,10]. Exact replay verified all8 reused old-filter recordings; changed-filter
trial excluded. Baseline is now best screened by roughness, not candidate08.
New immutable optimizer.continueSearchAfterFeasible=true permits search beyond
feasible historical baseline; old confirmation-first behavior unchanged. Eight
seeded observations leave2new evaluations under10total budget; confirmations
need review/extension afterward. Queued smooth_candidate09 gamma k=.90453500157774203,
d=.0022428161211506688,h=2.8192356276278843. New RUN_NEXT.md has exact instructions;
restore original filter before running. No live gains or hardware/model execution.

Smooth candidate09 now imported: original filter verified by exact replay0,
matching tuning/no edits/faults/pressure limits. Late lengths[2.52791,2.75200]mm
fail2mm; velocity.00351522m/s, fixed roughness.02207634m/s^2 pass caps but do not
improve qualified baseline smoothness. valid_infeasible. Nine observations/one
slot remain. Queued final smooth_candidate10 k=.25306817650802854,
d=.00069737816290234929,h=2.6573766973344188. See RUN_SMOOTH_CANDIDATE10.md.
Predicted meanroughness~.02030/std.00378 does not by itself predict improvement;
uncertainty influences EI. Stop/review budget after this trial; no extra runs
or confirmation authorized automatically. Original baseline still best screened.

Latest filter outcome2026-10-28: user ran filter_smoke_...12cycles at actual
candidate08 gamma(.75,.008,7.5). New filter verified by exact velocity replay.
Matched-gain fixed-analysis roughness worsens .0239697->.0340889m/s^2; qRMSE
[1.1725,1.3029]->[1.1894,1.4165]mm. User prefers restoring old filters; recommend
original measured.20/.30s,desired.10/.15s. Prepared paste-ready
pacMeasuredDerivativesOriginal.m; no SLX edits by assistant. Smoothness-first
draft now specifies old filter, still2mm coordinate ceilings and fixed raw-length
scoring. Changed-filter run is not mixed into old optimizer ledger. See
filter_review_v1/actual_arm_comparison/README.md. Historical prototype notes below
are superseded by this actual-arm outcome; do not install the experimental filter.

Latest direction2026-10-28: user prioritizes smoothness and approved2mm RMSE each
coordinate for a NEW campaign. Current1mm campaign remains immutable; no next
hardware candidate queued. Measured-filter replacement pacSmoothMeasuredDerivatives.m
is paste-ready/offline-tested only (.15s two-pole velocity,.35s acceleration).
Saved model measured filter .20/.30s replayed exactly; desired filter .10/.15s.
Measured acceleration is not used by PAC force. Proposed replacement reduces
non-fundamental content but adds ~.07-.08s velocity-estimate lag at10s period.
See experiment folder filter_review_v1/README.md. Next manual action: measured
filter-only3circle smoke with baseline gamma .5/.002/5; no SLX changes by assistant.
New objective/source policy supports roughness minimization from raw lengths
through a fixed offline filter, preventing changed live filters gaming the score.
pacTuneSmoothCampaignDraft proposes k[.1,1],d[.0005,.008],h[1,10], not an active
or physically validated search. New physical filter baseline required first.

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

First physical workflow protocol prepared: user saved PAC_autotuning_manual.slx
and selected no payload, phi=1 rad, T=10 s. Proposed sequence is 3-cycle smoke,
then two fresh 12-cycle baseline repeats after log review. Save originals under
passivity_adaptive_control/experimentations/autotuning/no_payload_phi1_T10_20261027;
see its README for names, timing assumptions and notes. Actual stable tuning
and reset/trim details still need verification from the first log. No trial has
been executed by the assistant; these are not yet formal confirmation runs.

First user-acquired smoke + two baseline recordings now reviewed offline; see
that folder's analysis_review_v1/README.md. Verified .2 bar activation, slew6,
lambda18/16.5, Ks15/10, gamma k.5/d.002/h5; all fixed across runs. Initial
estimates/history match. Baselines: no latched faults or pressure limiting,
99.17/98.83% valid coverage; last-three-lap thetaMAE .15202/.14836, phiMAE
.06248/.06234, lengthRMSE [1.568,1.889]/[1.529,1.857] mm. Learning continues
through lap12, with residual roughness; no steady-state convergence claim.
No final-three-lap parameter bound contacts. Raw-log/evaluator smoke check passes;
formal optimizer success remains unset pending targets/bounds/operator metadata.

Length-based acceptance implemented at user request: new evaluator schema
pac-tune-2 requires per-coordinate q2/q3 RMSE, vector velocity-error RMS,
velocity-error-rate roughness RMS and pressure-limiting limits; vector q RMSE
remains the objective. Angular/vector-RMSE acceptance limits are optional (NaN
disables); optional angular observability does not veto valid length acceptance.
Old pac-tune-1 snapshots retain their angular policy. Numerical limits remain
unset, not inferred from baseline. New baseline diagnostics are saved under
analysis_length_policy_v2: late velocity error ~.00300/.00297 m/s, error-rate
roughness ~.01934/.01918 m/s^2. These are observations, not target recommendations.

User subsequently agreed l2/l3 RMSE <=1 mm each, vector velocity-error RMS <=.01
m/s and velocity-error-rate roughness <=.025 m/s^2. Campaign-specific
agreedEvaluationConfig.m in the unloaded T10 experiment folder records these;
generic defaults stay unset. Pressure-limiting target, gamma search bounds,
budget and reset/trim/preconditioning details remain to be finalized before
optimizer-selected physical trials. No configured campaign or trial was launched.

User then approved gamma ranges k[.25,.75], d[.001,.008], h[2.5,7.5], 12 circles
scoring last3, zero limiting and two fresh confirmations, and attested full reset.
Prepared experiment-folder campaign_length_v1 from agreedEvaluationConfig and
both existing baselines; both are valid_infeasible on length RMSE. Explicit budget
is conservatively10 total observations including2 baselines/confirmations (8 left).
Archived model/current evaluator source snapshots; model history and physical
reset/calibration provenance are labeled as snapshot/user attestation, not measured.
Queued candidate03 gamma k=.645028694505006,d=.0010090123623767366,
h=5.8362570107825071. Full tuning and RUN_NEXT.md saved; no live gains applied.
150s schedule:startup20,tracking120,ramp-down5,pressurized hold5,then manual kill.
Do not batch without explicit controller/history/supervisor and physical reset work.

Candidate03 user-run log now imported as candidate03_r01: exact tuning, no live
edits/faults/limits, late qRMSE[1.42549,1.69001]mm/vector2.21092mm (~9% below mean
baseline), velocity.00314726m/s, roughness.0207654m/s^2. Both1mm length limits
still fail; other configured limits pass. Single-run comparison, not confirmation.
Queued candidate04 initial-design gamma k=.58034718095284454,
d=.0037139224909066579,h=4.4053683649622428; RUN_CANDIDATE04.md has instructions.
Three observations counted,7 remaining including confirmations. No gains applied
or hardware execution by assistant. Late roughness rising but below agreed bound.

Candidate04 imported as valid_infeasible: exact tuning, no edits or tracking
faults/pressure limits. Late qRMSE[1.69217,1.94034]mm/vector2.57456mm (16.4%
above candidate03), velocity.00302789m/s, roughness.02072852m/s^2. Only length
limits fail. Candidate03 remains best observed, not confirmed. Queued candidate05
initial-design k=.28570361215890433,d=.0077863906974842481,h=3.5952162270591348.
RUN_CANDIDATE05.md saved; four observations counted/six remaining. No gains applied.

2026-10-28 candidate05 imported valid_infeasible: exact tuning/no edits,
tracking faults/pressure limits0. Late qRMSE[1.83602,2.06862]mm/vector2.76589mm,
velocity.00372471m/s, roughness.03358536m/s^2 >.025. Reject on smoothness and
both length limits; roughness elevated fromlap7 throughlap12(.037364). No ratecap
contacts. Candidate03 still best observed. Initial-design candidate06 queued:
k=.42192531104075343,d=.0020323562614174683,h=2.5725609289380671.
Five observations/five remaining inclconfirmations. See RUN_CANDIDATE06.md.

Candidate06 imported valid_infeasible on2026-10-28: exact tuning/no edits,
no tracking faults/pressure limits/ratecap contacts. Late qRMSE[2.63375,2.90327]mm,
vector3.91990mm; velocity.00344364m/s and roughness.01935062m/s^2 pass their
limits. Smoother than05 but less accurate; candidate03 still best observed.
Queued candidate07 final initial-design k=.35438912779204967,
d=.0030892129738313034,h=6.8259366991000565. Six observations/four remaining
including confirmation. RUN_CANDIDATE07.md saved. No live gains applied.

Candidate07 imported valid_infeasible: exact tuning/no edits, no tracking
faults/limits/rate caps. Late qRMSE[1.23668,1.66263]mm/vector2.07213mm (best
observed,6.3% below03), velocity.00299571m/s, roughness.02258509m/s^2. Both
length limits still fail. Seven observations/6 unique locations complete initial
design. Candidate08 queued seek_feasible_region at upper box edge:
k=.75,d=.0079999999999999967,h=7.4999999999999991. GP mean2.56631/std.60053mm
does not predict improvement; weighted binary feasibility is not a physical
probability. Explicit boundary exploration, with05 roughness caveat recorded in
RUN_CANDIDATE08.md. Three observations remain including confirmations. No gains applied.

Candidate08 imported valid_infeasible: exact tuning/no changes, tracking faults/
pressure limits0. Late qRMSE[1.17252,1.30289]mm/vector1.75280mm (best observed,
15.4% below07), velocity.00299378m/s, roughness.02391251m/s^2. Both length limits
fail. Pooled smoothness passes butlap12 .026547 exceeds .025 individually; retain
caveat without changing agreed scoring. Eight observations/two slots left, no
accepted candidate; additional search plus2confirmations cannot fit remaining
budget. Raw next seek_feasible_region suggestion saved, NOT queued: predicted
2.69475mm/std.63161mm favors uncertainty, not closeness to the1mm constraints.
Review component-margin acquisition before more trials; current binary-label
algorithm treats all observed infeasibility alike. No algorithm/config change
made. REVIEW_AFTER_CANDIDATE08.md records findings and next-step recommendation.
No causal attribution to one gamma; multiple changed, plus different-day trial.

2026-10-28 research discussion: PAC_CONVERGENCE_RATE_RESEARCH_NOTE.md compares
Milad's paper Eq11/19/TableI to code. Paper Gamma is inverted (1e-5/1e-2 ->
update multipliers1e5/100), while our G=diag(scale^2*gamma). Nominal scalar poles
are -lambda and -(actual d+Ks)/m under exact matching/force delivery; these do
not assign adaptation convergence in the actual gated/sampled arm. Conditional
practical error/time-envelope derivation and prior art documented; no new theorem
proved, novelty established, controller change or campaign mutation.

User perceived05 quicker: follow-up shows first-lap l2 error lower than03
(3.258 vs3.967mm) but l3 higher(6.219 vs5.668mm). Approximate10s harmonic
lag late03 [.185,.226]s vs05 [.292,.321]s; no overall faster phase response
established. Current acceptance scores sustained last3circles, not first crossing.

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
