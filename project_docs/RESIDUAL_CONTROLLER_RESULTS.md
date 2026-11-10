## 2026-11-09 - Sine test1 pressure-vector wiring fault proven

Newhelperexecuted4241ticks;source replaywithactualA8andMC3nominal matches
pressure1.49e-14bar andallmode/age/flags. OutputsfiniteALLticks. OldSwitch
stilltreatshandoffpressureascontrol:>0selectrawcandidate,else rawsine.
DelayofthismapmatchesactualsavedselectedpressureEXACTLYincludingNaNs.
5realNDImissestrigger435NaNcommandrows/435mode6faults,434observerinvalid
withvalidsensors. Firstsensorloss145.90->NaNcommand145.95. Analoghandling
unknown;cannotclaimphysicalchamberpressureNaN. Thisisnotproperhandofftest.

Manualoutput1->UnitDelaydirect/bypassSwitch required. Sourceinput3also
MC3nominalratherthanoriginalrawsine(.091barsecondarydifference). KeepKill.
No runtimechange. Inverse868/3267scaled3bar;othermodelissuesstillpossible.
Offline analyzeSinePressureRun andcanonical Experiments/Model_Val/
sine_pressure_review_20261109_v1 retainexactreplays/plots/data/hashes.
NoSLX/model/hardware/trainingrun,commit,push. SaveddiskSLXstale;preserve
unsavedoperatorworkanddon'tinferexecutednewhelperfromdiskalone.

## 2026-11-09 - Cygnus confirms independent tip tracking loss

RawDesktoplog10369BXallCRCsvalid;after8startupfinalreports,23tipMISSING
reports/17episodesand0basemissing in10359runtimeobservations. Device-level
lossindependentofMATLABangularprocessing. Userconfirmsstandardsine/noMAT.
Motion-alignedP3windowscontainall23misses(12/11perrepeat);labelingisinferred
fromsequenceandmotion,notloggedpressuretimestamps. Misses~25-50ms updates;
nogrossvolume/systemstatusflags;uniquephysicalcauseunknown. Handoffsmoothing
doesnotrepairtrackerloss. Canonicalraw/CRC/CSV/events/plots/hashreport
Experiments/Model_Val/cygnus_sine_20261109_v1. Noruntime/hardware/modelchange.

## 2026-11-09 - Stateful pressure handoff source implemented

New residualCommandValid.m changes Booleanhelper to7input/4output pressure
handoff. BypassoldSwitch,retainUnitDelay.05 anddownstreamKill. Inputsactual
commandfeedback,clock,enable/reset pluscandidate/ready/originalfallback.
.10shold,.50sfades,shared6bar/s0..3barlimit;repeat-tickcache andclock/input
faultguards. NoGRUgateorhistorychanges. No fullobserverwarmuphold.

Syntheticand7501frozenrecordedinputcommand-pathchecks pass,sharedordinarystep
<=.30bar,0steps>.3 vsold30,max1.514. Counterfactualcommandsnotpredictionof
physicalphi/newNNstate;actualbuild/timing/hardwarebehaviorunverified.
Manual COMMAND_HANDOFF.md;artifacts command_handoff_validation_20261109_v1.
Current sinewaveoperatorrun leftunchanged. NoSLX/model/hardware/trainingrun,
commit,push.

## 2026-11-09 - First residual desired-pressure actuation test5 audit

Operator model1.59 completed375s.49orientationinvalid rows alltipMISSING;
9hostsamplegaps>.1s(max.636s),basealwaysvalid,noNNguard,statusage0.
37ready-losses(32missing,5gap) ->fallback.05s later,verifiedexactSwitch/
delayreplayandactualcommandmatch4.45e-16bar;Killalways1. Medianfallback
1.65s,max2.60s fromvelocity/historyrewarm.32channel-drops>.1bar,max1.514bar;
26stable-referenceepisodesfollowedbyphidips>.1rad. Sourceinvalidprecedesthe
fall;pressurecollapseisconsequenceandfurthertrackingdisturbance,notinitial
tipstatuscause. PhysicalNDIfailurecauseuniquenotknown;repeatsatcirclephases.

Phi1lengthRMSEnominaltest4[8.203,8.173]mm ->test5[3.875,3.923]mm;desired
referenceidentical,separaterunsnotcontrolledrepeats. Phi1readyMAE.157rad,
fallbackrewarm.464;phi2.5still[11.988,11.285]mm. Fullcorrectedpressure
scaling1012/6143,slew122. RawfallbackSwitchbypassesbranchslew;needsstateful
bumplesshandoff/sharedfinalratecap,notgammaorvalidityrelaxation. Proposalonly.
New analyzeDesiredPressureRun.m andcanonical desired_pressure_test5_review_
20261109_v2 underExperiments/Model_Val retainplots/events/CSV/MAT/hashesevidence.
No runtime/model/hardware/trainingchange/run,commit,pushduringaudit.

## 2026-11-09 - Saved V2 command-selection wiring audit

Read-onlySLX confirmscorrectSwitchtruecorrected/falseoriginalnominal,
validityfunctioncontrol,UnitDelayTs.05ICzeros3x1 beforevalves,coherent
postKillpressure/timefeedback. Function9validinputiswrongsource:A5baseNDI
orientationvalid. Required:BusSelectorvalidfromFunction8desiredPressureDiagnostic
output3 ->Function9input2. FinitecheckalreadyrejectsNaNs, butreadinessmustuse
diagnosticflag. No furthercodechange indicated. Audit/hash:
residual_controller/wiring_review_20261109_v1. Noload/compile/edit/run/commit.

## 2026-11-09 - Operator direct-actuation algebraic loop

V2diagnostics/savedwiring show inverse->valves->pressure/time->persistent
observer->inverseContext->inverse, involvingobserverports6/7 andinverseport1.
Manual ACTUATION_CAUSALITY.md describesvalid+finite selector with existing
protectednominal/startupfallback, thenUnitDelayTs.05/ICzeros3x1 beforevalves.
PhysicalKill/enable/stopremainpriority. Keepactualcommand/timestampfeedback
paired;commanddelay adds50ms andretainsforecastinputmeaning betterthansimply
delayingfeedbackalone. This iscausalplumbing,notcompleteclosed-loopqualification.
No source/modelchange,compile,hardware/modelrun,commit,push.

## 2026-11-09 - Operator forward-wrapper dimension failure

Connecting desiredpressure toactuators exposes compilation errorat
armS_single_dynamics [dq;acceleration], unresolved dX in MATLABFunction1.
Savedforward port sizes inherited-1; screenshotdoesnotidentifyfirstwrongshape.
New paste-ready RL_data_collection/rlForwardModelBlock.m explicitlynormalizes
pressure3x1/state4x1,torque2x1,RHS4x1 withsame3inputs7outputs andoriginal
constants/EoM/clipping. Shared dynamicsunchanged. Offline50cases/100row-column
variants matcholdRHS8.882e-15. UserSymbolsizes/manualguide FORWARD_BLOCK_SHAPES.
Report residual_controller/forward_shape_validation_20261109_v1.

Directactuation remainsseparate:desiredblockreturnsNaNduringcold/invalidctx.
Existingstartup/hold/fault/Killsupervision mustbe retained; simplewirechange
doesnotimplementsupervisedpressurehandoff. Noactualmodelcompile/hardwarerun,
SLXedit,commit,push.

## 2026-11-09 - Desired-motion subtraction source implemented

rlResidualDesiredPressureBlock.m now computes desiredq/v/a nominal demand,
then subtracts observedresidual. It reuses5ports/3outputs/diagnosticbus but
input2 is6x1 desiredMotion, requiring manualMux6 referencebranch. ActualGRU
features/history unchanged, no desiredstate historysubstitution or extra
recurrentcommits. This is frozenresidual desiredfeedforward, notexactinverse.
Oldsources and PAC/PDFB preserved;guide DESIRED_PRESSURE_BLOCK.md.

64/128synthetics and5433savedtest4rows pass (5118valid/315NaNinvalid).
Desirednominal=requested4.45e-16bar;desiredforce3.56e-15N;residualunchanged
1.07e-14N;independentmapping/clipping<1.34e-15bar. First70s correctedpeak
[.523,.468,.484]vsnominal[.242,.240,.242]bar,no slew. Fullrecordcorrected
scaled126/slew3027, modelbalanceerror med7.107N frozen/7.615N candidate;
limits andpressuredependence materiallyaffectlargerdemands. Notphysicalforce.

CoderC emitted,localstandalonecompiler existingIndexerror;independentGCC
764nativefixtures pass2.203e-13. Reports/plots/CSV/MAT/sourcehashes:
residual_controller/desired_pressure_validation_20261109_v1 and
desired_pressure_native_20261109_v1. No SLX/model/hardware/trainingrun,
commit,push. Operatorpaste/wire/compile required, outputs initiallyscopeonly.

## 2026-11-09 - Single-pass shadow is not a desired-trajectory inverse

User's test4 pressure comparison exposed a purpose/target mismatch. Collection
uses desiredq/v/a; current collection-policy shadow uses actualq/v+desireda.
It has no desiredlength/velocity input. Matching pressure allocation did not
make the targets equal. Recorded5118valid contexts confirm exacttargets and
source replay1.12e-15bar; collectiondesirednominal=requested to4.45e-16bar,
actualNI=requested exactly. No corrected slew/scaling in this recording.

Firstcircles desiredpeak2.600mm versus actual[.374,.180]mm. Corrected
shadow nearly reconstructs appliedp~.24bar; it does not compute pressure
to reach the larger desired motion. Offline change to desired force demand
with SAMEactualresidual gives peaks[.523,.468,.484]bar, but pressure-dependent
modelbalanceerror acrossrecord median1.478N confirms frozen correction is
approximate. This is not a physicalforce measurement or validated inverse.
Forwardmatched50ms residualRMSE[.226,.231]mm versus CV[.055,.058]mm;
good measured-state nextstep prediction alone cannot certify inversecontrol.

Next interface should accept shared nominalForceDemand_N2x1: existing
collection tauDesired for desired inverse comparison, complete PAC force
demand for PAC. Keep GRU actualstate/actualpressure/immutablecausalhistory;
do not feed desiredstates into that stream. Samepolicy/constraints and pre/
postlimit diagnostics needed. No runtime implementation changed thisaudit.
Report/figures/CSV/MAT/hashes: residual_controller/Experiments/Model_Val/
target_mismatch_review_20261109_v1. Offline analyzeSinglePassTargetMismatch.m.
Guide andPROJECT_STATE updated; no SLX/hardware/training/commit/push.

## 2026-11-03 - Residual PDFB/PAC source implementation and offline verification

Implemented residual_controller PDFB/PAC wrappers, frozen64/128 numeric packet
loader, bounded learned-pressure inverse, one-commit GRU state, typed diagnostic
bus and manual installation instructions. Original feedback/adaptive law/gains/
filters/references unchanged. No SLX was edited/opened/run. Saved PDFB base-rotation
inverse still had old bounds; new off/on pair uses corrected bounds and shared
supervisor/allocator, so a fresh new-wrapper nominal baseline is required.
PAC first28 model-off outputs exactly match existing wrapper; theta-theta0 now
logs remaining correction when model-on. Fixed-nominal activation/gravity
re-expression and trim alignment derived; no new passivity guarantee.

Canonical tests: residual_controller/validation_20261103_v4. Exact frozen-cell
agreement, FD Jacobian9.8e-11,48 inverse cases max8.82e-6N,804 full-wrapper
synthetic-kinematic ticks pass bounds/slew/hold/fault checks. MATLAB Coder emitted
both wrappers but hit existing toolchain error; independent GCC build and2576
native ticks passed (numeric<6e-12, pressure<3e-14bar). Strict activation-boundary
logging flag differences documented. This is component evidence, not physical
tracking or successful user-model build. Setup/diagnostic bus test passed.

Read residual_controller/README.md, MODEL_INTEGRATION.md and VERIFICATION.md.
Next: operator installation/build and model-off smoke, then fresh matched
model-on trial. Confirm actual live PDFB settings and timestamps first. No
hardware, retraining, HIL automation, new Git commit or push performed.


# Residual controller results

## 2026-11-06 - Operator test2 circle/helix confirm collection-policy block

Found circlesweep_FWD_INV_test_2.mat (user.at typo) and
helix_FWD_INV_singlePass_test_2.mat (actualhelixfilename). Model1.53completed
375s7501rows/290s5801rows, no saved errors; actualsaved options present.
6961/5257valid inversecontexts. Fullcollection-source replayoutputs1.34e-15bar,
diagnostics1.07e-14,flags0mismatch; independentmapping/clipping4.45e-16bar.
Actual/contextNI andforward/inverseresidual exact, networkevaluations1.
Newminimum-policy operator model execution now evidenced, no codefix indicated.

0.3startupartifact gone:0all-three-.3 vectors during979/852zero-commandrows.
First1.65s correctedcircle[.0200,.0168,0],helix[.0254,.0093,0]bar. Nominal/
corrected differallpre-slewrows, only3/0equalpost; maxprecorrection1.628/1.621bar.
Correctedslew774/6961(11.12%)/859/5257(16.34%);scaling64(.92%)/40(.76%).
Candidate-pressure modelerror median.06110/.10035N,P954.0416/3.2014N,
max8.4019/8.1804N. >1Ntail733/735 and750/751coincidewithlimits.
Frozenbalance is not nonlinear consistency; only1154/1321meet.02N after
reevaluation. Olditerative stillPACpolicy, notmatchedsolvercomparison.

NDItipMISSING16/18,baseinvalid0;helix1119valid-nearstraightthetaunobservable
allaccepted(phi<=.008654rad), notdropouts. No forwardnumericalguards.
Forwardall-scored lengthRMSEmm nominal/residual/CV:
circle[4.407,4.824]/[.848,1.283]/[.100,.114];
helix[4.526,4.650]/[1.085,1.281]/[.124,.137]. Prior-issued/error/rawqidentities
exact. ActualNI equalsoriginalcollectionrequestedpressure onvalidrows;
these remain shadowobservations, notresidual-onPAC/PDFB tracking evidence.

Canonicalreport residual_controller/Experiments/Model_Val/
collection_runs_review_20261106_v2: savedsettings,JSON/CSV/MAT,pressurefull/
startupPNG/PDF/FIG,source/rawhashes. Newoffline analyzeCollectionShadowRuns.
Initialv1ambiguousoverrun source markedincomplete; v2explicitNDIsource.
Next PACmath must retainpressure-limit/frozenresidual terms andgating/dropout
recovery. No hardware/modelrun,SLX/runtime/controller/NN edits,training,
commit orpush duringaudit.

## 2026-11-06 - Collection-policy single-pass replacement source ready

User authorized updated code. New paste-ready
rlResidualSinglePassCollectionShadowBlock.m replaces ONLY pressure conversion
with collection minimum-nonnegative rule plus maximum-pressure scaling.
Same5inputs/3outputs/ResidualSinglePassDiagnosticsBus and.05s schedule.
Single GRU atactualpressure, measuredq/v+desiredddq target, readiness/NaN and
immutablehistory unchanged. Rails/slew retained; compatibility variableCommon
and common/trim/PACs-a-cap slots ignored underthis different policy.
OldPAC-policy block/iterative/PAC/PDFB sources preserved. No SLX edits.

Offline canonical single_pass_collection_20261106_v1: bothwidth synthetic
feasibleforce/sign/zero-preload/scaling/rails/slew/invalid/port-invariance pass.
All5266valid native-inputhelixcontexts match independentminimum-policy
counterfactual1.34e-15bar; targetexact,force1.43e-14N. First1.60s corrected
[.022919,.011821,0],nominal[.009561,.038671,0]bar; zero-command all-three-.3
rows0/853. Correctedslew861/5266 vsold5099;nominal3336,scaled30corrected.
Candidate-pressure model mismatchmedian.1010N,not measuredforce/tracking.

Coder emitsC but existingstandalonecompiler failsIndexexceedsarraybounds;
independentGCC/native764rows passmax1.99e-13. GenericofflineCfixture/check
tools now accept optionalentryname,originaldefaultunchanged. Source/data/hash,
CSV/MAT/startupPNG/PDF/FIG andnativefixtures retained. Manualinstallation:
residual_controller/SINGLE_PASS_COLLECTION_SHADOW.md. Userpastesentirefile
includinglocalhelper intoexistingobservationblock, no rewiring/newbus required.
OlditerativeblockstillPACpolicy, so notmatchedsolvercomparison withnewvariant.
No hardware/modelrun,training,actualpressurecommands/controller changes,
commit orpush. ActualoperatorSLXcompile/full-loop/physicalinverse pending.

## 2026-11-06 - Operator single-pass helix:0.3bar is allocation/slew artifact

Reviewed helix_FWD_INV_singlePass_test.mat, model1.52,290s5801rows, no saved
simulationerror;5266jointvalidcontexts. Atfirstvalid1.60s actualNI=[0,0,0];
PAC nominalpre=[1.59348,1.62259,1.58392],singlepasspre=[1.61134,1.60024,1.58842]
bar; all3inverseoutputs=[.3,.3,.3]. ALL853valid zero-command rows have that
output. ReusedPAC nominalcommon1.6bar plus6bar/s*.05=.30bar peractualcommand
clipping causes this; not NN-learned deadzone or actualpressure starting.3.
Stateless shadow cannot accumulate a hypothetical ramp from its old output.

Singlepass=iterative3960/5266(75.20%),singlepassslewlimited5099(96.83%),
nominalslewlimitedall5266. Pre-slew singlepass differsfromnominalallvalidrows,
max.97574bar. Singlepassvsiterative maxaxis difference median0/P95.28666/max.6
bar. Candidate-pressure equationerror medianssinglepass.60208N/iterative.49672N;
nominal8.03121N. These are model algebra errors, not measured forces.

Trainingusedrecordedpressureforwardinputs,noPACcommon/slewpolicy. Saved
collectioninverse ZIP/XML usesminimum-nonnegative pressure then3bar scaling,
not PAC1.6common. Suppliedshadowwrappersintroduceddifferent allocation, which
obscures the requested comparison. Offline SAMEforce-target policy comparison
firstnominalmin=[.009561,.038671,0],singlepassmin=[.022919,.011821,0]bar.
5186force-feasible mappings preserveBp<=7.11e-15N, but pressure-dependentr
means this does not prove equal augmented/physical response. Target remains
measuredq/v+desiredaccel, distinct from collectiondesiredstate demand.

Alltargets/residuals/nominal/clipidentities exact.120contextsource-default
variablecommonreplays singlepass8.89e-16bar,iterative7.27e-9bar; fulloptions
not saved. NI uses nativefastergrid; timestampboundedalignment confirmscommand
age0 andexactcontextmatch, not a runtime timestampfailure. Initialanalyzerv1
strictgridassert markedincomplete; canonicalv2 handlesit.

Artifacts: residual_controller/Experiments/Model_Val/
single_pass_helix_review_20261106_v2 (README,CSV/MAT/JSON,full/startupPNG/PDF/FIG,
savedchart/source/inputhashes). Newoffline analyzeSinglePassHelixLog and
compareCollectionShadowAllocation only. Next proposedpassive comparison:
matchcollectionallocation and show pre-slew correction separately fromlimited
output. No runtimeblock/controller/settings/SLX/NN edits, hardware/modelrun,
training,commit orpush. Original inverse is iterative search, not integration.

## 2026-11-06 - Requested real-time motion/angle/length/error figure set

User clarification selected later runtime MATs circlesweep_test,helix_FWD_INV_test,
sinewave_test andsinewave_test_100g. 12figurestems in
residual_controller/Experiments/Model_Val/presentation_shadow_20261106_v2:
perrecord angles+errors andlengths+errors,plus completecircle159-169s and
helixfirstepisode5-65s. Each4320x2430PNG/PDF/SVG/FIG. Fullcurves/errorCSV,
MAT,metrics,metadatawindows,signalpaths,source/inputhashes and adjacentZIP.
All4rawq/loggederror/priorissuechecks exact; anglegeometrymatchesNDI<=9e-16.
Commonangularrows for comparative thetaMAE. LengthRMSEmm nominal/augmented:
circle[4.296,4.714]/[.748,1.220],helix[4.388,4.534]/[1.022,1.180],
sine0g[3.696,2.897]/[.482,.429],100g[4.505,3.733]/[.729,.809].
These are measured-state-initialized50ms predictions, not free trajectory or
controllertracking. No CV shown perrequest; previous CVbenchmark caveats remain.
Theta circular/wrapped/undefined nearstraight; missing/warmupNaNgaps preserved.
Payloadlabeluser/filename,no causalpayloadclaim. No fitting,smoothing,SLX/raw/
controller change,hardware/modelrun,commit orpush. v1 preliminary retained.

## 2026-11-05 - Single-pass subtraction shadow block implemented

User authorized a separate paste-ready inverse observation block. New
rlResidualSinglePassShadowBlock.m evaluates GRU ONCE at actual current pressure,
measuredq/v and immutable historyPre, subtracts that force from the existing
measuredq/v+desiredacceleration nominal demand, then uses bounded allocation
and the same actualcommand-reanchored slew box. No inverse iteration, ODE or
state commit. Existing forward/iterative inverse/PAC/PDFB/NDI/SLX unchanged.
New setupResidualSinglePassShadow returns singlePassOptions=[128;.02] and
ResidualSinglePassDiagnosticsBus. Manual ports/types/Selector/save instructions:
residual_controller/SINGLE_PASS_SHADOW.md. Pre/post-slew pressure diagnostics
separate correction from clipping; frozenBalanceSatisfied is NOT nonlinear
convergence or measuredforce accuracy. All observer readiness gates retained.

Canonical offline checks: single_pass_validation_20261105_v2. Bothwidths,
feasibletarget subtraction sign, immutable/repeated contexts, invalid/NaN,
bounds/slew/error decomposition pass.11964valid savedcontexts: force target
and nominalpressure exact, neuralforce1.07e-14N. Corrected/nominal pressures
differ EVERY valid pre-slew sample (max1.1135/.9767bar circle/helix), but equal
78.75/69.05% after clipping. Frozenforce-error medians.4197/.4427N; reevaluating
GRU at predictedpressure OFFLINE gives.4083/.5815N versus iterative.3772/.4841N.
This quantifies approximation, not superiority/physicalpressure validation.

Coder emits C but existing standalone toolchain fails Index exceeds arraybounds;
not an actual SLX build. Independent GCC compile and764native rows pass,
maxdifference1.92e-13 inclall diagnosticfields. Numeric P95.144/.132ms excludes
acquisition/scheduling/observer, not worstcase full-loop qualification.
Sources/inputhashes/CSV/MAT and PNG/PDF/FIG retained in v2 validation folder;
native fixtures/results in single_pass_native_20261105_v1. No hardware/model
execution, training, controller swap, protection change, commit or push.
Next operator action: manually add discrete.05s block, branch same context/ddq/
pressureoptions/common Boolean, new options2x1; new bus scalar, outputsScope
only. Inspect pre/post-slew pressure and frozen residual flags together.

## 2026-11-05 - Circle/helix inverse and NDI missing-transform audit

Operator circlesweep_FWD_INV_test.mat375s7501rows and helix_FWD_INV_test.mat
290s5801rows, model1.51, ReachedStopTime/no saved error. All17/28 orientation
invalid samples are raw tipstatus2=MISSING/positionnonfinite; base allvalid.
No valid-pair quaternion-only rejection, outvolume/partial/repeat/pairreject,
or eventcoincident overrun. Configuration invalid matches tip loss atNDIgrid.
Physical cause remains unproven. Public meas_position and correctedtip XYZ
hold preceding values at all17/28 losses; internal rawrelativeposition NaN.
Smooth XYZ does not establish better position sensing. Nearstraight valid
theta-unobservable9/796samples below0.5deg are distinct; observer accepts them
and every orientationvalid sample. True losses reset0.6+1s warmup.

Inverse equal79.76%/70.13%, but neuralforce median10.42/8.94N, not disabled.
Bothoutputs limited to actualcommand +/-0.30bar EACH tick (6bar/s*.05s),
not preceding hypotheticalcommand. Inverseconverged156/6908(2.26%) and
245/5056(4.85%); valid does not imply converged. Defaultvariablecommon replay
120contexts/file matches2.09e-9bar and convergenceflags exactly; fulloption
vectors notsaved. Diagnosticinverse targets measuredq/v plus desiredaccel,
not desiredlength/outerfeedback. Residual need not increase allpressures.

Confirmed supplied projectedsolver plateau:16selected equal-output contexts
have zero localprojector Jacobian;15better candidates through IDENTICAL
allocator/slew/rails/caps. Helix39.20s inverseerror4.9088->2.2984N within same
limits. Separate copied-option no-slew probes show nonzero correction up
~0.83bar; not proposed live limit removal or physical inverse validation.
Next: improve disconnected inverse/diagnostics offline before controller
deployment; keep genuinevalidity/readiness, display staleXYZ clearly and
operatorinvestigate tipMISSING without claiming a cause. Shared controller
allocator has same numerical limitation. No productionfix applied here.

Artifacts: residual_controller/Experiments/Model_Val/
inverse_ndi_review_20261105_v1/README.md; event/schema/sourcehash tables,
defaultreplays, sameprojector probes, PNG/PDF/FIG and immutablecontext MATs.
Added offline analyzeInverseNdiValidation/probeInverseProjectionPlateaus only.
No SLX/runtime/NDI/network changes, hardware/modelrun, training, commit or push.

## 2026-11-05 - Circle causality/future-target falsification audit

User questioned suspiciously close measured/predicted plots.375s7501rows,
model1.51;6767prospective targets all exactly previous-issued predictions.
IndependentNDI q/NI commands/error recalculation exact. Full source-default
causal replay reproduces prediction1.39e-16m,forces1.07e-14N,allflags exact;
exact options notsaved but defaults reproduce. Synthetic next target+5/-3mm
leaves old forecast unchanged and shifts error-5/+3mm. No future-target leak;
Nov3packedweightSHA unchanged, no retraining. Observed q resets each origin,
so absolute overlay is substantially easier than predicting dynamics freely.

Circle all/moving RMSEmm:res[.748,1.220]/[1.337,2.204];nom[4.296,4.714]/
[6.321,7.328];CV[.097,.111]/[.174,.199];persistence[.206,.256]/[.375,.464].
22invalidorientation,ready6806/7501,24segments; no numericguards. Velocity
maxstandardizedfeatures18.05/17.13;12.84%readyanyabsZ>5 (not support proof).
320conditional no-reset rollouts(88moving) using recorded futurecommands/dt,
no future measuredstate inputs, .05/.1/.25/.5s, zero guards. At.5smoving
res[1.660,2.906],nom[9.102,9.902],CV[1.748,2.042]mm. Residual does not
consistently beat CV; improve/check pressure-conditioned model behavior and
velocity coverage before treating as NMPC-ready. Logged one-step error is not
inverse/control/stability evidence. SavedScope thirdinput is separate nominal
model/length chain; keep distinct from matched one-step nominalcolumn2.

Artifacts, incremental-motion contrast plot, full replay/counterfactual, all
metrics and sourcehashes in Model_Val/circlesweep_test_review_v1. Added only
offline audit/plot functions and generic analyzer title handling. No controller,
SLX, raw or neuralweight changes, hardware/model run, fit, commit or push.

## 2026-11-05 - Sine-wave presentation exports for0g/100g

Read operator sinewave_test.mat and sinewave_test_100g.mat, model1.51,212s,
GRU128/horizon.05. Plot requested rawNDI l2/l3, matched nominal and residual
predictions only. RawNDI/observerq and context/NIcommands match; scored output
is previous issued output, error recalculation matches exactly. Validscore
coverage91.16%/93.94%; gaps preserved. Sine actualphase5; phase4relaxation.
Matching firsttrial cycle2/channel windows selected without considering errors.

All-scored RMSEmm:0g nominal[3.696,2.897],res[.482,.429];100g nominal
[4.505,3.733],res[.729,.809]. Active-sine:0g nominal[4.582,3.586],res[.594,.529];
100g nominal[5.556,4.587],res[.896,.994]. Bothmodelsmatched one-step conditions;
not nominal free-run versus residual one-step. No CV plotted per user request;
earlier CV baseline findings remain applicable. Payload labels user-supplied,
attachment/calibration timing not recorded; no isolated payload/control claim.

New offline plotResidualSinePresentation.m; fivefigures eachPNG300dpi3840x2160,
PDF/SVG/FIG,CSV/MAT/hash inventory in Experiments/Model_Val/
presentation_sinewave_20261105_v2; ZIP adjacent. Visual review fixed footer
overlap before final export. Originalraw/SLX/controller unchanged; no hardware,
model execution, fit, commit or push.

## 2026-11-04 - pressure_combine_test real-time forward observation

Operator residual_FWD_INV_validation1.50,191s ReachedStopTime/no error. Full
forward/context/diagnostic buses logged; no learned-inverse outputs. Independent
raw NDI q/currentNI pressure match exact; matched predictions equal one-tick
earlier issued output; recalculated errors exact.101context probes frozen-cell
force diff2.66e-15N, full three-column flow diff1.99e-17m. GRU128 horizon.05s.
Physical orientation3814/3821valid;base3821valid; ready3580; prospective3560
scores(93.17%),retimed3573. Seven invalid samples re-warmed; no numeric guard.
Two long poll gaps are startup only. Physical clockmedian.050807s; drift3.361s.

Prospective all-row length RMSE mm residual[.234,.198],nominal[1.007,1.569],
CV[.046,.047]. Moving vnorm>.001m/s1043rows: residual[.346,.300],nominal
[1.700,2.491],CV[.074,.074]. About76.8/87.4% residual RMSE reduction versus
nominal, but CV remains superior at50ms. Not a neural-inverse or physical
controller-performance conclusion. FeaturemaxZ4.118, no readyZ>5, forcespeaks
10.406/10.792N. Finalcommands0..3bar, Kill1. shadowOptions vector notsaved.

New offline analyzeResidualShadowLog.m; retained PNG/PDF/FIG including objective
max-speed zoom, tables, MAT, source/model/input hashes in Experiments/Model_Val/
pressure_combine_test_review_v1. Operator build/forward runtime now demonstrated;
worst-case compute time remains unknown. Recommend frozen longer causal rollouts
and disconnected inverse before residual-on PAC/PDFB; old PAC recovery policy
unchanged. No model/control edits, hardware execution, fit, commit or push.

## 2026-11-03 - Passive shadow block source and offline verification

Delivered paste-ready forward observer and separate stateless shadow inverse;
no SLX modification. setupResidualShadow prepares paths/hash checks/typed buses.
Manual ports, pressure/frame source locations and scope/save instructions:
residual_controller/SHADOW_VALIDATION.md. All outputs diagnostic only.

Raw inferredq, internally segmented tau.2velocity, .6s preprocessing then1sGRU
warmup; physical poll/frame separate from simulation command freshness. Gaps,
invalids, backward clocks and repeated frames reset observer segments without
latched actuator faults. Duplicate identical calls neither score nor advance.
0.05s issued predictions with.005s match tolerance separated from elapsed-time
retimed replay. No future measured state/command used for either prediction.
Inverse uses immutable current pre-update history and bounded/slew-limited
hypothetical pressure allocation; not PAC demand or physical force validation.

Canonical numeric shadow_validation_20261103_v4: exact propagator and sequential
cell; zero residual reduces to nominal; causality/wrapper/both widths/dropout
1/2/3 and timing/command/base-invalid regressions pass. Known-force inverse
fixed/variable-common consistency passes with context unchanged. firstTime
observer ready23s after controller fault20.8. Pre-fault prospective380samples:
augmented length RMSE[.694,.622],nominal[1.537,1.747],CV[.111,.098]mm. Mainly
trim/blend, only.65s circle; no circle-tracking qualification. Separate postfault
stationary metrics. All1935rows retained, with PNG/PDF/FIG, CSV/MAT and hashes.
Numeric forward observer P95.798ms excludes acquisition/Simulink/optionalinverse.

Coder emits C then known environment standalone compiler/index error. GCC
independent build/native1326ticks (64/128,widths,bothcommon,warmup/dropout/reset)
passes max5.13e-12 numeric difference, no flag exceptions. This is not the
user SLX build/hardware qualification. Source and artifacts indexed separately
from the previous implementation checkpoint in SHADOW_ARTIFACTS.csv. No model
edit, hardware, fitting, commit or push.

## 2026-11-03 - firstTime PAC history-reset shutdown

Offline saved-data/helper replay proves shutdown20.80s follows invalid tip
orientation20.70/20.75s and0.15s recovered-sample gap. maxGap0.10 resets GRU;
active-ready requirement plus invalidGrace0.10 latches status5. Pressure then
slews zero. Not an error-magnitude trip or a recorded solver error. Numerical
inverse converged immediately before dropout; physical accuracy remains open.
Residual enabled throughout, contract valid, common adjustment enabled; no
NI-command mismatch. Replay status/ready exact; max pressure/raw-force errors
8.53e-14bar/3.98e-13N. Nominal supervisor seeded atlastgood recovers tostatus3.
Prior tests missed two-invalid-then-recover boundary. No controller fix here.
Common preload retained, variable adjustment optional with matched comparisons.
Next: passive forward/inverse validation and explicit dropout/time policy;
physical timestamps differ from logged validity-latched simulation time.
Saved SLX lacks executed residual edits; save actual runtime source snapshot.
Report, plot, CSV, MAT and reproducible audit source:
residual_controller/Experiments/PAC/20261103/firstTime_review_v1/README.md and
residual_controller/analyzeFirstResidualPac.m. No model/hardware run by agent.
