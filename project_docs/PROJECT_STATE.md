Last Updated: 2026-11-10
Current Milestone: Sine test1 executes the new handoff but retains the old Switch, interpreting pressure as Boolean control and injecting NaN commands. Correct the manual command-path wiring before further actuation or evaluating the handoff's physical performance.

# Current project state

## 2026-11-10 - Residual integration and independent IMU Git checkpoint

User requested commit/push of the accumulated project work and the newly added
Yost-Simulink-Direct package. Keep the IMU package in a separate commit for
selective cherry-picking; no destination branches have been specified. Checkpoint
scope, local-only artifact policy, and offline checks are recorded in
project_docs/GIT_CHECKPOINT_20261110.md. Saved SLX files are archived unchanged;
no hardware/model run or live sensor access was performed. The command-path
wiring issue above remains unresolved by publication.


## 2026-11-10 - Literature supports in-volume orientation/quality loss

Userrequestedonlineevidence fornearFGtrackingloss. ClosestmatchMandal2016
thesisprintedp84Fig6.2showsAuroraBadFitcertainorientationwhileinsidevolume.
Scott/Liu/WuAIM2024continuumrobotpaperPDFp5reportsbriefBadFitinterruptions;
geometryadjustmentwasattempted,notuniversalfix. OlderNDIV2manufacturer
guide(mirroredScribd) discussesmissingbeforevolumeedges/interference/motion
andtoolcablesnearFGcable. OfficialNDIbrochure environmentalaccuracycaveat.
Sources/locators/limitationsretained project_docs/NDI_TRACKING_LOSS_LITERATURE.md.

Inference:ourtowardFGmotionchangesdistanceANDorientation;notproveddistance-
onlyfailureoroutofvolume. No universalfailedangle/minimumdistanceassumed,
andV2docsnotexactPlanar2.0/V3qualification. No runtime/hardware/modelchanges.

## 2026-11-09 - Rotation test moves missing-transform association P3 -> P2

NEWDesktoplogSHA0f9f5079...,usernominal180degrotation,P1now+X;userexplicitly
confirmsorderP3/P2/P1.5279BXallCRCpass;final132.216s5277replies,6startup,
5271runtime. Baseallvalid,tip36missing/26episodes~25-50ms,all66.795-91.749s
inSECONDgroupP2. P3first0,P1last0;volume/partialflags0,systemstatus0.
Before rotationP3/P2/P1missing32/0/0;after0/36/0. NotfixedP3valveidentity
or+Xalone(P1positiveXtracks). Bothfailinggroupsmove towardlessnegative
trackerZ;supportsphysicalpose/orientationassociation,notZ-only/outvolumeproof.
Basequaternionchange~171deg supportsreorientation,notprecisemechanical180.

EarliermanualP3withoutvalvesreproducedloss;actuationnotnecessary. Sensor/
field/cable/motioncauseunresolved;do notweakenvalidity. ControllerSwitchNaNbug
remainsseparate. NewofflineclassifyCygnusStepGroups.py;report/raw/counts/
before-aftercomparison/annotatedplots/hashes cygnus_rotated180_steps_20261109_v1
underresidual_controller/Experiments/Model_Val. No runtime/SLX/deviceedit/
access,hardware/model/trainingrun,commit,push.

## 2026-11-09 - Reversed-order pressure steps localize losses to P3

NewDesktoplogSHA1976de55...,userP3->P2->P1individualmanualpressuresteps,
longinter-PMApauses.4978BXallCRCpass;final125.496s/4976reports,3startup,
4973runtime. Baseallvalid,tip32missingin30episodes(~25-50ms),all9.09-37.24s
inFIRSTP3group. ApproxresponsesP3[8.97,37.32],P2[57.10,83.81],P1[96.25,
121.22];missing32/0/0. Labelsfromuserorder/motiongapboundaries,noactual
pressuredataoronsyncclaim. No volume/partialflags,systemstatus0.

MissingrecursonP3transitions,plateausvalid. P3firststillfailslaterP2/P1none,
soassociationnotmerelylastchannel/elapsedtime. Priorno-valveP3manualalso
reproducedloss;valveactuationnotnecessary. Direction/pose/orientation/cable/
motion-dependentcauseunresolved. Validendpointholdsdonotvalidateintermediate
poseholds. Runtime/GRUvaliditygatechangesnotjustified.
Canonicalraw/events/PMAcounts/annotatedplots/hashes cygnus_pressure_steps_
P3_P2_P1_20261109_v1 underresidual_controller/Experiments/Model_Val.
Noruntime/SLX/deviceedit/access,model/hardware/trainingrun,commit,push.

## 2026-11-09 - Manual P3-direction motion reproduces tip loss without actuation

NEWDesktoplogSHAa6a8cd6c...,userreportsmanualP3-direction/no-valve-actuation.
3035BXallCRC/lengthpass;final3033reports/76.037s,2startupthen3031runtime.
Base0Amissing0,tip0Bmissing30/28episodes(26single,2double~25-50ms),allbase
valid/pairedframes,volume/partialflags0,systemstatus0. Runtime30.55-69.92s.
ThereforevalveactuationNOTnecessaryfortrackingloss;doesnotruleoutadditional
interferenceinotherconditionsoridentifyphysicalcause. Nooutvolumeclaim.

Somebase-frameadjacent-validposesoverlapearlierP3failures(nearestdistance
min1.12mm,median12.2,rangebroader);notmatchedtrajectoryoruniquesmallregion.
Earliergenericmanual0missesdoesnotcontradictthisnewP3manualreproduction.
Nextmanualdiscriminatorcouldholdaffectedposevskeepmoving;notperformedhere.
TrackerlossindependentofseparateSwitch/NaNcommandbug;bothissuesstillmatter.
Canonicalraw/events/CSV/PNG/PDF/FIG/hashes cygnus_manual_P3_20261109_v1 under
residual_controller/Experiments/Model_Val. Genericmanualplotterupdatedcounts.
No runtime/SLX/devicechange/access,model/hardware/trainingrun,commit,push.

## 2026-11-09 - Sine residual actuation test1: retained Switch causes NaN cycle

sinewave_FWD_INV_test_1_fwdPressure.mat model1.59/212s4241ticks,newhandoff
logs4outputs. Fullsequentialnewsource replay pressure1.49e-14bar,modes/ages/
limitflagsexact withcollection-inversenominal andactualA8feedback. Handoff
outputsALLfinite. BUTcommanddelaydoesnotgethandoffpressure. ItgetsoldSwitch:
perchannel ifhandoffpressure>0 choose rawcandidate else rawsine, then1tickdelay.
ThisreconstructionmatchesALLfinitevalues0barerror andNaNpatternexactly.

145.90realNDItiploss->handofffinitehold(mode4)->oldSwitchselectsNaNcandidate
->145.95NI-commandNaN->handoffmode6zero->nextfallbackfinite->coldcandidate
NaNselectedagain.435nonfiniteNI-commandrows;435mode6faults. Only5real
tipMISSING,434observerinvalidrowswithvalidsensors. Thereforecorruptedcommand
feedbackratherthanhundredsofNDIdropouts preventsrewarm. Actualanaloghardware
response toNaNunknown;NIproxiesnotmeasuredchamberpressure. Nohardwareaction.

Requiredmanualfix:Function9output1pressureCommand(oldwirenameuseResidual)->
UnitDelaydirect, BYPASSoldSwitch. Keep.05delay/Kill. Input3originalprotected
rawsinecommand,notMC3nominalreconstruction(.091barsecondarydifference).
Mode2/3arenotselectorbooleans. Do notjudgepressurehandofffromthismiswiredrun.
868/3267validinverseallocationsscaled3bar;remainingmodel/constraintissues
mayremainafterwiringfix. Frozenpressuremapprobesrho<.258notfullstabilityproof.
Canonicalofflineaudit plots/CSV/MAT/replays/hashes:
residual_controller/Experiments/Model_Val/sine_pressure_review_20261109_v1.
New analyzeSinePressureRun.m;COMMAND_HANDOFFclarified. On-diskSLXstaleearlier
Booleanhelper;actualnewsourceinferredbyexactreplay. Preserveoperatorunsaved
work. No runtime/SLXedit,model/hardware/trainingrun,commit,pushduringaudit.

## 2026-11-09 - New no-pressure manual Cygnus capture has no runtime loss

UseroverwroteDesktoplog.txt withmanualfast/no-pressuremotion, SHA7253f838...
(differentfrompreservedearliersine56204593...).2448BXCRC/lengthpass;final
2446replies/61.449s,first2startupmissing,ALL2444laterbase/tippairsvalid.
RuntimeMISSING0/0,repeatframes0,identicalconsecutivetransforms0,volume/partial
flags0,systemstatus0. Newcapturedoesnotreproduceearlier23P3tipmissingreports.
Userinsidevolumeclaimconsistentwithflags;nooutvolumephysicalcauseestablished.

Runtimequerygaps2x~50msand2x~175ms;bothlonggapsattheveryend61.10-61.45s
justbeforestop,notprovenmidrunfreezeorhardwaremissing. GUIrefresh/holding
behaviornotestablishedfromseriallog. ManualposecoveragediffersfrompriorP3;
cannotconcludepressurecauseslossorspeedaloneexplainsit. Trackerkeepsvalid
transformupdateshere. Parsernowhandleszero-missingfiles;newmanualplotter.
Canonicalraw/CSV/gaps/flags/plots/hashes cygnus_manual_20261109_v1 under
residual_controller/Experiments/Model_Val;NDI_DIAGNOSTICS_RESULTS updated.
No runtime/SLX/devicechange/access,model/hardware/trainingrun,commit,push.

## 2026-11-09 - Independent Cygnus raw log verifies P3-associated tip loss

Desktoplog.txt archived/audited offline:10369BXrepliesallheader/payloadCRCpass.
Finalsession10367reports,8startup;afterstartup10359reports,base0Amissing0,
tip0Bmissing23 in17episodes,1-2reportedupdates(~25-50ms),.222%runtime.
Allmissesbasevalid,systemstatus0,noout/partialvolumebits. Trackerexplicitly
reportsMISSING0x02;notjustMATLABangle/lengthvalidity. TimefromfirstfinalBX
usesreportedframefrequency320.264Hz,transformupdates40.033Hz,elapsed260.56s.
Misses155.722-209.302s;separate43longqueryintervalsnotcountedastransformloss.

Userconfirmsstandard3cycle/2repeatP1->P2->P3sine source,noMATpressurelog.
Motion-fitoffset14.625s alignssixgroups;ALL23missesininferredP3windows
[153.873,177.435]12reports and[187.435,210.997]11reports. Association
sequence/motioninferred,notindependentsynchronizedpressureproof. Poseson
loading/unloadingintermediateorientation,notonlymaxbend. Physicalcausenot
uniquelyidentified. Strongindependenttrackerloss evidence;handoffsmoothing
addressesconsequence,notrootNDIfailure. No validityrelaxation/interpolation.

Newoffline analyzeCygnusBxLog.py,plotCygnusSineAudit.m. Rawlog/decodedCSV/
events/PNG/PDF/FIG/sourcehashes canonical residual_controller/Experiments/
Model_Val/cygnus_sine_20261109_v1. NDI_DIAGNOSTICS_RESULTS updated.
No runtime/SLXedit,tracker/hardware/model/trainingrun,commit,push.

## 2026-11-09 - Updated residualCommandValid pressure handoff ready

Userrequestedcodewhilecollectingsinewave. Added residual_controller/
residualCommandValid.m, now7inputs/4outputs; output1pressure3x1 notBoolean.
BypassoldSwitch ->existingUnitDelay.05->protectedvalves. Inputs correctedp,
diagnosticvalid,originalfallback,actualA8command,t_sim,enable,reset. Actual
feedbackanchorssharedfinalslew;physicalKill/enable/stopremainpriority.
Defaults.10sgracehold,.50sfallback/reentryblends,6bar/s,0..3bar.
NoGRUgate/historychange. Repeatedtickscached,incorrecttickrestartmode7,
disable/reset/faultsemanticsdocumented. Hold/fadetimerssimulationtime,not
independenthardwarewatchdog. Guide COMMAND_HANDOFF.md exactportsandmodes.

validateResidualCommandHandoff syntheticstartup/dropout/reentry/invalid/
clock/disable/reset/actualfeedbackanchoringpass.7501frozenoldtest5inputs plus
simulatedcommanddelay showmaxordinarystep.30bar,0steps>.3 vsold30/max1.514.
Notnewarmresponseprediction:nominal/candidate/validityfrozen,GRU/plantnotrerun.
ReportplotsCSV/MAT/hashes command_handoff_validation_20261109_v1.
Actualmodelbuild/timing/physicaltestunverified;currentoperator sinewaverun
notmodified. NoSLXedit,model/hardware/acquisition/trainingrun,commit,push.

## 2026-11-09 - Operator test5 residual actuation improves tracking; fallback dips

circlesweep_FWD_INV_test_5_fwdPressure.mat model1.59completed375s/7501rows.
Onlyout saved.6143ready,1358notready.49NDItipMISSING/orientationinvalid,
baseinvalid0;9physicalintervals>.1s,max.6357s;13overrunsallorientationvalid.
No numericalguard;commandage0.37ready-loss episodes:32tipmissing,5timinggap.
Fallbackdurationmed1.65/max2.60simseconds due.6velocity+1GRUrewarm/restarts.
Selection+UnitDelayreconstructionexact;postKillsaturationNIcommand matches
selected4.45e-16bar,Kill1throughout. Eachfirstfallbackapplied.05s after
readinessloss.32eventsdrop> .1bar,max1.514bar;26stable-referenceepisodes
followedby>.1radphidip. Initialinvalid/timinggap precedesfallback;later
pressurestepcannotexplaininitialmissingflag,butdeepensphysicalbendingloss.

Example182.30tipMISSING ->182.35P1 2.516->1.002bar ->phi~1->.22rad
duringrewarm,ready184.70.37events33duringcircles. Onlycandidatebranchslew
limited;fallbackSwitchbypassesit. CurrentA5basevalidwiredtohelperstillwrong
readinesssource,butallbasevalid1andfinitecheckmatchedcontextvalidhere.

Prior test4 commonreferenceidentical0m. First3nominalvsresidualon lengthRMSEmm:
phi.2:[1.812,1.799]->[1.275,.917];phi.8:[6.961,6.975]->[4.328,3.220];
phi1:[8.203,8.173]->[3.875,3.923]. Separaterunsnotcontrolledrepeats.
Phi1MAE ready.157 vsrewarm.464rad;fullphi2.5RMSE[11.988,11.285]mm remains
poor. Fullcorrectedscaled1012/6143,slew122. Do notclaimall-amplitude success.

Nextproposal only:statefulcommandhandoff,briefboundedhold,rate-limitedblend
tofallbackandback;sharedfinalslewafterselection,physicalKillpriority.
Do notbypassGRUreset/readinessorholdoldrotatingcommandsindefinitely.
Canonicalreport/plots/events/CSV/MAT/sourcehashes:
residual_controller/Experiments/Model_Val/desired_pressure_test5_review_20261109_v2.
Newoffline analyzeDesiredPressureRun.m. V1incompleteKillalignmentmarked.
No runtime/SLXedit,model/hardware/trainingrun,commit,pushduringaudit.

## 2026-11-09 - V2 Switch/delay wiring checked

Savedv2read-onlyaudit: correctedpressure->Switchtrue, originalnominal->false,
validityfunction->control,thenUnitDelay.05/ICzeros3x1->valves allcorrect.
PostKill/saturation pressure+timestamp samplingfeedback retainscoherence.
Oneinputcorrectionneeded:Function9validcurrentlyFrom6/A5=NDIbaseorientation
valid,not desiredPressureDiagnostic.valid. BranchFunction8output3bus,select
valid,andreplaceFunction9input2. FinitecheckalreadyrejectsNaNwarmupbutbase
validisnotmodelreadiness. Audit residual_controller/wiring_review_20261109_v1.
NoSLXedit,load,compile,run,hardware,commit,push;fullmodelcompileunverified.

## 2026-11-09 - Direct-actuation algebraic loop identified

Operatorv2 screenshot andread-onlysavedSLX confirm desiredpressure->valves->
coherentsampledpressure/time->persistentobserverinputs6/7->context->inverse.
ZOH/RateTransition didnotbreaksame-tickfeedthrough. Manualguide
residual_controller/ACTUATION_CAUSALITY.md recommendsvalid+finite selector
(existingprotectednominal/startupfallback) thenUnitDelayTs.05,ICzeros3x1,
beforeactuators. KeepKill/enable/stoppriority;rawNaNpredictioncannotdrivevalves.
Command-pathdelay retainsactualpressure/timestampcoherence andexplicitlyadds
50msactuation/reference delay. No feedbacktimestamp-onlydelay,historybypass,
orclosed-loopstabilityclaim. No source/SLXedit,compile/model/hardware run.

## 2026-11-09 - Forward block dimensions and direct-actuation boundary

User connected desired/single-pass output toactuators; compile errorreported
armS_single_dynamics concatenation and MATLAB Function1 unresolveddX.
SavedSLXread-only:forward symbols inheritedsize-1,IntegratorIC4x1; exact
firstwrongoperandnotestablished. Added paste-ready RL_data_collection/
rlForwardModelBlock.m, same3inputs7outputs/constants/EoM, explicitpressure3x1,
state4x1,torque2x1,derivative4x1. Shared dynamicsunchanged. ManualSymbolsizes
in FORWARD_BLOCK_SHAPES.md. Numeric50cases100row/column variants pass,
originalRHSdifference8.882e-15. Actualmodelcompile notperformed.

Important integration boundary:desired-pressure source returnsNaN while
cold/invalid; directvalveconnection lacksstartup/validity/hold/faultpolicy.
Need supervisedtransition preservingphysicalKill beforeactualactuation.
No validitybypass orautomatichardware/modelrun. NoSLXedit,commit,push.

## 2026-11-09 - Desired-motion pressure block source ready

User authorized desired-length pressure prediction. New paste-ready
rlResidualDesiredPressureBlock.m replaces observed-motion demand with full
DESIREDq/v/a nominal force.5inputs/3outputs, existing diagnostic bus; input2
NOW6x1[l2d;l3d;dl2d;dl3d;ddl2d;ddl3d], branch existing collectionMux6 output.
GRU stillactualq/v,p,dt,immutablehistory,1eval/no commits. Frozen-residual
single-pass approximation, no exact pressure-dependent inverse claim.
Collectionminimum-pressure mapping,rails/slew/readiness/NaN gates retained.
Manual DESIRED_PRESSURE_BLOCK.md; historical blocks/PAC/PDFB unchanged.

Numeric64/128synthetic and5433test4rows(5118valid/315invalid) pass. Newnominal
equalscollectionrequest4.45e-16bar; desiredforce3.56e-15N,recordedactualresidual
1.07e-14N;independentpre/post1.34e-15/1.12e-15bar. First70s correctedpeak
[.523,.468,.484]bar vsnominal[.242,.240,.242],no slew. Fullrecordcorrected
scaled126,slew3027;larger demands hit constraints anchored toactualcommand.
Modelbalanceerror afterlimits retained (not physicalforce measurement).

Coder emittedC but localstandalonecompiler existingIndexexceedsbounds failure.
IndependentGCC/native764fixtures passmax2.203e-13. Canonical
residual_controller/desired_pressure_validation_20261109_v1 and
desired_pressure_native_20261109_v1 retainreports,plots,CSV/MAT/sourcehashes.
Operator must paste fullfile, change input2 size/wire, keep discrete.05s,
reuse ResidualSinglePassDiagnosticsBus and log prelimit/flags. Initial manual
installation observationonly; no valve/command/Kill changes. No SLX/model/
hardware/training run,commit,push. Actualmodelcompile/timing remainsoperatorwork.

## 2026-11-09 - Circle test4 inverse comparison uses different state targets

circlesweep_FWD_INV_test_4.mat: model1.54, operator StopCommand271.6s,
5433 rows/5118valid contexts. Only out saved. Collection nominal at desired
q/v/a reproduces tauDesired3.56e-15N and requestedpressure4.45e-16bar;
actual NI=requested=context exactly. Current single-pass target instead
uses measuredq/v+desiredacceleration, no desiredq/v input. Full source replay
matches outputs1.12e-15bar; this is an observed-motion reconstruction, not
the desired-trajectory inverse the user intended to compare.

Firstcircles13-62.95s: desiredpeak[2.600,2.600]mm,measured[.374,.180]mm,
trackingRMSE[1.777,1.806]mm. Actualpeakp[.242,.240,.242]bar versus corrected
[.258,.259,.252]; correctedvsactualRMSE[.00660,.00873,.00694]bar.
No corrected slew/scaling across5118valid rows; similarity not clipping here.
Offline desiredforce minus SAME frozenactualresidual gives peakp
[.523,.468,.484]bar: illustrates missing target, NOT validated requiredpressure.
Pressure-dependent residual at this candidate leaves modelbalance error
median/P95/max[1.478,2.233,3.569]N, not physicalforceerror.

Matched50ms forwardRMSE nominal/res/CV mm[3.751,3.807]/[.226,.231]/
[.055,.058]. Local measured-state prediction accuracy does not establish
desired inverse accuracy/deadzone identification. Next: separate shared
nominalForceDemand_N input (collection tauDesired, or complete PAC demand),
subtract actual-state/history residual; maintain constraints and diagnostics.
Do not replace causal measured GRU stream with desired states.

Report/PNG/PDF/FIG/CSV/MAT/source hashes:
residual_controller/Experiments/Model_Val/target_mismatch_review_20261109_v1.
New offline analyzeSinglePassTargetMismatch.m; guide clarification added.
No runtime/NN/SLX edits, model/hardware run, training, commit or push.

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

## 2026-11-06 - Residual real-time motion/error presentation exports

User explicitly chose circlesweep_test,helix_FWD_INV_test,sinewave_test and
sinewave_test_100g. Created separate angle/error and length/error figures
perfullrecord, plus firstcompletephi-reference1circle159-169s andhelix5-65s
firstepisode.12stems PNG4320x2430/PDF/SVG/FIG in residual_controller/Experiments/
Model_Val/presentation_shadow_20261106_v2, adjacentZIP; README/FIGURE_INDEX.
Predictions are previouslyissued50ms matchedones, not free runs or controller
tracking; errorprediction-minusmeasured, theta circular/nearstraight masked,
sensor/warmup gaps preserved. Allrawq/error/priorissue identities exact;
measuredanglegeometry<=9e-16rad. Full CSV/MAT/metrics/sourcehashes retained.
No fitting,SLX/controller/raw changes,hardware/modelrun,commit or push.

## 2026-11-06 - Requested PAC candidate comparison figure

Created one six-panel presentation figure for candidates1/5/8/9: theta/phi
tracking, signedangularerrors and l2/l3errors. Candidate1 maps to firstbaseline
r01;9 maps smooth_candidate09, explicitly labeled. Common recorded20-140s
window shown0-120s; references match,no smoothing/invalidfill. Allfour full/
late metrics reproduce archived evaluations<1e-9. OutputPNG/PDF/SVG/FIG and
CSV/MAT/sourcehashes under passivity_adaptive_control/experimentations/
autotuning/no_payload_phi1_T10_20261027/
presentation_candidates_1_5_8_9_20261106_v1. See README/FIGURE_INDEX.
Offline plotting only; no SLX,raw,controller,hardware,commit or push changes.

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

## 2026-11-05 - Circle-sweep leakage concern audited

circlesweep_test.mat375s/model1.51;7501rows,6767prospective scores. Full causal
observer replay matches issued/matched predictions1.39e-16m,force1.07e-14N,
allstatus/ready/masks exact. FutureNDI target+5/-3mm at87.10s leaves prior
prediction unchanged, errorchanges-5/+3mm. Nov3packetSHA unchanged;101context
probes also reproduce. No future target leak found. Current measuredq DOES
initialize each new50ms forecast; one-step overlays are not whole-circle
free simulations. Even persistence beats residual at50ms here.

AllRMSEmm residual[.748,1.220],nom[4.296,4.714],CV[.097,.111],persistence
[.206,.256]. Moving2042rows:res[1.337,2.204],nom[6.321,7.328],CV[.174,.199],
persistence[.375,.464].320no-reset anchors,88moving, conditional recorded
futurep/dt only: .5s moving res[1.660,2.906],nom[9.102,9.902],CV[1.748,2.042].
No guardfailures; not a live forecast or uniform CV advantage reversal.
VelocityfeaturesmaxZ18.05/17.13;12.84%readyanyZ>5, normalized-range diagnostic
not proof of support violation. Review velocity coverage/model pressure
sensitivity before NMPC. SavedScope input3 is separate model/length chain,
not matchednominalcolumn2; logged residual input is correctly prior-issued.
Report/data/figures/sourcehashes: residual_controller/Experiments/Model_Val/
circlesweep_test_review_v1. Offline analyzers only; no SLX/controller/raw-data
edits, hardware/model runs, training, commit or push.

## 2026-11-05 - Sine-wave presentation plots, no payload and100g

Operator sinewave_test.mat and sinewave_test_100g.mat both212s,model1.51,GRU128,
logged0.05s horizon. New plotResidualSinePresentation.m exports five16:9figures:
fullrun perpayload,matchingcycle2/channel close-ups perpayload,channel2cycle2
side-by-side. Only requested measuredNDI, nominal and residual+nominal curves;
mm units, aligned previously issued predictions, shared scales and invalid gaps.
PNG3840x2160/300dpi,PDF/SVGvector,FIG and data/metrics/sourcehashes retained in
residual_controller/Experiments/Model_Val/presentation_sinewave_20261105_v2.

IndependentNDI/contextNI match, previous-issued identity, recalculated errors
all exact. Scoredrows3866/4241unloaded,3984/4241loaded. Actual sine phase5,
phase4relaxation; fixed error-independent cycle2 windows12.9–20.7,80–87.8,
147.15–154.95s forchannels1/2/3. All-scored RMSEmm nominal/residual:
unloaded[3.696,2.897]/[.482,.429],100g[4.505,3.733]/[.729,.809]. Active-sine
metrics separated in report. Labels fromuser/filenames; calibration/payload
attachment timing not embedded. No free-run/control improvement inferred.
No SLX/controller edits, hardware/model run, fitting, commit or push.

## 2026-11-04 - Operator pressure-combine forward-model validation

pressure_combine_test.mat (residual_FWD_INV_validation1.50) completes191s with
no simulation error. GRU128,0.05s issued predictions verified one prior tick;
raw NDI q/NI commands match observer exactly. Frozen-context101probes reproduce
neural force/forward prediction to2.66e-15N/1.99e-17m. No accidental current-q
comparison.3814/3821orientationvalid,baseallvalid;3560prospective scores93.17%.
Seven invalid NDI samples excluded/re-warmed, no numeric-guard status. Two
>0.1s physical poll gaps at startup only; typical physical interval.050807s.

All-scored length RMSE mm: residual[.234,.198],nominal[1.007,1.569],CV[.046,.047].
Moving(normv>.001m/s),1043rows: residual[.346,.300],nominal[1.700,2.491],
CV[.074,.074]. Residual beats nominal substantially, but CV beats residual at
this short horizon. Not tracking-error/force-error/inverse-accuracy proof.
Only forward/context buses logged; no learned-inverse diagnostics. Operator
model build/run is now evidenced, actual worst-case block timing unmeasured.
Source shadowOptions itself not saved; save options on next trial. Pressure
0..3bar, feature maximumZ4.118, neural force peaks10.4/10.8N.

Report/plots/data: residual_controller/Experiments/Model_Val/
pressure_combine_test_review_v1. Next recommendation: frozen longer-horizon
causal rollout comparison and disconnected learned-inverse observation before
further residual-on PAC/PDFB. PAC dropout recovery still unresolved. No
SLX/controller changes, hardware run, fit, commit or push during this audit;
only reusable offline analysis and docs/artifacts added.

## 2026-11-03 - Passive residual validation block code ready

User clarified add-block means source only; all SLX files untouched. New
rlResidualShadowBlock.m and optional rlResidualInverseShadowBlock.m under
residual_controller, setupResidualShadow.m and SHADOW_VALIDATION.md provide
exact ports/buses/manual wiring for RL_dataCollection_v2. Forward uses raw
orientation-inferred lengths, independent training-matched causal velocity,
physical pollStart and tip frame, final post-Kill commands and separate sim
command timestamps. Prospective0.05s predictions and retimed replay/errors
are distinct. Strict invalid segmentation and0.6+1s warmup; no control output,
fault latch or automatic actuator rearm. Inverse uses immutable pre-update
context, fixed trained nominal LHS and diagnostic-only pressure predictions.
No PAC recovery fix or protective-limit change included.

Final numeric verification shadow_validation_20261103_v4 passes both widths,
exact trained flow/cell, causality, wrapper agreement,1/2/3-invalid recovery,
duplicate frames/ticks, gaps/clocks, stale pressure and inverse consistency.
firstTime1935rows replay:1825warm predictions; observer ready again23s after
recorded PAC fault20.8. Pre-fault prospective length RMSE aug.694/.622mm,
nom1.537/1.747mm,CV.111/.098mm; mostlystartup/trim/blend, not tracking proof.
Postfaultstationary metrics separate. Observer-only numeric P95.798ms.
Coder emitted C but hit known standalone toolchain error; independent GCC
compilation and1326native ticks pass, maxdifference5.13e-12. Actual SLX build
and full real-time timing untested. Preserve original calibration/base for
fixed trained gravity; no interpolation, fitting, hardware, commit or push.
Next: user installs forward block first, verifies source clocks/NI branches,
logs predictions/errors/validity during an existing manually run experiment.

## 2026-11-03 - First PAC residual runtime fault diagnosed offline

Operator firstTime.mat: residual enabled throughout, variable common enabled.
Tip orientation invalid at20.70/20.75s; recovery20.80 gives0.15s accepted-sample
gap, exceeding residual maxGap0.10s. History resets/ready=false; active-anchor
ready guard makes recovery unusable and invalidGrace0.10s latches status5.
Inverse was converging immediately beforehand; not a tracking-error trip.
Offline pressure/helper replay matches status/ready exactly and pressures
to8.53e-14bar. Original-supervisor counterfactual statuses4/4/3 versus4/4/5.
This exposes a two-missing-then-recover gap in prior synthetic coverage.
Retain common preload; variable-common choice optional but matched off/on.
Need explicit recovery design and shadow forward/inverse checks, not bypassing
readiness or merely increasing limits. measurementTime is latched simulation
time, not physical NDI elapsed time; settle separate/aligned clock contract.
Saved PAC_residual.slx still has old PAC chart; operator runtime edits appear
unsaved. Detailed report/figures: residual_controller/Experiments/PAC/20261103/
firstTime_review_v1/README.md. No controller/model changes or hardware run here.

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



## 2026-11-03 - Residual controller branch handoff

User requested committing the completed residual-learning work and creating
codex/residual-controller-integration. Handoff: RESIDUAL_CONTROLLER_HANDOFF.md.
Start PDFB integration, then PAC after resolving adaptive compensation overlap;
perform offline checks before preparing operator-run matched comparisons.
Defer HIL learning automation until after the initial controller comparison.
GRU128is the short-horizon candidate; GRU64retained for stronger circle free runs.
Source/results commits preserve operator SLX files/staged rename, local raw
recordings/derived datasets/checkpoints and unrelated paper-draft files.
Versioned/local artifact inventory: residual_dynamics_learning/VERSION_CONTROL.md.
No controller implementation, hardware, SLX edit, artifact migration or push
is included in this checkpoint/branch operation.



## 2026-11-03 - HIL paper reviewed: proposed supervised residual-model loop

Reviewed Marquez et al2023 local HIL/TD3 paper. It learns a locomotion POLICY,
not residual dynamics; its reported gains do not establish model accuracy.
Proposal: reuse hardware-data/replay/batch-update structure for supervised GRU
model learning, contiguous recurrent sequences and held-out prediction tests,
with PAC/PDFB retained. Start between trials; neural weights frozen within each
comparison. No actor-critic or ROS rewrite required for this proposed version.
Detailed source/limitations/sequence: project_docs/HIL_RESIDUAL_MODEL_LEARNING_NOTE.md.
Primary next step remains consistent PAC/PDFB model integration/comparison;
no implementation, hardware run, new fit or automatic deployment authorized.



## 2026-11-03 - Three-seed GRU capacity comparison completed

User authorized quick offline1x64/1x128/2x64 comparison. Three seeds113/211/307, matched original1200-update stage plus1200 old-only continuation updates; unchanged physics/loss/masks/normalizer.
Six new stage1 and eight new continuation fits; exact existing64-unit fits reused by hash. Old TRAIN gradients only; old/Nov3 holdouts select checkpoints and architecture; Nov2 evaluation excluded from selection.

| Architecture | Parameters | Mean score (mm) | Sample SD | Seed113 | Seed211 | Seed307 |
|---|---:|---:|---:|---:|---:|---:|
| 1x64 | 14162 | 0.58531 | 0.02495 | 0.58716 | 0.55949 | 0.60928 |
| 1x128 | 52882 | 0.48684 | 0.02226 | 0.46154 | 0.49553 | 0.50344 |
| 2x64 | 38930 | 0.53814 | 0.05209 | 0.50814 | 0.50799 | 0.59828 |

Predeclared winner: **gru128**. Frozen representative candidate: candidate_force_model.mat.
SHA256 `8f2e0bb55ea9399b62abd5abe0b8f6dd15ab0470a400ad3e593e1edcd8f9f545`.
The existing fine-tuned GRU64 research release remains unchanged. This candidate is an offline result, not deployed or controller-qualified.

Representative seed113 fixed in advance, not chosen by best seed. Current fine-tuned GRU64 release retained; separate frozen candidate archived.
GRU128 improves 4/10 fixed long examples relative to GRU64, while others worsen; these correlated examples are not independent statistical trials. Circle free RMSE worsens on Nov2 0.623->0.930mm and Nov3 1.005->1.366mm. It improves several sine/helix examples. The short-horizon winner is therefore a candidate for measured-state model augmentation, with mixed long-horizon accuracy. Retain the current GRU64 reference; do not infer improved physical tracking from short prediction scores.
All nine cases have common-anchor evaluation; three representatives have whole-record one-step, fixed long motion predictions, and208-case bounded numerical inverse checks.
Report: residual_dynamics_learning/force_model/capacity_comparison/campaign_20261103_v1/REPORT.md.
Figures FIG/PDF/PNG, predictions/checkpoints/source/hashes and verification retained. No physical force/closed-loop qualification; fresh confirmation and PAC/PID FB integration pending. No hardware/SLX/migration/Git changes.


## 2026-11-03 - Force GRU64 fine-tuning completed

User authorized fine-tuning. Two matched1200-extra-update arms from exact frozen
GRU64/B seed113: old October TRAIN only;16old/16Nov3 stream replay. LR1e-4,
fresh Adam, unchanged architecture/.5s horizon/burn-in/loss/filters/normalizer.
Nov3 whole-trial split:28train/9selection,2264new train chunks,80new holdout
anchors; original128old selection anchors/splits/59,630-row normalizer unchanged.
No interpolation. Both candidates improve old and new internal scores.
Predeclared .5old+.5Nov3 score (each .7oneStep+.3halfSecond), old regression
guard10%, chose old_only BEFORE Nov2 evaluation: composite.71126->.58716 versus
mixed.59696. Old internal.54411->.41699; Nov3.87840->.75734. Mixed improves
new score more (.68348), but old score less (.51044). Keep both candidates;
this is one optimization seed, not evidence of universal old-only superiority.

Selected research model: residual_dynamics_learning/force_model/fine_tuning/
campaign_20261103_v1/release/force_model.mat, SHA256 1c38edc36ba6ff4a41da8f974a88ca1c5ea775655426da91124a35b255a9c870.
Full report: same campaign/REPORT.md. Selected gradients use older TRAIN data;
Nov3 holdouts were used for checkpoint selection, not gradients in that arm.
Nov2 external one-step macro vector RMSE .53588->.36143mm; mixed.44580mm.
Nov3 heldout one-step .93307->.75485mm; mixed.69078mm. Values use fixed original
seed113, not earlier three-seed mean. CV remains better one-step.
Full circle q-vector RMSE: Nov2 .9440->.6229mm; Nov3 1.2636->1.0047mm.
Selected model improves allfive Nov2 motion examples, but some Nov3 long
windows regress slightly; mixed replay has larger long-window regressions,
despite better new one-step scores. Allnegative results retained. Nov3 helix
comparison uses HELDOUT trial2 (4.0294mm original), not prior TRAIN trial1
(2.8787mm); do not mix those comparisons.30 complete long model/window paths.

75one-step and30long-window table rows independently recomputed<5e-15mm;
source/model/data/splits/normalizer and original releases preserved. Selected
inverse self-consistency208/208, max equation error3.57e-5N, state-commit and
bounds pass; no physical force/control qualification. Seven FIG/PDF/PNG bundles,
full CSV/MAT/checkpoints/optimizer states/source hashes archived. Current
controller/PAC adaptation integration plan still pending; no SLX/hardware,
migration, commit or push performed. No fitting or checkpoint selection on Nov2.

## 2026-11-03 - Primary goal clarified: model augmentation for PAC and PID FB

User explicitly intends to use the improved dynamics MODEL with existing PAC
and PID feedback controllers and test whether tracking improves because the
model is more accurate. This is the primary integration target. Detailed plan:
project_docs/RESIDUAL_CONTROLLER_COMPARISON_PLAN.md. Four primary conditions:
PAC nominal / PAC augmented; PID FB nominal / PID FB augmented. Compare model
off/on within each controller with matched gains, references, filters, pressure/
common-pressure allocation, load, calibration, initialization and duration.
Freeze neural weights; recurrent state evolves. Match PAC adaptation settings
and reset policy. Report tracking/convergence/smoothness/limiting and independent
model accuracy, not just one-step prediction error. A separate retuned-controller
study can follow but cannot isolate model benefit by itself.

Before integrating, reconcile the residual's fixed nominal K0/D0/B*p definition
with each actual controller. PAC already estimates K/D/hysteresis/bias; derive
consistent compensation/adaptation before stacking the neural force term.
Pressure-dependent residual requires bounded inverse calculation with fixed
pre-current hidden state during candidate search and one selected-state commit.
Verify actual PID/integral path (repository also uses PDFB/PDI/FB names), actual
pressure/common-pressure operating coverage and runtime settings. Original
passivity is not automatically inherited. This is a recorded design/test plan,
not a completed controller derivation, SLX change or hardware authorization.
No continuation training, controller edits, hardware, commit or push this turn.

## 2026-11-03 - November3 frozen-model plots and continued-training recommendation

User requested same circle/sine/helix plots for November3 and asked which newer
day to use for continued training. Complete seven FIG/PDF/300dpiPNG bundles,
MAT/CSV/source/hashes: residual_dynamics_learning/figures/
motion_predictions_20261103_v1/frozen_model_predictions_v1/REPORT.md.
Same frozen FORCE GRU64/B seed113, no fitting, changed normalizer or imputation.
Separate measured-state one-step and free recorded-input predictions, nominal
shown separately. New helper reproduces ALL saved November2 circle one-step
anchors/outputs<1e-12m; all15 initial free steps match measured replay<1e-10m.
All15 model/window free paths complete, four helix labels stay NaN/unscored.

Free q-vector RMSE mm, GRU/CV/nominal: circle1.264/1.948/11.625;
sineP1 1.718/9.767/6.742; P2 1.379/15.735/9.700;
P3 1.688/19.562/7.670; helix2.879/13.169/11.808. Correction to an intermediate
commentary range: sineP2 is1.38mm, not1.69mm. Selected source intervals:
circle phi_reference1,159-169s; sine[6.6-19.2,102-114.6,197.4-210]s;
helix6.6-64.95s (59.047physical seconds). These are chronological examples,
not independent repeated trials. One-step CV still better; reconstructed PCC
XYZ is not sensor XYZ. No controller tracking or unique physical deadzone claim.

Logged Nov3 circles10s versus Nov2 20s; sine~12.575s versus~31.4s; helix10s/turn
and peak requested phi2.5 versus20s/turn,phi1.5. Speed, amplitude and date all
change, so larger circle/helix errors cannot be attributed to speed alone.
Actual59,630 old TRAIN feature rows audited: some new circle/helix q/v extend
beyond old observed extrema (~2-3.5% per affected feature); q2 minimum old
-14.98mm versus newcircle-22.18mm/helix-20.73mm, positivev2 old.01322m/s versus
new~.0209m/s. Large training z-scores alone do not prove outside-range input.

Recommended NEXT offline experiment, not executed: fine-tune frozen GRU64 on
old October19/27 TRAIN plus November3 TRAIN, retaining older data replay and
whole-trial/condition internal holdouts. Keep November2 evaluation-only for
frozen before/after comparison, acknowledging it is previously inspected, not
blind confirmation. Lower learning rate, same architecture/loss/history and
original normalizer initially; declare fresh optimizer state (release weights
are not a full Adam resume). Compare extra-old-only updates against old+Nov3
at matched budgets, to separate new data benefit from just optimizing longer.
Do not combine imputation/filter/capacity changes in the first continuation
experiment. If Nov2 enters training it loses held-out status; use another
cohort/fresh day for independent evaluation. A future fresh repeat/day is needed
for final confirmation. Neither current day establishes2s-circle performance.
Sources/model/raw/old Nov2 plotting archive unchanged; no hardware, SLX,
training, migration, Git commit or push performed.

## 2026-11-03 - Requested November2 motion plots; November3 NDI gap audit

User clarified main plotting cohort is November2. Frozen FORCE GRU64/B seed113,
not the older position-correction GRU. Seven FIG/PDF/300dpiPNG sets plus full
MAT/CSV/source/hashes: residual_dynamics_learning/figures/
motion_predictions_20261102_v1/REPORT.md. Circle phi_reference1,259-279s;
one sine period/channel; helix6.6-124.95s (119.85physical seconds), loading and
unloading after observed warmup. Separate exact saved measured-state one-step
and free recorded-pressure predictions, nominal shown separately. Fifteen free
model/window paths complete; first steps reproduce frozen results<1e-10m.
Length-vector RMSE mm, GRU/CV/nominal: circle.944/2.253/11.481;
sineP1 2.336/9.481/6.070; P2 1.778/17.480/9.656;
P3 1.796/21.731/7.559; helix.763/4.740/8.972. Chronological examples exceed
.5s training and prior1s evaluation; this long-window GRU advantage coexists
with CV winning short-term. No universal or controller tracking claim. PCC XYZ
is reconstructed from inferred q, not independent sensor XYZ. Six missing
labels remain NaN/unscored; free model carries on known pressure/time without
measurement injection. One-step retains original reset/recovery warmup.
Source/model/raw/export hashes unchanged; previous release preserved.

All five November3 files exported/audited separately, no fitting or new model
score evaluation.57 invalid NDI rows/55 episodes:53 one-frame, two two-frame.
Separate explicit preview fills53 one-frame gaps at<=.12s valid-to-valid
physical span, keeping raw values/validity and imputed masks, q1=-q2-q3; both
longer episodes remain missing. Artificial-mask proxy reconstruction vector
RMSE.040-.098mm acrossfive records (max individual one-frame error.4855mm).
Actual dropout truth remains unknown; overlapping masked blocks are correlated.
Artifacts/source/CSV/FIG/PDF/PNG in figures/motion_predictions_20261103_v1/
short_gap_audit/REPORT.md. Preview never entered training or frozen evaluation.
Proposed future training experiment: bounded offline inputs, explicit imputed
flag, measured-only targets, long-gap resets, unchanged measured-only/causal
validation. Two-sided interpolation uses the future endpoint: masking only
filled targets does not remove leakage through recurrent inputs/state.
Alternative: gap-aware physics carry with loss only at real observations.
No interpolation policy deployed, retraining, SLX/hardware, migration or commit.

## 2026-11-03 - Streaming force-GRU sequence completed

Explicit sequence executed: six architecture/horizon fits, one burn-in comparison,
four additional seed fits. New observed-state TBPTT primary path and independent
short free-prediction state; no future target inputs. Old October19/27 whole-trial
splits and normalizer unchanged; all finalists frozen before November2 replay.
Selected B,64-unit single layer,0.5s horizon,.5s burn-in;
fixed representative seed113, chosen config by three-seed mean older score.
Two-layer32+32 did not win under the common1200-update budget; not a general
depth result. November2 macro one-step 0.5629 +/- 0.0684mm; .5/1s path
0.3613 +/- 0.0148/0.4272 +/- 0.0183mm. Matched CV one-step
0.0621mm. CV remains better at all tested horizons through1s. All numerical
guard failures 0; per-phase/coordinate/
force-rate/seed results retained. New mode differs from parent release, so do
not attribute all changes to width. No post-test tuning or new November3 pooling.
Offline bounded, pressure-dependent model inversion: 128/128 known
targets feasible, max equation residual 4.743e-05N, median/p95 latency
5.19/13.93ms excluding acquisition/control scheduling.
Trial hidden states remain fixed; one selected-state commit passes. This is
self-consistency, not physical force or closed-loop tracking qualification.
Code/state/gradient/horizon/split/hash checks and independent saved-metric
recompute pass. Prior release unchanged. Reports/weights/CSV/MAT/FIG/PDF/PNG:
residual_dynamics_learning/force_model/streaming/REPORT.md. No SLX/hardware,
migration, commit or push.

2026-11-03 user proposes two GRU layers and <=1s training rollout for current-
force pressure correction. Proposed test plan recorded in
residual_dynamics_learning/force_model/NEXT_EXPERIMENT_PLAN.md; no new training
or runtime edits here. Align measured-state streaming/recurrent-state update
contract first; current release's warm-up+free prediction differs from live
measurement-conditioned estimation. Six-model screen:1x32,1x64,2x32 crossed
with .5/1.0 PHYSICAL second training horizons, explicitly weighted one-step
and short-rollout length loss, matched data/objective/update budgets. Then
burn-in .5vs1s and repeated seeds for finalists. Select using development data,
retain November2 evaluation only. Pressure-dependent residual still requires
bounded inverse/common-pressure policy and one committed state update/sample,
not one per pressure candidate. No hardware/SLX/model replacement/migration.

2026-11-03 user plans LATER migration of heavy residual-model/data artifacts
outside OneDrive to C:\Lab\Research\Residual modeling, preferably while no
agents/training are active (possibly at night). Deferred: no move or scheduled
task executed now. Future migration must verify destination paths,
copy/check hashes and preserve provenance/reproducibility before cleanup.
User is collecting additional circle speeds/phi and step rise times for later
model development. Files/settings not yet audited; do not assume the new runs
are complete or pool them automatically. Retain current October19/27 training
and November2 evaluation-only contract for the frozen release. For a later
campaign, declare training/held-out runs/conditions before fitting; compare
capacity/history/loss settings on development data, not November2 scores.

2026-11-03 explicit user authorization executed: NEW generalized-force GRU32
trained using October19+October27 only; November2 reserved for evaluation.
Canonical model/report: residual_dynamics_learning/force_model/release_20261103_v1.
It outputs two forces in N on the RHS of the corrected nonlinear mechanical
equations; lengths result from integration, with no learned position jump.
New old-TRAIN-only normalizer,2,883 training windows,160 whole-trial internal
selection windows; day2 step only3 train windows and combinedzero due gaps.
Integrated-length loss includes residual-force-rate/magnitude regularization.
No force/derived-acceleration truth labels or new-data normalization/fitting.
Compared two completed fits using older holdouts only; lowerlr.0003 epoch38
selected, old macro2.1103mm. Exact release frozen before evaluation.
November2:50,078 anchors; macro one-step nominal4.7352/CV.06235/forceGRU2.1510mm.
At2/5s path: nominal7.0438/6.6372, CV.8927/3.4062, forceGRU1.8375/2.1267mm.
All540 matched rollout paths avoid numerical guard. Force model improves poor
nominal and CV at5s, but weak short-term accuracy prevents qualification for
inverse/PAC/PDFB. Do not hide that limitation or equate new force scores with
old position-correction GRU. OI-11 retains pending development. Current pressure
enters the residual, so inverse allocation must solve that input dependence.
Differentiable exact-PCC polynomial coefficients, implicit SDIRK solver and
first-order numerical solver-Jacobian training checks pass; no claim of physical
force identification/stability. Reports, weights, old-only split/source hashes,
CSV/MAT/FIG/vectorPDF/PNG retained in force_model. No post-test refit, old-data
overwrite, SLX/runtime/hardware/commit/push. November2 previously inspected:
held-out evaluation, not blind confirmation. Earlier proposal-only status below
is superseded for first force-model training, not controller deployment.

2026-11-03 user asks to retain the spike-reduction idea for later. Record in
RESIDUAL_DYNAMICS_LEARNING_PLAN.md: diagnose nominal versus residual origin;
prefer multi-step integrated-length training and modest temporal smoothness
regularization of the intended force residual. Optional causal output low-pass
must be evaluated with its delay, not merely smooth displayed curves. No
retraining, filter/runtime edits or new model implemented for this discussion.

2026-11-03 requested SineWave4 visual comparison completed, same five existing
predictors. CV+GRU grouped separately from nominal/ridge/kinematic-ridge.
One-period free-running windows chosen first chronologically within mature
uninterrupted segments: P1 trial1 40.8-72.25s, P2 trial3 214.6-246.05s,
P3 trial6 527.35-558.8s. No measured-state resets within a window. q2/q3
path RMSE mm, P1/P2/P3: CV9.604/16.748/20.685; GRU.951/3.899/2.401;
nominal+ridge3.376/2.724/4.714. GRU has spikes; ridges better on P2 here.
Full-record one-step plots explicitly separate from free-running plots.
Four FIG/PDF/PNG sets, source/CSV/MAT/hash/verification/report retained under
residual_dynamics_learning/figures/sine_predictors_20261103_v1. All15 complete;
first rollout steps match existing frozen predictions, independent score
recompute passes. No training/newpredictors/controller/hardware/commit/push.
Current learned outputs are position corrections; force-residual work pending.

2026-11-03 user requested circle visuals for the five listed predictors only.
Retained under residual_dynamics_learning/figures/circle_phi1_predictors_20261103_v2.
Logged phi_reference=1 selects trial3,five20s circles,259-359s; measured phi
.052-.230rad. First complete circle259-279s (20.23physical seconds) predicted
from one measured history, no state resets. q2/q3 path RMSE mm: nominal11.489,
constantvelocity2.253, nominal+ridge3.825, kinematic+ridge4.750, nominal+GRU1.096.
GRU closest length/derivedXY for this onecircle, not universal winner. XY is
PCC reconstructed from inferred/predicted lengths, not independent sensorXYZ.
Separate one-step plots across allfive circles are explicitly teacher-forced.
Source/CSV/MAT/hashes and FIG/PDF/PNG retained. No fitting/newpredictors,
controller/hardware/Simulink/commit/push. Generalized-force residual still pending.

2026-11-02 fresh offline validation complete: canonical output is
residual_dynamics_learning/validations/fresh_20261102_v2/REPORT.md.
StepWave/SineWave4/combinePressure/circleSweep/helix supplied exports exactly
match current raw MAT re-exports. 50,078 guarded anchors, frozen October19
normalizer/weights unchanged; no fitting or new model selection. Mean per-run
q2/q3 one-step vector RMSE: nominal 4.9183 mm, constant velocity .06235,
ridge .05976, kinematic ridge .05961, GRU .34447. At 2/5 s, ridge path RMSE
.7268/1.6204 mm versus constant velocity .8930/3.4076; GRU 1.1265/1.5112,
kinematic ridge .6795/1.8064. Sixty deterministic windows per horizon/model,
no measured-state resets, known recorded pressures/query times. All 1,440
propagations stay within the inherited numerical guard; no safety guarantee.
One-step sine ridge is slightly worse than constant velocity; phase-level
and per-run scores retained. Physics-prior advantage remains mixed; inverse
pressure-only AR still matches full inverse ridge, confirming shortcut concern.
Mode0 frozen context retained; corrected mode1 nominal separately tested.
Recorded model versions differ (1.46 step/combined, 1.47 sine/circle/helix).
These are next-position corrections, not force residuals. Next implement
history-conditioned generalized-force augmentation with integrated-q training
and explicit corrected nominal version; no further calibration detour or
controller integration justified yet. November2 is now evaluated, not untouched
confirmation for later changes. Raw exports, diagnostics, CSV/MAT/source hashes,
FIG/PDF/PNG retained; v1 incomplete attempt marked DO_NOT_USE. No training,
SLX/runtime/hardware execution, commits or pushes by assistant.

2026-11-02 SineWave4 reviewed: new same-BXdiagnostics logged/aligned all12721
samples, replycounter+1each.10isolated.05s tip losses, allch3, rawtipstatus2
(MISSING), base1; masks4098tip/4096base. Volume/partial/brokenflags0, system0,
freshpairedframes, no overrun atloss. Replies87->55bytes, coherent32byte tip
transform omission. CRCstillunchecked. Lastvalidphi.648-.884rad,4loading/
6unloadingevents. Underlying field/orientation/cable/motion cause not unique.
Sevenoverruns/five>.10s gapselsewhere;12theta-unobservablegoodorientation frames
are separate. Keep invalid/gap guards, no threshold relaxation or whole-run
discard required solelyfor10gaps. Forwardphi peaks2.5268all3consistentwithbounds
fix; otherinlinecontroller installations stillnotestablished. Report/CSV/
FIG/PDF/PNG/source: residual_dynamics_learning/audits/validity_sine4_20261102_v1.
No hardware/model run, fitting or runtimechanges by assistant; logging status
pending below is superseded by this actual user-recorded observation.

2026-11-02 NDI diagnostics implemented after explicit user request. Driver/port
cache retains same-reply raw status/flags and availability, no extra BX query.
BaseRotation appendsport23, Position_only port20: ndiDiagnostics32x2 double
(base,tip). Existing ports/pose holds/validity thresholds unchanged. Helper/
schema/decoder and named optional CSV export added; CRC received, not verified.
32 synthetic reader/packet cases, metadata, original-output equality (except
host clocks), incomplete/recovery and no-extra-query pass. Actual legacy
SineWave3 observations unchanged. No serial/DAQ open, reader setup/preflight
or Simulink/hardware execution. Operator must refreshclass/ports with model
stopped and log newlastport asndiDiagnostics; physical integration/timing
unverified. See ndi_diagnostics/README.md, validation_20261102_v1 and
project_docs/NDI_DIAGNOSTICS_RESULTS.md. Actual dropout cause awaits new log.
No commit/push. Supersedes earlier diagnostics-proposal status below.

2026-11-02 SineWave3 quality follow-up:14 isolated acquisition losses/12721sensor
samples (13channel3,1channel1), same tip-invalid/base-valid/nonfiniteXYZ chain,
fresh paired frames/overrun0, one.05s interval each. Neighbor poses spanphi.186-
2.202rad; several recurring regions rather than one unique location.13ch3 losses
versus19inSineWave2 is not a controlled generator-move effect. Detailed raw NDI
status not logged, exact reason unresolved. User asks diagnostic-output need:
recommend same-BX raw sensor/port status, decoded volume/broken bits, quaternion
norm and host rejection mask, associated with existing frames/times. Proposed
only: no runtime reader/driver/SLX changes or hardware run. Current code discards
packet/system details; frame pairing alone does not verify packet integrity.
Report/CSV/FIG/PDF/PNG/source under residual_dynamics_learning/audits/validity_sine3_20261102_v1.

Standing user reminder2026-11-02: carry the equal three-actuator bound correction
into future work on other PAC and PDI/FB controllers (including PDFB/feedback-
linearization variants when they use the same q=[l2;l3],l1=-l2-l3 model).
Before porting/comparing, audit embedded old formulas and match forward/inverse/
controller force terms. +/-35mm is already the limit; include l1, do not just
increase limits. External-source correction does not mean all saved models are
migrated. D-12 in DECISIONS.md and OI-10 in OPEN_ISSUES.md retain the rule and
pending installation/build/validation work. Preserve declared legacy model
versions/results. This reminder authorizes documentation, not more SLX/hardware
work or unrelated controller edits.

2026-11-02 authorized correction: bound already+-35mm; now shared helper
armS_single_boundForce maps allthree actuator forces includingl1=-l2-l3 intoq.
Forward RHS, fixed/adaptive PAC and common-pressure inverse source use corrected
force. Paste-ready RL_data_collection/rlInverseModelBlock.m keeps3inputs/5outputs;
collection SLX still has old embedded inverse until operator pastes/rebuilds.
No SLX edited or hardware/Simulink run. Bound strength/mu/gains/pressure settings
unchanged. Offline mapping/potential/symmetry/forward-inverse/PAC tests passed.
3bar synthetic equilibria now2.5268rad inallchannels vslegacy3.0051/2.5268/2.5268.
Standalone Coder post-emission toolchain check failed; actual build unverified.
Report, transcripts, CSV/FIG/PDF/PNG and source hashes under residual_dynamics_learning/
validations/actuator_bounds_20261102_v1. Installation/version guidance in
RL_data_collection/ACTUATOR_BOUND_CORRECTION.md. Current RHS defaultboundMode1;
legacy datasets/frozen contexts/calibration remainmode0 explicitly; no retrain,
relabel or old-output overwrite. New residual baseline must be versioned.

2026-11-02 nominal bound asymmetry identified from fresh sine diagnostics:
des_phi peaks3.005/2.527/2.527rad are simulated, not measuredphi. Current RHS
penalizesl2/l3 only, omits reconstructedl1=-l2-l3. Channel1l1 reaches39.067mm
without boundforce; channels2/3 extendingcoordinate32.849mm sees~6.02N softbound.
No channel-specific phi clip. Possible symmetric actuator-force mapping logged,
not implemented; any correction needs versioned baseline/revalidation. Report,
CSV/FIG/PDF/PNG/source under audits/bound_asymmetry_20261102_v1. No new simulation,
training/hardware/helper/SLX edits. Operator repeating after generator movedback;
treat as setup change and recalibrate/preflight; new run not yet reviewed.

2026-11-02 fresh sine quality audit: reviewed20261102/SineWave.mat andSineWave2.mat
for configurationValid dropouts only; no fitting/prediction scoring.15/20 losses
out of12721sensor samples each, all isolated.05s simulation intervals, all tip
sensor invalid/base valid with nonfiniteXYZ/orientation invalid, paired fresh
frames andoverrun0.14/15 and19/20 during channel3; one earlier channel1 loss each.
Current XYZ-IK helper correctly rejects nonfinite input; detailed tracker status
not logged, exact physical cause unresolved. Preserve flags/mask invalid windows;
no validity/solver changes needed from this evidence. Report/events/nativechecks,
hashes andFIG/PDF/PNG in audits/validity_20261102_v1 under residual_dynamics_learning;
figures indexed in FIGURE_INDEX.md. Other20261102 recordings not reviewed here.

2026-11-02 fresh-test proposal: current frozen next-position residual predictors
can be tested offline on new operator-collected recordings. Recommend initial
StepWave/SineWave repeat, unchanged established collection protocol, complete out
and physical/simulation timing/validity/calibration metadata. Freeze weights,
training transforms and baselines; evaluate one-step and longer autonomous
physical-time rollouts before any retraining. Fresh-cohort evaluator adapter
still needed: current preparation manifest and rollout cohort are hardcoded.
No current-model PAC/PDFB integration or learned-pressure arm test proposed.
This is a proposed data test, not an accepted/started hardware schedule. See
RESIDUAL_DYNAMICS_RESULTS.md for limits and later expanded waveform coverage.

2026-11-02 discussion/priority update: user correctly noted the nominal model's
inaccuracy was already known. Calibration was a limited offline diagnostic to
test simple parameter/lag explanations, not a prerequisite to make the nominal
model accurate. Stop extending calibration as the main workstream. Next keep
an explicit frozen nominal baseline and learn missing generalized-force dynamics
from motion/pressure history, preferably fitting integrated length output error
rather than treating differentiated acceleration as measured force. Placement,
architecture and training method still need implementation/verification; earlier
step3 position-transition predictors are not this force-residual model. Compare
nominal+residual with nominal and simple/data-driven baselines on held-out longer
rollouts/loading-unloading before PDFB/PAC integration. No new training here.

User requests retained paper-useful figures. Index project_docs/FIGURE_INDEX.md;
versioned bundle residual_dynamics_learning/figures/20261102_model_diagnostics_v1
retains calibration/screening diagnostics, source tables/code/reports, PNG/PDF/FIG
exports and original example PNG with provenance/captions/claim limitations.
Future important analysis should retain plots and numerical inputs in this form.
Local logs/artifacts are not yet committed; commit/push only when requested.
The older calibration entry's next-step priority below is superseded by this one.

2026-11-02: executed existing-data static/dynamic identification offline. Canonical
report residual_dynamics_learning/calibration/greybox_20261102_v2/REPORT.md.
October19 only:15 StepWave training holds fit structured positive K and bias;
10 training windows fit D, then shared effective lag at fixed delay grid.
30 disjoint-trial validation windows, full nonlinear ~2s propagation: original
mean vector pathRMSE5.136mm -> mechanics3.369 -> mechanics+memory3.295.
Matched persistence1.338 and constantvelocity1.055mm remain better. Static
validation4holds improves7.604->5.583mm, still substantial equilibrium mismatch.
Both fitted D reach250 upperbound; training-only profile favors even higher D,
consistent with suppressing model motion rather than identified physical damping.
Delay/T trade off; selected.2s grid edge/T.01066s is not identified actuator lag.
Pressure commands only, approximate physical clock mapping; fixed force gain,
mass and geometry. Oct27 not loaded here; prior screening already evaluated it.
Simulator original-RHS/leakage/delay tests and disjoint-trial checks pass.
No SLX/controller/hardware/commit/push changes. Next inspect equilibrium error
versus common/differential pressure/loading and timing before more capacity or
arm integration; possibly design informative operator-run loading/unloading data.
Do not use fitted values as PAC tuning. Incomplete calibrationv1 quarantined.

2026-11-01 step3: trained separate forward/inverse ridge, MLP32/16 and GRU16
onOct19 train anchors(stride3=11,785), validation13,977; fixed small20epochCPU
budget/seed73. Model transforms/selection frozen before external evaluation.
Canonical results residual_dynamics_learning/benchmarks/screen_20261101_v2;
v1 superseded numerical ridge warning, fixed double-QR refit reproduces exactly.
Forward macro validation RMSEmm nominal4.9845,CV.06128,ridge.05802,MLP.41469,
GRU.30156. External CV.05193,ridge.08939,MLP2.4797,GRU2.8055; kinematic ridge
ablation.04887 -> nominal prior benefit not established. Short autonomous
10step validation rollouts ridgefinal~.33mm versusCV~.37, not full-circle proof.
Inverse external pressure RMSEbar ridge.000626,MLP.02779,GRU.02764, but pressure-
history/common-pressure-only AR.000424 and no-pressure-history ridge.40082.
Thus tiny inverse errors reflect a strong command-predictability shortcut,
not demonstrated missing mechanics or new-motion inverse capability. Bounds
audited unclipped. No candidate deployed, SLX/hardware/commits unchanged.
Next investigate prior calibration/actuator-memory and motion-conditioned
identification/excitation before arm tests. Oct27 now evaluated; not untouched
confirmation for future iterations. Full protocol/weights/tests/report saved.

2026-11-01 residual step2 completed offline: valid package is
residual_dynamics_learning/datasets/two_day_v2. Ten raw hashes checked,
day-separated exports, physical/simulation clocks, strict invalid/gap/trial
segments, causal derivatives/warmup, anchor indices and train-only normalizers.
73,406 shared forward/inverse candidate anchors:35,350train +13,977validation
fromOct19 whole-trial split;24,079 external-test anchors fromOct27, no fit there.
residualBatch enforces forward causality and inverse pressure-history k-1;
inverse desired-preview/common-pressure conditional semantics documented.
Nominal helper matches existing RHS and feasible inverse roundtrip. Candidate
residual labels use one explicit zero-activation reference and deriveda, not
measuredforce. Recorded maps revealOct19 step/sine/combined inverse activation
.8 versus0 elsewhere; not evidence of physical deadzone change.
Day2 step/combined physical dt~.0885/.0877s and many>.10s gaps leave46/26 windows,
so not representative external tests. Hold thresholds unchanged; no architecture
winner or trained model yet. validateResidualNominal/Data/recorded-map checks
pass. Initialtwo_day_v1 invalid/incomplete, explicitly quarantined DO_NOT_USE
after helper sign error caught/corrected. No raw changes, hardware, SLX edits
or commits. Step3: baseline/memory/MLP/GRU preliminary offline benchmarking.

2026-11-01 residual step1: logged user's six-step plan and completed interface/
target specification. Audited ten raw RL recordings separately: root2026-10-19
and20261027 subfolder2026-10-27, dates confirmed by saved acquisition metadata,
five waveforms each. Different model names/layouts and durations; no pooling.
Thirty final NI command channels identified with indices/paths, ten source
hashes, current saved SLX source snapshots extracted read-only. Specification
residual_dynamics_learning/MODEL_SPECIFICATION.md defines forward2-output
acceleration residual and inverse3-pressure interface with explicit common
allocation (recommended2 differential residuals + prescribed common pressure).
Nominal inverse minimum-offset allocation differs from arbitrary recorded
common mode. Prevent current-pressure/future-measurement feature leakage;
inverse desired preview permitted, forward future observations are targets only.
Shortlist memory baseline/windowedMLP/GRU, LSTM conditional; no winner selected
without training. Step2 data alignment/quality/splits and nominal verification
remain next. No raw changes, model training, SLX edits, hardware runs or commits.

2026-10-28 branch transition: user requested relevant checkpoint and new branch.
Committed PAC code, models, campaigns/raw data/reviews and documentation as
5532a5d on codex/pac-autotuning-framework. Now codex/residual-dynamics-learning.
Research objective: learn missing dynamics as MODEL augmentation, validate
offline predictions, then compare original/augmented PDFB and PAC with manual
experiments. Not a separate residual feedback controller. Continuous online
tuning deferred. See RESIDUAL_DYNAMICS_LEARNING_PLAN.md for first audit, target
identifiability, train/test protocol and comparison plan. No residual model
fitted/deployed, hardware run or remote push. Unrelated RL changes (including
staged rename), manuscript files and backups remain outside checkpoint.

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
