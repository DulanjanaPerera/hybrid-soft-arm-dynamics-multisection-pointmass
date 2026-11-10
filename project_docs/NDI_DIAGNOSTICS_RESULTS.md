# NDI diagnostic-output implementation

## 2026-11-10 - Literature context for in-volume missing transforms

Onlineprimarysourcesfound: Mandal2016thesisp84Fig6.2directlydocuments
AuroraBadFitatanorientationwhiletipinsideworkingvolume; Scott/Liu/Wu2024
continuumrobotpaperp5briefBadFitdatainterruptions. OlderNDIauthoredV2guide
mirroredonline discussesmissingbeforevolumeedges/interference/motion/cables.
Geometrychangesinpapersareobservations/attempts,notuniversalfixes orproof
thatcloserFGalwaysworse. Officialbrochureaccuracydependsenvironment/tool.
Exactlinks/sectionsandversionlimits:NDI_TRACKING_LOSS_LITERATURE.md.
Ourarmmovesdistanceandorientationtogether;causalitynotseparated. Nooutvolume,
fixedangle/clearancethreshold,damageorfirmwarediagnosisestablished. No runtime/
NDI/hardware/modelchange,deviceaccess,commit,push.

## 2026-11-09 - Arm rotation transfers tracking loss from P3 to P2

Usernominal180degarmrotation,P1now+X,explicitlyconfirmedP3->P2->P1order.
NewDesktoplog6335092bytes/SHA0f9f5079...,5279BXallCRC/layoutpass;final
5277replies132.216s,6startupthen5271runtime.0Abaseallvalid,0Btip36missing
in26episodes~25-50ms,all66.795-91.749s inP2secondgroup. P3/P1zero.
No volume/partialbits,systemstatus0. BeforemissingP3/P2/P1=32/0/0,
after=0/36/0. CurrentP1rawpositiveXmotiontracks;lossnotfixedP3valveor
generic+Xeffect. Bothfailinggroupsbend towardlessnegativeNDIZ;thisis
association,notuniqueZcause,outvolumeclaim,orcontrolledmatchedtrajectory.
Basequaterniongeodesicchange171.17deg supportslargeorientationchange.

CombinedmanualP3-no-actuationlossesruleoutvalveswitchingasnecessarytrigger;
rotationresultsfavortracker-relativepose/orientation/externalgeometryover
PMAidentity. Cable/field/motionfactorsunresolved. Keepobservervaliditygates.
IndependentcontrollerSwitchNaNcommandproblemremainsseparate.

Reusableoffline classifyCygnusStepGroups.py andcanonicalraw/events/PMAcounts/
before-aftertable/PNG/PDF/FIG/hashes:
residual_controller/Experiments/Model_Val/cygnus_rotated180_steps_20261109_v1.
NoruntimeNDI/controller/SLXedit,device/model/hardware/trainingrun,commit,push.

## 2026-11-09 - Pressure steps P3first/P2/P1: losses remain P3-specific

NewDesktoplog5986190bytes/SHA1976de55...,userreportsmanualindividualpressure
stepsP3->P2->P1withlongpauses.4978BXallCRC/lengthpass;final4976/125.496s,
3startupthen4973runtime. Base0Amissing0,tip0Bmissing32/30episodes(~25-50ms),
allbasevalid/framepaired,volume/partialbits0,systemstatus0. Explicittracker
MISSINGnotjustquerygaps. ObservedresponsewindowsP3[8.97,37.32],P2[57.10,
83.81],P1[96.25,121.22]s. P3all32losses9.09-37.24s;P2/P1zero. Group
labelsuserorderplusmotion/restsegmentation,notloggedpressuretimestamps.

P3lossesontransitions,stepplateausvalid;doesnotestablishintermediateposes
trackwhenheld. P2/P1rawspeedindicatorscanoverlapP3missingintervals,not
universalvelocitylimitidentification. P3firststillfails,soeffectdoesnot
merelyfollowlastchannelorlongduration. PriorP3manualwithoutvalvesalsolost
tip:actuationnotnecessary,physicaldirection/pose/motion/cablecauseunresolved.
Nooutvolumeassertion. No observationgateweakening/interpolationproposed.

Canonicalreport/raw/events/groupcounts/PNG/PDF/FIG/hashes:
residual_controller/Experiments/Model_Val/cygnus_pressure_steps_P3_P2_P1_20261109_v1.
Genericofflineplotternowoptionallyannotatesmovement_groups.csv. No runtime
NDI/controller/SLXchange,device/model/hardware/trainingrun,commit,push.

## 2026-11-09 - Manual P3-direction/no-actuation test reproduces tip MISSING

User'snewDesktoplog3640632bytes/SHAa6a8cd6c...: manualP3-directionmovement
withoutpneumaticvalveactuation.3035BXallCRC/lengthpass;final3033/76.037s,
2startupthen3031runtime.0Abaseallvalid;0Btip30MISSING(.990%)in28episodes,
26one-update/2two-update~25-50ms. Allmissesbasevalid/pairedframe;noout/partial
volumeflags,systemstatus0. Noentire-systemserialCRCfailureidentified.
Misses30.55-69.92sfromfirstfinalBX. Reportedframe320.264Hzvsupdate40.033Hz
handledcorrectly;querygapsdistinguishedfromexplicitmissingstatus.

ValveactuationisNOTnecessaryunderuser-reportedtestconditions. Thisrulesout
valveswitchingasthesolecause,notallpossibleadditionalactuationinterference.
Direction/pose/orientation,cablemovement,motionqualityremainunresolved;
notproofoutofvolumeoruniversalspeedlimit. Someadjacentvalidbase-frameposes
overlapearlierP3loss(min1.12mmnearestneighbor),manualcoverswiderregion.
Previousgenericmanual0lossesandnewP3manual30lossesaredifferentcaptures.
Potentialnextmanualtest:holdanimplicatedposevscontinuousmotion;notrunhere.

Canonicalraw/events/plots/CSV/posecomparison/hashreport:
residual_controller/Experiments/Model_Val/cygnus_manual_P3_20261109_v1.
Genericmanualplotternowderivescountsfromcapture. NoruntimeNDI/controller/
SLXedits,device/model/hardware/trainingrun,commit,pushduringaudit.

## 2026-11-09 - Fast manual/no-pressure capture: no runtime MISSING

NewDesktoplog.txt2973106bytes/SHA7253f838... differsfrompriorP3sinecapture
archivedseparately. Usermanualrandomfastmotion,noactuation.2448BXallCRCs/
lengthspass. Finalsession2446reports/61.449s;2initialbothmissing,2444later
bothvalid. Runtimebase/tipMISSING0/0,duplicates0,exactconsecutivetransforms0,
pairedframesallmatch,system/volume/partialflags0. No evidenceoutofvolume;
userstatesinsidevolume. Initialsetupdisabled/missingexcludedfromruntime.

Tworuntime50msqueryintervals,tworuntime175msintervals;longonesarefinaltwo
replyintervals61.10-61.45s immediatelybeforestop. Cannotcallthemexplicit
sensorMISSINGorprovenmidrunfreeze. App/serial/stopschedulingcauseunknown;
unqueriedintermediateframevalidityunknown. Mostqueries25ms,allreturned
runtimeposevaluesfresh. GUImeasurementfreezebehaviornotverifiedbythislog.

EarlierP3sine23explicitmissesremainvalidfinding;newmanual0doesnotreproduce.
Rawposecoverage differs(tipZmanualmax-222.6mm,adjacentprior-missZ-220..-169mm,
baseframepositionnearlysame). This isnotoutofvolumeclaim,pressurecausalproof,
orspeedlimitidentification. Furtherphysicalcausalityremainsunresolved.

ParseranalyzeCygnusBxLog.py nowoutputsvalidempty-eventfilesforzeromissings;
plotCygnusManualAudit.madded,nostandardsineassumption. Canonicalmanualreport/
rawlog/CSV/gaps/flags/PNG/PDF/FIG/hashes:
residual_controller/Experiments/Model_Val/cygnus_manual_20261109_v1.
No device/runtime/SLXchange,model/hardware/trainingrun,commit,push.

## 2026-11-09 - Cygnus6D rawBX log: missing tip in inferred P3 phase

UserprovidedDesktoplog.txt(12.416MB),confirmsstandardcollection sine block,
noMATcommands. Offlinehexdecoder10369BX:0001repliesallheader/payloadCRCpass,
recordlengths/handlestatusboundariesvalidated. CygnaCRrawsensorpacketsignored,
notmistakenfortransforms. PrimarydefinitionsPlusToolkit/ndicapi.h handle
status1valid2missing4disabledandBXtransformlayout. No directtrackeraccess.

Finalsession10367queries/260.560s;8initialreportsbeforefirstbothvalidexcluded.
10359runtime:0Abasemissing0,0Btipmissing23(.222%),17episodes,maximum2updates.
Reportedtransformrate40.033Hz,framecounter320.264Hz(8ticks/update). Individual
misses~25-50msto nextvalidreportedtransform. All23basevalid/pairedframes,
systemstatus0,volume/partialbitsclear. MISSINGisexplicittrackerresponse,not
MATLABthetaobservabilityorlengthbound. No uniquephysicalcausefromflagsalone.
43queryintervals>8ticks(max88/.275s) aredistinctfrommissingstatus;do notcount
unqueriedupdatesasknowntrackingfailures.

Misses155.722-209.302s fromfirstfinalBX. Standard6groupsP1,P1,P2,P2,P3,P3,
omega.8,3cycles,10srelax;motion-fitcommontimeoffset14.625s. InferredP3windows
153.873-177.435(first12misses),187.435-210.997(second11):all23inP3noneP1/P2.
Thisisconfirmed-sequence/motionassociation,notactualsynchronizedpressuredata.
Relativeorientationchangeplottedisnotcalibratedphi;rawXYZnottipoffsetcorrected.
Adjacentvalidrotationatmisses.527-1.047radbelow~2.4peak,suggestsintermediate
posepassesratherthanexclusivepeakanglefailure. Pose/field/cable/motion/EMI
causalityunresolved;no arbitraryvaliditygatechangesproposed.

IndependentCygnuslossconfirmsroottrackingissueexistswithoutMATLABvalidity
processing. PriorPAC/desiredpressurefall-backcanamplify25mslossinto1.6swarmup
andphi dips;commandhandoffonlyaddressesconsequence. Missingtransformcontains
neitherfreshpositionnororientation;publicheldXYZmustnotbeconsideredfresh.

Newoffline RL_data_collection/analyzeCygnusBxLog.py,plotCygnusSineAudit.m.
Canonical residual_controller/Experiments/Model_Val/cygnus_sine_20261109_v1:
archivedsource_log,alltracking/missing/eventsCSV,alignmentJSON,PNG/PDF/FIG,
sourcehashes. ParserCRCknownvector/corrupt/truncatedrejectionchecks pass.
No NDI/runtime/SLXsourceedit,model/hardware/trainingrun,commit,push.

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

## 2026-11-02 - Motion speed is a documented possible MISSING cause

The NDI-authored Aurora V2 User Guide, Revision 3, printed page 52, explicitly
lists excessive tool motion among possible reasons for a MISSING transform,
alongside interference, system errors and leaving the detection region.
Source: https://www.scribd.com/document/1010309716/Aurora-V2-User-Guide
The reviewed sections do not specify a numerical translational/angular speed
ceiling. Their measurement-rate specification is not such a speed ceiling and
does not establish the installed unit's exact configuration.

This makes motion a plausible explanation for SineWave4, not a demonstrated
cause. Its pose/direction clustering does not separate speed from field,
orientation or cable effects. Proposed operator-run comparison: repeat the same
channel-3 excitation at twice the period with amplitude, generator and cable
placement fixed; compare missing rates per passage near matched measured poses,
then hold near the affected pose if feasible. Fewer losses at slower motion
would support motion dependence but would not uniquely identify the mechanism.
No hardware run, runtime setting or validity gate changed for this review.

## 2026-11-02 - SineWave4 confirms logged tracker MISSING transforms

Actual user recording contains schema1 matrix; all frames/flags/timestamps align
and reply counter increments one per sample.10/12721 isolated tip losses, all
channel3, lastvalidphi.648-.884rad. Every rawtipstatus2=MISSING, base1, masks
4098/4096; no volume/partial/broken/overrun flags at those events, system0.
Packet body87->55bytes matches omitted32byte transform. This establishes the
rejection stage, not a unique physical cause; CRC not checked. Additional7
overruns/5>.10s gaps elsewhere and12near-straight theta-unobservable samples are
separate. No new runtimecode/threshold/hardware execution. Source/event tables,
plots/verification in residual_dynamics_learning/audits/validity_sine4_20261102_v1.
Operator diagnostic logging now observed; broad hardware timing/other branches
not thereby qualified. Future short controlled pose/cable/field tests proposed,
not run. Retain masking; no need discard all data solelyfor10one-frame losses.

## 2026-11-02 - Same-BX diagnostics added; offline verification passed

User authorized updating NDI functions to expose rejection causes. No SLX edited,
serial/DAQ device opened, hardware model run or reader setup/preflight executed.

Driver/PortHandle retain current reply membership, raw sensor/port status,
availability and packet/system metadata. measureTipPositionAll optionally
returns an eighth diagnostic matrix; getLastDiagnostics reads cache only.
Omitted/disabled status words remain unknown instead of reusing old flags.
Existing acceptance, masking, byte-read order and exception policy are retained.

Both readers append32x2 double ndiDiagnostics, base/tip columns: BaseRotation
port23, Position_only port20. Finalizer records frame/quaternion/calibration
conditions and existing final validity results. Existing output order, types,
shapes, held estimates and threshold rules remain. No new control/pressure gate,
device query or per-step disk logging. CRC is received only; crc_checked=0.

ndiDiagnosticsSchema returns field/bit tables; decodeNdiDiagnostics returns
named sensor/condition tables. Optional exporter support produces64 named fields
and schema/reason CSVs without overriding prior validity masks. Old logs supported.

Verification: ndi_diagnostics/validation_20261102_v1.32 reader cases replay exact
saved/current step bodies against actual decoder with synthetic transport. Legacy
non-wall-clock outputs match for normal/recovery/missing/disabled/omitted/
volume/partial/broken/repeated/mismatched/quaternion/calibration/nonfinite/
straight/wrap cases. Metadata, one-poll/cache-no-query, incomplete exception/
recovery, named export and unknown-stage decoder pass. Actual legacy SineWave3
observation table remains identical. No real-time/hardware qualification claimed.

Operator: save current out, refresh class/ports with model stopped, log the new
last port as ndiDiagnostics, typically through a Terminator connection. See
ndi_diagnostics/README.md. Model compilation/installation and physical timing
remain pending. A new recording is needed; no reason codes recovered from old
Boolean logs. No raw/model-weight changes, commit or push.
