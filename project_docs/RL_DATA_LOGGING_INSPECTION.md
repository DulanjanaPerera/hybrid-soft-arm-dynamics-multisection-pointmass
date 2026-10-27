# Open-loop RL/system-identification logging inspection

Date: 2026-10-19
Model: Dynamic_model_Val_FWDnINV_v2.slx (saved file)
Scope: read-only archive and MATLAB load_system/parameter inspection. No update,
compile, simulation, hardware access or SLX modification. Open-loop/manual data;
controller mode toggles are outside the requested inspection.

## Confirmed saved logging

SignalLogging=on, logsout Dataset, ReturnWorkspaceOutputs=on, output variable out.
SaveTime=on (tout); no global point limit, decimation1. Active logged output
ports individually have DataLoggingDecimateData=off and DataLoggingLimitDataPoints=off.
Displayed MaxPoints5000 is inactive, NOT a truncation cap. StopTime30s.
LoggingToFile=off: workspace logging does not automatically persist a trial MAT.

- Real arm thetaOrientation/phiOrientation and orientationValid/thetaObservable.
- Paired sensor validity, frame, hostTime/pollStart, overrun.
- R_NDI_from_arm, gravity_arm, baseOreintationValid and nonBendingRotation_rad.
- Raw-position/error outputs raw1..raw8; sensor XYZ/fitResidual/configurationValid.
- Orientation-inferred three lengths and filtered two-coordinate l/dl/ddl.
- Requested manual pressure vector from root Mux.
- Final post-Kill/post-saturation P1_NI_DAQ/P2_NI_DAQ/P3_NI_DAQ before voltage gains.
- Kill switch; command-voltage gain outputs logged but unnamed.
- Simulation clock, sampled clock, derived measurementTime_s.

## Important interpretation and gaps

1. Duplicate pressureCommand_bar: root Mux is actual manual request; root MATLAB
   Function3 is inverse-model prediction and is NOT connected to valve command.
   Never train by first-name-match selection; resolve BlockPath or rename later.
   Model thetaObservable also duplicates measured sensor thetaObservable.
   length, frame, overrun and X have duplicate names too; some are redundant,
   some physically different. Preserve block provenance in export.
2. des_theta/des_phi come from forward dynamics driven by requested pressures,
   not an independently commanded physical target. Treat as predicted theta/phi.
   For open-loop identification, measured response+applied pressure is primary;
   do not mistake prediction as observed state or desired policy action target.
3. No pressure sensor inputs. NI_DAQ signals are commanded pressure setpoints in
   bar, not measured chamber/regulator pressure. System identification can fit
   command-to-motion dynamics, not separately identify pressure response without
   additional measurements/assumptions. Gain10/9 maps bar to command volts.
   Channel mapping P1/P2/P3 -> NI channels1/2/3. Saved saturation0..3bar;
   this is model configuration, not a newly authorized experiment range.
4. NDI SamplePeriod=.05s (20Hz), reference hold=.05s; NI output=.01s (100Hz).
   Other-chat statement50Hz is not this model's configured NDI rate. Compiled
   times/rate compliance were not tested. Logging repeated values faster does
   not create new physical observations. Keep frame/freshness/timestamps.
5. hostTime and pollStart are monotonic host elapsed-time estimates for sensor
   transactions; measurementTime_s is generated from sampled simulation time
   and validity, not an independent hardware acquisition timestamp. Alignment
   of NI writes and NDI wall time must be checked in a short saved trial.
6. Raw1..8 are XYZ plus tracking error for two sensors (XYZ converted/masked by
   reader); raw BX quaternions and complete raw pose packets are not exposed/
   logged. Reader consumes quaternions internally. Existing calibrated angles
   and R_NDI_from_arm suffice for current features; raw base/tip quaternions are
   recommended for future orientation recalibration/reprocessing.
7. Sensor position is rotated relative tip-minus-base XYZ (mm->m), with a later
   hold-last-finite block feeding meas_position. Use configuration/orientation
   validity to reject held/repeated values. Full arm-frame rotation alignment
   for task-space labels should be documented; no independent strain sensing.
8. Persist session/trial ID, date, module/port mapping, calibration file,
   pressure units/conversion, sensor mounting, base/payload conditions and manual
   excitation phase/reset notes alongside MAT. No dedicated trial manifest or
   phase-ID source was found in this model. No controller-mode input requested.

## Proposed next discussion

Model is broadly instrumented for command-to-motion acquisition, but names,
raw-orientation retention and trial metadata need agreement before a large batch.
A short saved raw trial is required to verify actual coverage, monotonic times,
fresh frames, channel mapping and wall-time synchronization. Model inspection
alone does not validate a training dataset. No logging edits performed.

## Offset-correction helper added
moduleTipXYZBlock.m consumes existing Sensor/position (after rotateX90), calibrated R_NDI_from_arm, measured orientation theta/phi and current orientation/base validity. Outputs moduleTipXYZ_m and tipXYZValid in the current module-base frame. Keep original sensor XYZ too. Assumes axial offsets .055/.050 m; lateral mounting errors are not calibrated by this helper. No model wiring changed.

## 2026-10-19 - First RL open-loop30s trial validated
forward_test_RL_30s.mat read offline. Requested/NI-pathP1identical,exactwaveform
1.5*(1-cos(.2*t))bar(range0-3),P2=P3=0. Period31.416simulationseconds,
30scontains.955cycle. pressureDesired_bar nowdistinctinverseprediction;
pressureUsed_bar max2.2 ismodelmax(p-.8,0),NOTphysicalcommandclip.
601NDIsamples20Hz,600valid,invalid8.30s;pairedframes601,norepeats,basevalidall.
10overruns;hostelapsed31.484536vs30simulationseconds,hostdtmedian.05629,
min.021985max.347293. Needhost/simulationmappingandrecordedcommandtimes;
NIwritewalltimestampsnotpresent. Do notassumeuniformphysical50msormarkall
signalswithsamewallclock. Onlyphysicallograte20Hzdespite100Hzrepeatedlogs.
moduleTipXYZValid600;recomputedhelpermatchesloggedXYZexactlyonvalidsamples,
invalidXYZNaN. moduleendpointbounds[-.12238,-.001387,.09710]to
[-.002866,.014492,.17554]m;phiOrientationmax1.88401rad. Rawpositions/errors,
Rbase,gravity,flags,host/poll/frame,commands/volts present. Rawquaternionsabsent.
ConstantP2/P3/killlogonly2points0..18.86;notmissingacquisitionpersebutexportmust
handleconstantlifetimesexplicitlyorpersistmetadata. Othermodel/sensorname
collisionsandunnamedsignalsremain;exportbyprovenance. Datasetcontainsoutonly,
noembeddedtrialmanifest. Goodplumbingtrial,notcompleteRLtrainingbatch. No model
orsourcechanges, nohardwareexecution.

## 2026-10-19 - RL_dataCollection circle feedforward wiring audit
Read-only archiveinspection RL_dataCollection.slx saved14:46. circlesweepout1/2
->Mux3theta/phi->MATLABFunction config2length(reduced2coords)->DesiredlengthsZOH.05
->MATLABFunction2lengthDerivatives(withsampledclock.05)->Mux6[l;dl;ddl]
->MATLABFunction3inverse->pressureDesired_bar->PneumaticValves; correctfeedforward
statepath,notmeasuredlengthfeedback. Inversedeadzone0,max3,K650,D40,gravityfixed
[0;0;-9.81],stiffinputunused. CRITICAL From3 tagA1 feedsinverseorientationValid;
GotoA1 carriesmeas_theta,notorientationValid. ReconnectexistingNDIout4(current
orientationValid_control)orproperdesiredvalidgate; theta0wouldfalsegatepressure.
StopTime342truncatescircles620;set625forcircleor530forhelix. Helixblockpresent
butnotconnected;circlecurrentlyactive. Forwardpredictionstill subtracts.8,so
notanexactinverse-forwardvalidationpair;doesnotalteractualinversecommand.
DesiredratesfilteredtauV.10/tauA.15,notanalytic,anglevelocityswitchesatcircle
start/endcanproduceboundedderivativetransients. Desiredreferenceangles logged
as theta_from_block/phi_from_block; des_theta/des_phi remainforwardpredictions.
Reusedlength/d_length/dd_lengthnamescollideforwardanddesired. Requestrename
inverseoutputpressureRequest_barforclearactionprovenanceifuserchooses;NIcommands
remaindatasetaction. No SLX/sourceedited,compileorhardware execution.

## 2026-10-19 - Saved RL tag wiring correction verified
RL_dataCollection.slx nowGotoA1fedbyNDIout4orientationValid_control;From3A1feeds
inversevalidity. MeasuredthetaNDIout1usesnewGoto/Frommeastheta. Correct separation.
Savedlogginginstrumentationretainsmeasuredtheta/phi/inferredlengths,sensorXYZ,
moduleTipXYZandvalid,finalNIpressures,timestamps/frame/validity. StopTimestill342;
completecircle625orhelix530required. No model/sourceedits duringinspection.
Userfocusphysicaltheta/phi/lengths/XYZ/pressure;retainingvalidity/timestampsneeded
tointerpretsamples,notadditionalcontrollertrainingfeatures.
