# Cygnus 6D/Aurora raw log: tip MISSING during the P3 sequence

2026-11-09. Source C:/Users/dperera/Desktop/log.txt,12416162bytes, SHA256
56204593a7734975c1db959d60607583936e5ea6b81f0c47891e328f5c073b60.
The original source is retained as source_log.txt. This is a device-protocol
log, not a MATLAB signal log; user confirms no matching MAT file exists.
User then confirms the standard sinewave block used for collection/validation.
No tracker/hardware/model access or source modification was performed.

## Directly verified findings

Decoded10369 BX:0001 replies. Header and payload CRC16 pass for every reply;
all lengths, tool counts and payload boundaries parse consistently. Separate
CR raw-sensor replies (0xA7C4 packets) were not misread as BX transforms.
The tracker reports handle-status0x01valid/0x02missing/0x04disabled;0x02means
the device cannot supply a transform, not merely unobservable azimuth phi/theta.
Primary interface definitions/documentation:
https://github.com/PlusToolkit/ndicapi/blob/master/ndicapi.h
(NDI_HANDLE_VALID/MISSING/DISABLED and ndiGetBXTransform/Frame documentation).

There are two BX-bearing sessions:2setup replies in session2,10367replies
in final session3.8initial final-session reports precede the first jointly
valid base/tip reading. Configuration/startup missing/disabled tools are
reported separately and excluded from runtime counts.

| After initial acquisition | 0A(base) | 0B(tip) |
|---|---:|---:|
|Runtime reports|10359|10359|
|MISSING transforms|0|23|

Handle assignment follows the existing project acquisition settings
ndiSupportFnc.acquisitionSettings:0A=base,0B=tip. Both are6D tools according
to the log initialization. All23runtime tip misses have the base valid and
matching base/tip frame numbers.23reports comprise17separate consecutive
missing episodes, maximum2reported updates each. At the reported40.033Hz
transform update rate these span approximately25-50ms to the next valid
reported transform. This is observed report duration, not proof no additional
unqueried frames were missing.

Runtime tip-missing fraction.2220% of all logged runtime reports. Misses occur
155.7215-209.3023elapsed seconds. Trackerframe49926first missing,67086last.
Full per-report list with source line numbers: missing_tip_reports.csv;
episode grouping/bounds: missing_episodes.csv. Original line55122 is the first
runtime tip-MISSING BX reply, with a valid0Atransform and0Bstatus0x02.

Both out-of-volume and partial-volume bits are clear at all23misses; system
status is0for every final-session BX packet. These flags do not identify a
unique physical cause and do not rule out pose/field/cable/interference issues.
Checksum success rules out CRC corruption in the recorded BX packets, not
all possible serial transport/acquisition faults elsewhere.

## Time contract and query gaps

The log reports Track Frequency40.033Hz and Frame Frequency320.264Hz.
These are different: frame counter normally advances8ticks per returned
transform update. Elapsed time is(frame-frameFirstFinalReply)/320.264.
Origin is the first BX reply in final session, not MATLAB time or wall-clock
start of the sine experiment. Total final-session elapsed260.560s.

There are43reply intervals longer than8frame ticks; largest88ticks(~.275s).
They correspond to65unqueried nominal update slots. This is separate from
explicit MISSING status: absent queries cannot be classified as good or bad
transforms. No repeated or mismatched paired frame counters appear here.

## Association with P3

The standard source RL_data_collection/rlPressureSineBlock.m is P1thenP2thenP3,
3cycles per trial,2repeats per channel,omega.8rad/s,peak3bar,5sinitial zero,
10srelaxation. Its6three-cycle motion groups are visible in the raw poses.

There are no actual pressure values, channel markers or synchronized command
timestamps in this text log. We therefore aligned this **known timing template**
to raw tip-to-base orientation change. A common offset was fitted with
separate response amplitudes for each channel; bestoffset14.625s,orientation
fitRMSE.134rad. This offset includes physical actuator response lag and is
not an independently measured clock synchronization.

| Inferred P3 excitation window(elapsed seconds) | Tip MISSING reports |
|---|---:|
|153.873-177.435(first repetition)|12|
|187.435-210.997(second repetition)|11|

All23runtime misses fall within these P3windows. None fall in inferredP1/P2
excitation. Channel association is strongly supported by the confirmed sequence
and six-group motion pattern, but remains an inference without pressure logs.
Do not present the plotted pressure template as an observed pressure trace.

The total relative-orientation plot is2acos(abs(dot(qRelative,qReference))),
using a1-5sstraight reference. It is NOT the calibrated continuum-arm phi:
mounting offsets, twist separation and project sensor compensation are not
reconstructed here. Missing transforms remainNaN; no interpolation was used.
Raw XYZ is sensor position in the tracker frame, not offset-corrected arm-tip
position. Adjacent valid samples around the17episodes have relative rotation
change approximately.527-1.047rad(median.847), below this motion's~2.4radpeak.
Thus losses recur while passing particular intermediate poses on loading/
unloading, not exclusively at the maximum bend. This does not establish
whether motion speed, pose/field location, wiring or another effect is causal.

## Implication for the controller project

Cygnus independently receives missing tip transforms during the standard
pressure excitation. The issue is not solely created by MATLAB's orientation
validity logic, length bounds or the residual-selection Switch. A missing
transform withholds both raw position and orientation; smooth held XYZ output
elsewhere is not evidence of fresh position tracking during those frames.

The earlier residual-on experiment separately showed that brief invalidations
trigger1.6sobserver re-warm and abrupt fallback-pressure dips. That explains
how a~25-50ms tracker loss can become a much longer physical tracking dip.
Command-handoff smoothing addresses that consequence; it does not repair
the original tracker loss or justify making invalid transforms valid.

## Retained files/reproduction

-summary.json,sine_alignment.json,tracking.csv,missing_tip_reports.csv,
 missing_episodes.csv,motion_alignment.csv.
-tracking_and_sine_sequence.PNG/PDF/FIG,source_log.txt.
-source_snapshot,protocol_reference.json,input/artifact SHA256manifests.
-parser_checks.json: CRCknown-vector, corrupted/truncated-packet rejection.

```
python RL_data_collection/analyzeCygnusBxLog.py SOURCE_LOG FRESH_OUTPUT_FOLDER
```

Then offline MATLAB: plotCygnusSineAudit(outputFolder), with
RL_data_collection on the MATLAB path. No model load/run or hardware call.
No runtime/controller/NDI-source edits,training,commit orpush.
