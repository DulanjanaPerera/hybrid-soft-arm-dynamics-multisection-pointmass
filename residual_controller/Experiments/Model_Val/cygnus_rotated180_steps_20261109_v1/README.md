# Arm nominally rotated180degrees: missing transforms move from P3 to P2

2026-11-09. New Desktop/log.txt6335092bytes,SHA256
0f9f50794d0c91262b9719f455f0260de4d74e7b3d1add897dfde06c37c3f974.
Header18:13:31CST; lastwritten18:16:43. User reports rotating the arm180degrees
and thatP1nowbends in+X. User explicitly confirms activation order remains
**P3 ->P2 ->P1**, with separate pressure-step groups. No pressure values or
synchronized channel markers are encoded in the tracker log.

## Direct result

All5279BXreplies pass header/payload CRC and layout checks. Final tracking
session5277replies/132.216s;6initial acquisition reports precede the first
jointly valid pair. Runtime5271:base0Avalid5271/missing0;tip0Bvalid5235/
missing36(.683%). All36tip misses have the base valid and matching frame
counters. No out-of-volume or partial-volume flags; systemstatus0throughout.

| PMA, confirmed order | Approximate response window(s) | Tip MISSING |
|---|---|---:|
|P3 first|20.68-54.26|0|
|P2 second|66.67-about91.8|36|
|P1 last|100.62-125.27|0|

All losses occur66.7949-91.7493s. They form26consecutive-missing episodes:
16single and10double reported updates(~25-50ms to the next valid report).
These are actual handle-status0x02missing transforms, not merely an unqueried
frame. TrackFrequency40.033Hz andFrameFrequency320.264Hz are distinguished;
time=(frame-firstFinalBXFrame)/320.264. No invalid transforms are interpolated.

## Comparison with the unrotated step experiment

Same user-confirmedP3/P2/P1sequence, previous report
../cygnus_pressure_steps_P3_P2_P1_20261109_v1/README.md:

| PMA | Before rotation | After nominal180degree rotation |
|---|---:|---:|
|P3|32missing|0|
|P2|0|36missing|
|P1|0|0|

The loss does not remain attached to theP3actuator/valve identity. Current
P1 has strong raw tracker+X displacement(median about+159mm over baseline)
and no loss, so a generic+X-only explanation is also insufficient.

Both failing groups bend toward less-negative raw trackerZ: previouslyP3,
currentlyP2. CurrentP2median raw displacement over initial baseline is
[-66.5,+133.4,+163.1]mm; currentP3is[-99.8,+195.0,-125.5]mm and tracks,
currentP1is[+159.3,+142.8,-11.1]mm and tracks. This supports a direction/pose
or sensor-orientation relationship in the tracker/field setup rather than
a fixed badP3valve. It is not proof thatZalone causes the loss or that paths
occupy an identical3Dregion; pressure amplitudes/trajectories were not logged
as controlled matched inputs. It is not an out-of-volume diagnosis.

The baseline base-sensor quaternion changes by about171.17degrees geodesic
angle relative to the previous test, confirming substantial reorientation.
This is not an independent precise measurement of the user's mechanical
axial rotation; mounting/base tilt and tracking errors can contribute.

Adjacent valid tipZ(min/median/max) at missing events:
previous[-243.77,-222.30,-170.58]mm,current[-250.27,-193.73,-146.79]mm.
The ranges overlap, but missing-transform positions themselves are unavailable.
Raw tracker coordinates are not offset-corrected continuum-arm tip coordinates.

## Combined interpretation

Earlier no-valve-actuation manualP3-direction movement also reproduced30tip
misses. Thus valve switching is not necessary. The rotated test further shifts
the association from a PMA identity toward physical pose/orientation relative
to the tracker setup. Sensor/field quality, externally constrained cable motion
and motion-dependent estimation remain unresolved candidates. Neither cleanCRC
nor clear volume flags identifies a unique cause. No observation validity gate
should be relaxed merely to suppress these reported losses.

The separate residual-actuation oldSwitch/NaN command bug and long GRU re-warm
remain separate controller issues; they cannot explain these independent
Cygnus tests. No hardware/controller/model fix was attempted in this audit.

## Reproducibility and evidence

Valid raw-tip displacement>30mm from a1-5sbaseline identifies response samples;
samples separated by<=8s are grouped. Midpoints of long rests60.463s/96.171s
set broad classification boundaries. User order assignsPMA labels. Counts
2411/1421/1439sum5271runtime reports and0/36/0missing. End of the30mm response
threshold may precede the last small return motion; it is not exactpressureoff.

New reusable offline classifier:RL_data_collection/classifyCygnusStepGroups.py.
Run parser to a fresh folder, then classifier with --order P3 P2 P1 and
--compare earlierFolder, then plotCygnusManualAudit(folder,description).
All operate on saved logs/CSV only. No device or Simulink access.

Retained source_log.txt,summary.json,tracking/missing/eventsCSV,
movement_groups.csv,movement_classification.json,rotation_comparison.csv,
annotated manual_tracking.PNG/PDF/FIG,source snapshots and SHA256manifests.
No runtimeNDI/controller/SLX edit,pressure/hardware/model/trainingrun,
commit orpush by the agent.
