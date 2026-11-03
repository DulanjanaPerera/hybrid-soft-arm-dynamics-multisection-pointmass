# SineWave4: first recorded same-reply diagnostic review

2026-11-02. Raw source RL_data_collection/20261102/SineWave4.mat, recorded
model RL_dataCollection_v2. Read-only offline export/native-log/diagnostic
review; no hardware/model execution, source gate changes, learned prediction
evaluation or fitting.

## Main result

The new diagnostic output is present, schema1,32x2. All12721 diagnostic frames,
sample-valid outputs, orientation flags and physical timestamps match the
original logged channels; reply counter advances exactly1 per sensor sample.
This recording demonstrates operator-installed diagnostic logging. It does not
validate every diagnostic branch or future hardware timing.

There are10 isolated tip losses, all during channel3 sine excitation:
10/12721=0.078610%, orientation valid99.921390%. No channel1/2 losses.
Each lasts0.05s simulation time; physical poll-to-recovery intervals are
0.050038-0.050868s.50 false native solver rows correspond to10 sensor losses.

Every loss has:

| Field | Value |
|---|---|
| Base sensor status | 1, valid |
| Tip sensor status | 2, MISSING |
| Tip raw port word | 63 =0x0000003F |
| Out-of-volume / partial / broken | 0 /0 /0 |
| Base reject mask | 4096, pair_sample_invalid |
| Tip reject mask | 4098 =4096+2, pair_sample_invalid + sensor_status_not_valid |
| Tip raw position / quaternion norm | unavailable/nonfinite /NaN |
| Frames | finite, changed, paired |
| Overrun at loss | 0 |
| System status | 0 |
| Reply fields consumed | complete |
| CRC checked online | 0, not verified |

The raw tracker/driver status says the tip transform is missing before our
quaternion/angle/frame gates could reject a returned transform. Quaternion norm
NaN is a consequence of no transform, not a measured bad quaternion.
The [Combined API implementation documentation](https://github.com/PlusToolkit/ndicapi/blob/master/ndicapi.h)
describes MISSING as inability to compute the tool transform and documents
port initialized/enabled/volume fields. The code/status interpretation is not
a unique diagnosis of the physical mechanism.

## Packet consistency

Normal replies have reported body length87bytes (12711samples), missing-tip
replies55bytes (10samples). The32byte difference is exactly four quaternion
float32 values, three position values and one indicator value omitted by the
missing-transform branch in the current decoder. Header start word is42436
(0xA5C4) throughout; received header CRCs are8979/9019 for the two lengths.
This is coherent missing-transform framing, supporting a device-level result
rather than simply our angle calculation setting its flag false.

The port word also equals0x3F on4705 valid-tip samples; it alone is not a
missing-transform indicator. Other valid-tip words are0x0002003F. We retain
uninterpreted higher bits rather than infer undocumented meaning. The
initialized/enabled/in-port bits remain set.

All reply-read-complete fields are1, but online CRC verification is not
implemented. These observations do not prove full body integrity or rule out
every parser/communication mechanism. No checksum flag is retroactively changed.

## Where it happens

Loss times(s):430.55,430.85,447.60,479.00,492.95,510.35,566.05,583.20,614.65,614.75.
Four on pressure increase, six on decrease, all trials5/6.
Last valid phi before losses0.64775-0.88365rad. Example:
430.55s: precedingphi0.70875, P3=1.71077bar, status2; recovery430.60s.

Loading losses occur P3~1.58-1.80bar; unloading~0.884-0.926bar. Similar bend
region at different pressure on the two directions supports pose/history
association rather than a single pressure threshold. It does not distinguish
field/orientation/cable/motion effects. The exact missing-frame pose is unknown;
only actual neighboring valid samples are retained, not interpolated.

Generator setup/calibration/cable changes were not recorded as controlled
factors. Previous sine loss counts15/20/14 versus10 here are descriptive, not
causal evidence that moving the generator or adding diagnostics fixed tracking.
No reported volume flag is NOT proof that the sensor is inside every effective
detection region: MISSING is generic and no coordinates were computed.

## Other checks and next action

- Seven overrun flags and five host gaps>.10s occur elsewhere: gaps at0.05,
  0.10,37.40,39.70,415.55s; none coincide with these10 losses. Preserve existing
  gap segmentation/warmup/history guards in offline preparation.
- Twelve orientation-valid samples have thetaObservable=false near straight;
  they are separate from acquisition losses.
- Logged forward phi peaks are2.526838/2.526838/2.526837rad across channels,
  consistent with the corrected three-actuator bound implementation. This does
  not independently establish installation of every inline inverse/controller.

Do not loosen validity/quaternion/volume checks or discard the whole recording
because of10 one-frame gaps. Exclude invalid/gap-crossing windows and retain
original data; full training/evaluation eligibility is a separate audit.

If investigating physical cause, use a short operator-controlled reproduction
through this moderate channel3 bend range and compare slower motion/stationary
hold and cable/field arrangement one factor at a time. An independent
[NDI ToolBox check](https://www.ndigital.com/software/) can help separate vendor
tracking behavior from custom-decoder behavior; it must not compete with our
reader for the same COM port. MISSING alone will not name a unique mechanism.
No such test was run or queued by the assistant.

## Artifacts

SineWave4/ retains full export/provenance/field dictionaries. Native intervals,
invalid_sensor_events.csv, missing_transform_events.csv, diagnostic_summary.csv,
packet_size_distribution.csv and verification.txt retain every event and check.
SineWave4_validity and same_reply_missing_diagnosis have FIG/vectorPDF/300dpiPNG.
The latter shows P3/measuredphi, one raw status/validity event, and last-valid
relative XZ poses. Markers are preceding measurements, not missing-pose estimates.
Grey spatial path is display-decimated; numerical source remains full resolution.

Source snapshots/hashes include analysis/decoder/schema, current raw/status
driver and readers, and original MAT hash. Analysis script asserts this run's
recorded findings; it is not a generic future-log classifier. No raw, controller,
SLX, weight, calibration, threshold, commit or push change.
