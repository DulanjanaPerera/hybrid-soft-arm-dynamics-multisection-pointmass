# NDI same-reply diagnostics

Implemented 2026-11-02 after repeated tip 0B losses in the channel-3 sine trials.
The diagnostics identify acquisition/host-rejection conditions; they do not
infer a unique physical cause for a missing transform.

## Reader output and operator installation

Existing output order, shapes and validity/hold logic are retained. One real,
fixed-size **32x2 double matrix** is appended:

| MATLAB System class | New output port | Signal name |
|---|---:|---|
| AuroraReader_s2_BaseRotation | 23 | ndiDiagnostics |
| AuroraReader_s2_Position_only | 20 | ndiDiagnostics |

Columns are **base, tip**, in the configured BaseSensorId/TipSensorId order.
There is still one BX query per reader step. No additional polling or per-step
file writing is introduced. Calibration/preflight and existing pressure/control
connections are unchanged. Raw driver data are captured before output masking.

With the model stopped and current recordings saved, reload/update the MATLAB
System class/ports in your MATLAB session. If a loaded class definition prevents
the update, use a fresh session after saving workspace data. Keep existing
connections; connect the new last port to a Terminator, name its signal
**ndiDiagnostics**, and enable signal logging so it is saved in out.logsout.
The project root must be on the MATLAB path. No SLX file was edited by the
assistant. Actual model port refresh/compilation and hardware timing still
require operator verification; no sensor or hardware model was run here.

## Most useful rows

| Row | Field | Interpretation |
|---:|---|---|
| 3 | reply_counter | Host BX query count; may include preflight polls |
| 4 | handle_present | Whether this handle occurred in this reply |
| 5 | sensor_status | Numeric raw status:1 valid,2 missing,4 disabled; NaN if absent |
| 6 | port_status_word | Raw current complete port-status word, stored exactly as double |
| 7 | port_status_available | Whether that word was present in this reply |
| 8-10 | volume/partial/broken flags | Existing decoded status bits; NaN if unavailable |
| 11 | raw_position_finite | Before driver replaces unusable positions with NaN |
| 12-13 | quaternion norm/valid | Before driver output masking; norm tolerance unchanged0.02 |
| 14 | driver_usable | Existing driver acceptance result |
| 15-16 | frame/error | Same-reply frame and tracker indicator |
| 17-20 | frame/sample checks | Existing finite/change/pair/sample conditions |
| 21-22 | calibration/orientation | Existing calibration availability and pair orientationValid |
| 23 | reject_mask | Host condition bitmask; decode below |
| 24 | system_status | Raw same-reply system word |
| 25 | crc_checked | Always0: checksum verification is NOT implemented |
| 26-29 | received packet metadata | Start, length and received header/body CRC values |
| 30 | reply_read_complete | Expected fields consumed, not proof of integrity |
| 31-32 | poll/host times | Existing acquisition timestamps, duplicated per column |

ndiDiagnosticsSchema() returns the complete fields and bit-value tables.
Unknown or omitted status remains NaN rather than reusing the prior reply's
flags. Decoded volume/broken bits use the same PortHandle rules as before.
Frame change retains the existing ~=previous rule, including counter wrap.
An overrun remains a separate original output and is not a new rejection gate.

## Reading the reason mask

    [sensors,conditions] = decodeNdiDiagnostics(D);

D is one 32x2 sample. sensors contains named fields; conditions lists all set
bits. Examples:

- tip sensor_status2: tracker reports a missing transform; not enough to decide
  field interference versus another physical cause.
- tip sensor_status1 with out_of_volume1: returned pose rejected by the
  existing volume rule; raw quaternion diagnostics may still be available.
- driver_usable1/sample_valid1 but quaternion_valid0: orientation-stage rejection.
- repeated/mismatched frames: existing host freshness/pairing rejection.

| Value | Condition |
|---:|---|
| 1 | Handle absent from reply |
| 2 | Sensor transform status not valid |
| 4 | Out of volume |
| 8 | Partially out of volume |
| 16 | Sensor-broken flag |
| 32 | Nonfinite position despite transform status1 |
| 64 | Nonfinite frame |
| 128 | Repeated frame |
| 256 | Pair frame comparison rejected, including NaN |
| 512 | Quaternion rejected at the orientation stage |
| 1024 | Calibration unavailable |
| 2048 | Orientation estimator rejected despite valid sample/quaternion prerequisites |
| 4096 | Pair orientation unavailable because either sample is invalid |

The mask reports conditions, not mutually exclusive physical diagnoses.
Base may have pair_sample_invalid because the tip failed. Conversely base mask0
does not establish pair orientation validity; inspect row22 and the tip too.
A raw driver-only snapshot has no reader-stage mask yet (NaN), which the
decoder labels reader_stage_unavailable rather than interpreting it as success.

The driver preserves same-reply system/header/CRC values but does not validate
checksum, header or length beyond the previous field reads. No packet-acceptance
rule was relaxed or strengthened. Existing incomplete-read exceptions continue
to propagate; the cache marks partial replies incomplete. A diagnostic output
is not guaranteed for a step that throws before returning. getLastDiagnostics
reads the cache only and never issues another query; readSensorStatus still
polls and must not be used to assemble same-frame diagnostics.

## Export and compatibility

exportRlTrainingData detects a logged ndiDiagnostics signal at the physical
MATLAB System path, validates schema/shape and exports named columns, including:

    ndi_base_sensor_status
    ndi_tip_sensor_status
    ndi_tip_out_of_volume
    ndi_tip_reject_mask

It also writes ndi_diagnostic_fields.csv and ndi_rejection_reasons.csv.
Existing acceptance masks continue to use the original logged flags, not the
new diagnostics. Recordings without the new signal remain supported.

For a fresh saved MAT:

    result = exportRlTrainingData(recordingPath);
    O = readtable(fullfile(result.folder,'observations.csv'));
    bad = isfinite(O.ndi_tip_reject_mask) & O.ndi_tip_reject_mask~=0;
    O(bad,{'time_sim_s','ndi_tip_sensor_status','ndi_tip_port_status_word','ndi_tip_reject_mask'})

No diagnostics can be recovered retroactively from the old Boolean-only logs.

## Verification and boundaries

See validation_20261102_v1/verification.txt, case_results.csv and before/ for
pre-change sources. Tests used a byte-accurate synthetic transport with the
real decoder, cached getter, finalizer, both saved/current reader step bodies,
output-metadata probes and synthetic MAT export. The 32 reader cases cover missing,
disabled/omitted handles, volume/partial/broken flags, repeated/mismatched
frames, bad quaternion, missing calibration, nonfinite position, straight pose,
counter wrap, incomplete-packet/recovery and normal measurements.

Existing reader values/types/holds match saved pre-change bodies in those
cases; host clock values are excluded, and real-time pacing is not qualified.
No serial object/DAQ is opened, reader setup/preflight executed, or model run
by the tests. Constructor transport injection is an offline test seam; normal
one-argument AuroraDriver_2 construction is unchanged.

Future physical tests should first reproduce a short known channel3 motion
with these fields logged. Preserve full raw out, setup/calibration information
and new filenames. Status codes may identify the rejection stage without
identifying a unique underlying hardware/field/cable cause.
