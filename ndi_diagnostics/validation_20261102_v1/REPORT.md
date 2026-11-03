# NDI diagnostic verification: 2026-11-02

Offline only: no serial/DAQ device opened, System reader setup/preflight called,
SLX edited or hardware model run. Real driver code uses an injected synthetic
byte transport; initialization/detection communicate only with mock responses.
CRC fixture values are received data, not checksum validation.

## Checks passed

- 32 cases,16 per reader: exact pre-change/current step bodies, matching packets.
- Existing XYZ/error, sensor/frame, held/recovered orientation/length,
  observability/nonbending and base attitude/gravity/validity outputs match.
  Host clocks excluded; pacing reset in fixture, no real-time qualification.
- Missing/recovery, volume/partial/broken, repeated/mismatched frame, bad
  quaternion, missing calibration, disabled/omitted handle, nonfinite position,
  straight pose, counter wrap and normal measurements.
- Existing output metadata plus appended32x2 double/fixed/real port contracts.
- One BX per sample; cached getter sends none. Raw status/quaternion data
  precede output masking and belong to the same reply.
- Disabled/omitted current status stays unavailable/NaN, not stale.
- Device-valid position/invalid quaternion distinguished from acquisition loss.
- Incomplete packet raises unchanged IncompleteBX exception; cache marks
  incomplete and recovers with the next full reply.
- Raw-only reader stage remains unknown; named rejection decoder passes.
- Synthetic MAT with/without diagnostics has identical old observation values.
  Matrix flattening maps64 named base/tip fields and exports dictionaries.
- Updated exporter preserves all12721 original SineWave3 observation rows.

case_results.csv, verification.mat and verification.txt retain outcomes.
before/ and before_source_hashes.csv preserve original sources; source_snapshot/
and source_manifest.csv preserve corrected production/test sources.
reader_replay/ contains extracted standalone test bodies, not runtime readers.
Synthetic MAT/export fixtures and derived legacy export are not new experiments.

A supplemental decoder test command initially failed shell quoting before
execution; corrected strcmp invocation passed. Reader/source data were unaffected.

## Reproduce

From repository root, copy the old comparison sources into a NEW folder:

    addpath('ndi_diagnostics');
    copyfile(fullfile('ndi_diagnostics','validation_20261102_v1','before'), ...
        fullfile(newFolder,'before'));
    validateNdiDiagnostics(newFolder);

Transcript also records supplemental decoder and actual legacy export checks.
No hardware/model execution required. See ../README.md for logging instructions.

## Limits

Actual System-block port refresh/compilation, acquisition overhead and hardware
timing remain unverified. Headers/length/CRC are received, not newly validated.
Status narrows the rejection stage, not necessarily a unique physical cause.
New operator data is needed to diagnose actual channel3 tip losses.
