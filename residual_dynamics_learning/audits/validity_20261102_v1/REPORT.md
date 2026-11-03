# Recorded configuration-validity diagnosis

2026-11-02. Sources: RL_data_collection/20261102/SineWave.mat and SineWave2.mat,
both recorded model name RL_dataCollection_v2. Offline read-only inspection;
no Simulink/hardware run, prediction evaluation, fitting or normalization change.
These fresh recordings have only been inspected for quality here.

## Finding

Every recorded configuration rejection is explained by invalid tip acquisition
and a nonfinite relative-XYZ input. No configuration-invalid sample occurs with
otherwise valid tip/orientation measurements. The data show isolated losses,
not a continuously false flag.

| Recording | Sensor rows | Invalid sensor intervals | Invalid fraction | During channel3 |
|---|---:|---:|---:|---:|
| SineWave | 12721 | 15 | 0.117915% | 14 |
| SineWave2 | 12721 | 20 | 0.157220% | 19 |

All losses last one sensor interval, exactly0.05s of simulation time. Physical
poll-start to next recovery poll spans0.050024-0.050904s for SineWave and
0.050021-0.050942s for SineWave2. All invalid samples have:

- base sensor valid1, tip sensor valid0;
- orientationValid0 and nonfinite relative sensor XYZ;
- finite equal base/tip frame IDs that advance from the preceding measurement;
- overrun0 and tip error output ndi_raw8=NaN.

SineWave's first loss is14.05s during channel1, then14 occur438.10-611.35s
during channel3. SineWave2's first loss is14.15s during channel1, then19 occur
429.70-619.70s during channel3. All are phase5 (sine excitation). They occur at
several P3 command values, including near zero, so this is not solely a failure
at the maximum-pressure point. Concentration during channel3 is an observation,
not proof of a specific physical sensor failure mechanism.

Example SineWave2: at429.70s tip valid0 -> relative XYZ NaN -> configurationValid0;
at429.75s all recover. Base remains valid throughout the event. The same chain
holds for all35 invalid sensor samples across the two recordings.

## Source trace and limits

Current saved RL_dataCollection_v2.slx routes sensor-relative XYZ into
xyz2ThetaPhi, which calls ndiSupportFnc.configuration. That helper immediately
returns valid=false on a nonfinite input; otherwise it checks the estimated
theta/phi/residual for real finite values. It does not enforce a fit-residual
accuracy threshold. The current chart holds last displayed XYZ-IK estimates
when invalid while retaining configurationValid=false. Thus held smooth outputs
are not evidence of fresh successful estimation. This flag is distinct from
orientationValid, though the recorded losses coincide here.

Driver source hasUsablePosition requires valid sensor status, finite position,
and no out-of-volume, partial-out-of-volume or sensor-broken bits. The MATs do
not log those individual status bits, so the raw acquisition rejection cannot
be narrowed uniquely to missing tracking versus volume/other hardware causes.
NaN tip error and finite paired frames support unavailable tip transforms;
they do not identify a physical cause. Current saved source snapshots are
retained as implementation evidence, not proof of historical source-byte identity.

No repeated/mismatched-frame or flagged-overrun explanation is supported for
these particular drops. No evidence here requires modifying the XYZ solver,
relaxing validity checks, or discarding the complete recordings. Keep invalid
frames masked and exclude windows crossing them. If investigating the physical
cause, focus on tip acquisition during channel3 motion and record detailed NDI
status bits. This report does not qualify the full datasets for training/testing;
timing/domain/coverage checks still apply.

## Artifacts and reproduction

SineWave/ and SineWave2/ hold existing exporter outputs and provenance dictionaries.
invalid_sensor_events.csv lists all35 events with pressures, flags, frames,
trial/phase/channel and physical times. native_configuration_events.csv lists
native flag intervals/recovery times; validity_summary.csv gives totals.
Native configuration runs at variable solver times:80 and107 logged false
rows correspond to15 and20 sensor losses, respectively, not80/107 independent
sensor failures. No event is filtered away based on model prediction quality.

SineWave_validity and SineWave2_validity have editable FIG, vector PDF and300dpi
PNG. Each shows configuration flags, P3 commands and the first channel3 event
at full sensor resolution. Red markers identify all dropouts. These are quality
diagnostics/paper candidates, not control or model-performance figures.

diag source: diagnoseRecordedRlValidity.m. After running exportRlTrainingData
for the two source files into fresh per-run folders, call
diagnoseRecordedRlValidity(auditFolder). Native log checks use a1ns numerical
timestamp comparison tolerance because nominally coincident solver/sensor times
differ by about5.7e-14s at one boundary. This is not lag correction or filling.
Initial diagnostic/export development failures are preserved in verification.txt;
the final pass checks every native configuration rejection against independently
logged tip/orientation and XYZ channels. Both final plots were visually checked.
source_manifest.csv contains hashes for raw MATs and current model/source
snapshots. No source MAT, SLX, acquisition code, model weights or controller changed.
