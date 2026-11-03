# SineWave3 recurring acquisition losses and diagnostic-output proposal

2026-11-02. Raw source RL_data_collection/20261102/SineWave3.mat; older sine
exports from audits/validity_20261102_v1 used for comparison. Raw recordings
are unchanged. Offline quality/pose analysis only, no hardware/Simulink run,
learned prediction evaluation, parameter fit, threshold relaxation or runtime
reader/driver modification.

## What the recording establishes

12721 sensor samples,14 isolated invalid intervals (0.11005%);13 occur during
channel3, one during channel1 at14.10s. Every native configuration rejection
matches tip valid0/base valid1, orientationValid0, nonfinite relativeXYZ and
NaN tip error. Frames remain fresh/paired and overrun0. Each loss lasts0.05s
simulation time (roughly0.050-0.051s poll-start to recovery).78 native solver
rows are false; those are14 sensor losses, not78 separate failures.

SineWave/SineWave2/SineWave3 have total losses15/20/14 and channel3 losses14/19/13.
This is an uncontrolled setup/repeat comparison, not proof of improvement from
generator relocation. Generator displacement and recalibration were not supplied
as explicit metadata. Initial base raw-position medians differ by about7.4mm
between SineWave2/3 in NDI coordinates; that does not uniquely measure generator
motion because the arm/reference frame may also have changed.

## Is it one specific location?

It is strongly associated with the channel3 motion, but the available neighbor
measurements do not support one single tip position. Exact invalid-frame pose
is unavailable and is never interpolated. Last valid phi before channel3 losses
ranges0.186-2.202rad in SineWave3:

- Near0.19rad:449.75,481.35,512.85,585.80,617.25s.
- Near0.47rad:461.00,492.35,533.75s.
- Near2.12-2.20rad:435.40,475.35,570.50s.
- Two additional first-cycle losses429.75/430.05s follow phi0.397/0.521rad.

For11/13 channel3 losses, other valid samples more than1s away occur within3mm
of that preceding valid XYZ position (21-116 examples for these11 cases).
This is a diagnostic neighborhood count, not an independent statistical test
or proof that the missing pose itself was valid elsewhere. It argues against
treating the whole neighboring region as continuously untrackable.

Position, orientation, cable configuration and motion phase are correlated in
this experiment. A repeated pose-dependent acquisition issue is plausible;
generator distance alone is not identified as the cause. The [NDI-authored
Aurora V2 User Guide](https://www.scribd.com/document/1010309716/Aurora-V2-User-Guide),
sections5.7 and12, describes multiple causes for missing transforms, including
detection-volume limits, interference, motion and cable placement. This primary
manufacturer document is hosted as a mirror; it is context, not evidence of
which condition occurred in this arm.

Current driver reads per-handle status/flags but output logs collapse them to
valid/invalid. It also reads/discards system status and CRC fields without a
verified packet-integrity output. Fresh matched frame IDs do not independently
prove packet integrity. We therefore cannot uniquely distinguish tracker
status, a host rejection, or a more specific hardware/field/cable cause here.

The logged simulated des_phi still has legacy peaks3.005/2.527/2.527rad; this
file does not demonstrate execution of the subsequently corrected bound model.
Those model limits do not directly set raw tip-acquisition validity.

## Next diagnostic extension (proposed, not implemented)

User asks whether NDI code needs diagnostic output. Yes: preserve the reason
from the SAME BX reply before host validity filtering, with one fixed-size
appended diagnostic output while retaining existing output order/meaning.

Recommended per-sensor fields:

- Raw transformation/sensor status code and raw complete port-status word.
- Decoded out-of-volume, partially-out-of-volume and sensor-broken bits.
- Raw-position finiteness, quaternion norm and orientation acceptance.
- A host rejection-reason mask: device/port rejection, repeated/mismatched frame,
  nonfinite data, quaternion/calibration rejection, as applicable to current checks.
- Associate with existing sensor ID, frame number, poll/host time and error
  output. Preserve raw system/packet diagnostics where available; do not label
  CRC as validated unless an actual check is implemented/tested.

No second BX query to fetch diagnostics; read cached values from the pose's
reply. Avoid per-step disk I/O and log through the normal model data path.
Keep pressure, validity, hold/fault and calibration behavior unchanged by this
diagnostic addition. Verify reason masks offline with synthetic status/frame/
quaternion cases before operator installation. No runtime NDI code or SLX
has been changed by documenting this proposal.

Do not spend another full sine run merely to obtain the same binary flags.
The useful next operator test should capture reason codes during a short
reproduction of the channel3 motion. A separate NDI ToolBox test can also expose
vendor status, but it must not compete with our reader for the same serial port.
No experiment is started/queued here.

## Retained artifacts

SineWave3/ exporter tables/provenance; invalid_sensor_events.csv,
native_configuration_events.csv, validity_summary.csv; cross_run_summary.csv;
dropout_previous_valid_poses.csv with all channel3 neighbor poses/commands.
Diagnostics functions/source copies and raw/source hashes are retained.

SineWave3_validity has FIG/vectorPDF/300dpiPNG. dropout_pose_comparison has
FIG/vectorPDF/300dpiPNG: XZ projection of the legacy sensor-relative plotting
frame, grey every tenth valid channel3 sample, colored markers every actual
last-valid sample before a loss. Color is preceding phi. No dropout position
is fabricated. Other-valid-sample counts use full3D XYZ, not projected distance.
Runs are not registered into one fixed world frame; do not infer generator-
volume boundaries from those relative coordinates.

Initial3D raster/FIG attempt is retained separately as3d_attempt; its vector
export failed. Final2D rendering and full numerical tables are the canonical
plot. No analysis values were changed to resolve the export failure.
These are quality-diagnostic paper candidates, not model/control-performance
figures. See verification.txt and source_manifest.csv for provenance/checks.
