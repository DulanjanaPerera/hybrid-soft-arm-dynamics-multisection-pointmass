# NDI diagnostic-output implementation

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
