# Continuous PAC tuning: Simulink implementation

Research status2026-10-28: experimental checkpoint; continuous tuning deferred
in favor of established manual runs for residual dynamics model research.
This document preserves integration details, not unattended hardware validation.

R033 review update: invalid completed-motion scores now produce explicit
rejection reason, retain an attempt/checkpoint, and send terminal reply to the
runtime's next request ID. Stop replies do not wait for stationarity. Accepted
observations and completed-motion results are reported separately. Coverage
policy remains95%; R033's88.4% was excluded from GP training.

Implemented 2026-10-28. This supersedes the integration gaps in the earlier
host-only `pacTuneOnlineStep` prototype. The new runtime uses actual reference
phase transitions for timing; it does not use that prototype's host timer.

## Created artifacts

- Root `PAC_autotuning_online.slx`: new wired copy of the inspected saved
  `PAC_autotuning_manual_v2.slx`. The source model is preserved.
- `online_sessions/phi1_T5_scale5_initial/pacOnlineSession.mat`: session,
  bounds, targets and full fixed tuning snapshot.
- Same folder `PAC_online_tuning_blocks.slx`: standalone library with
  `Online controller` and `Settings vector` MATLAB Function blocks.
- `pacOnlineRuntimeStep`: pure fixed-size runtime state machine and scoring.
- `pacOnlineWorker`: separate MATLAB-process BO worker, localhost transport
  and audit logs. No optimization runs inside the control task.
- `pacOnlineInstallModel` / `pacOnlineBuildLibrary`: reproducible construction
  functions; refuse to overwrite outputs. Neither compiles/runs arm models.

The wired arm copy is saved without updating/compiling/running it. The generated
controller library is compiled/tested using new synthetic models only. Hardware
integration, timing, signal validity and zero-command shutdown still require
operator-run smoke verification. No arm experiment has been performed here.

## Initial settings

User-supplied baseline: gammaK=0.0515, gammaD=0.0005999, gammaH=1.49999.
Period 5 s, phi 1 rad, ten scored-motion circles per candidate, score final three.
The 5 s period is the current operational assumption; change before compilation
if another speed is intended. Budget ten completed candidates, including exact
gain repeats. Two consecutive passing batches at identical gains finish the
session. They are operational checks with preserved learning, not independent
reset confirmations.

Initial search ranges (engineering defaults, not proven stable ranges):

| Gain | Lower | Baseline | Upper |
|---|---:|---:|---:|
| gammaK | 0.02575 | 0.0515 | 0.103 |
| gammaD | 0.00029995 | 0.0005999 | 0.0011998 |
| gammaH | 0.749995 | 1.49999 | 2.99998 |

Each proposal is additionally limited to a factor of 1.5 from the last gains.
The GP uses the observations inside that local region; the session budget uses
the whole ledger. The feasibility score is heuristic, not a safety certificate.

The copied tuning chart reads the saved snapshot: hysteresis scale 5, max rate
10 and upper bound 50; the remaining settings are the user's supplied arrays.
Lambda=[18;16.5] and Ks=[15;10], verified in the saved source model, remain fixed.
The old saved source tuning chart used scale 1/rate 5/upper 20; only the new
copy is changed to the selected scale-5 configuration. Mechanical/common
adaptation operator switches must be on and freezeHysteresis off. A runtime
settings mismatch triggers terminal shutdown instead of being ignored.

Objective: fixed-raw velocity-error-rate RMS divided by 0.025 m/s^2.
Acceptance: l2 and l3 RMSE each <=2 mm, vector velocity-error RMS <=0.01 m/s,
roughness <=0.025 m/s^2, coverage >=95%, no tracking limiting or faults.
The score uses independent fixed backward-difference/first-order filters:
measured tau=0.20 s, reference tau=0.10 s. It is a roughness proxy, not measured
physical jerk. Online accumulators were checked against independent offline
`pacTuneFixedVelocity` reconstruction. Position and consecutive-pair coverage
both contribute to the reported minimum coverage.

## Run procedure (operator only)

1. Open and review the new online model. Preserve your physical kill switch.
   Verify sensor/PAC connections, 0.05 s execution, Simulink pacing OFF
   (Desktop Real-Time I/O is incompatible with simulation pacing),
   pressure slew-down behavior and the selected tuning snapshot. Automatic
   candidate changes do not require individual approvals after session start.
2. In a **second MATLAB process**, from this repository root:

   ```matlab
   addpath('passivity_adaptive_control', ...
       'passivity_adaptive_control/autotuning');
   session = fullfile(pwd,'passivity_adaptive_control','autotuning', ...
       'online_sessions','phi1_T5_scale5_initial');
   pacOnlineWorker(fullfile(session,'pacOnlineSession.mat'));
   ```

   The prepared session also includes `start_worker.m` (run in the second
   process) and `open_model.m` (open/review in the model process). Neither script
   starts the simulation or enables pressure. They add the required paths.
   The model opener also calls `pacOnlineConfigureBuild`, routing Simulink
   cache/code generation to short process-specific folders under local TEMP.
   This avoids Windows' 260-character build-file limit when MATLAB's current
   folder is the deeply nested session directory. Call the helper in the model
   process before Run; worker process settings do not affect it. Existing
   generated folders and experiment logs are preserved, not deleted.

3. In the model's MATLAB process, add the same paths and start the model
   manually. The new copy defaults operator pressure enable to **Off**. Set the
   dashboard switch On **before Run** when ready for the test. Off is a terminal
   session stop, not an automatic re-arm; if left Off at Run, the model shuts
   down without starting an experiment. The saved source dashboard referred to
   a nonexistent `Kill` block; the new copy's dashboard is explicitly rebound
   to `Online operator pressure enable`. The source model is unchanged.
   Do not start it using an assistant/hardware automation command.
   The first integration test should be observed by the operator; unattended
   candidate iteration is the implemented mode, not evidence of verified
   unattended hardware safety.
4. Keep normal `out.logsout` logging enabled and save the complete run. The
   worker creates `batch_NNN.mat` checkpoints and `worker_log.mat` inside the
   session directory. Worker packets are UDP audit data and may be incomplete;
   retain the normal model log for full-resolution replay.
5. Do not restart against the same completed session or edit its configuration
   during operation. Logs/lock prevent accidental overwrite. After an abnormal
   termination inspect the lock and logs; there is no automatic resume.

The worker waits up to 30 minutes for the first model packet, allowing preflight
and first compilation. After receiving packets, a 10 s telemetry gap ends the
worker. The runtime's shorter heartbeat check still supervises active operation.
An exited worker saves its attempt even with zero packets; do not delete that
log to bypass the overwrite check. Prepare a fresh session/model instead.

The Java NIO UDP endpoints bind only 127.0.0.1 (ports 32151/32152). No network
service or API key is needed. This is **normal desktop simulation only**, with
an extrinsic MATLAB transport function; it is not a deployable generated-code,
accelerator or Simulink Real-Time target implementation.

Before enabling pressure, run `pacOnlineConnectionCheck(configFile)` in MODEL
MATLAB with the matching worker running in the other instance; require PASS.
R03 includes check_connection.m for this. Stop the old worker before restarting
updated code. Worker startup clears any model bridge retained from its prior
role; the model preflight also clears its own bridge before/after checking.
Probe packets do not start the live telemetry timeout or enter trial logs.
Model UDP initialization failure is explicit warning/reason13. Never run model
and worker in the same MATLAB instance. Keep the current model stopped during
the probe; it sends no hardware command. Save retries under distinct filenames.

## Runtime sequence

Startup preserves the existing pressure/trim/feedback sequence and waits for
fresh valid measurements, PAC status 3, a stationary dwell, worker heartbeat,
and matching settings. Then:

1. Straight hold: accept an ID-matched, bounded proposal.
2. Apply all three gain groups together; verify actual gamma after one-sample
   delay. Hold for 2 s. Estimates/history/filter states are not reset.
3. Quintic bend approach for 5 s; smooth angular acceleration for 2 s.
4. Ten constant-speed circles, scoring the last three; speed ramps are excluded.
5. Smooth angular deceleration for 2 s; quintic bend withdrawal for 5 s.
6. Stationary dwell and BO request for the next candidate, or terminal stop.

Mechanical adaptation is gated off during straight holds and gain application,
and enabled during movement. Common-pressure allocation retains its independent
switch. History continues to evolve through the existing PAC implementation.
The circle reference phase resets only at zero bend, where it produces no
length discontinuity. Fixed length filters are not reset between candidates.

On terminal completion or runtime failure, reference withdrawal runs, then
pressureEnable=false invokes the existing PAC pressure slew-to-zero behavior.
The new Stop Simulation gate waits until commanded pressures are below 0.01 bar
for one second. This confirms commands, not physical pressure sensing. Existing
PAC fault handling and the external kill switch remain independent.

## New block wiring / signals

The installer wires these automatically in the new copy. For manual use, the
controller has four inputs:

| Input | Contents |
|---|---|
| u (12x1 double) | `[t;q2;q3;dq2;dq3;measurementTime;measurementValid;PACstatus;PAClimited;operatorPressureEnable;modelValid;baseValid]` |
| actualGamma (10x1) | Actual tuning gamma vector, Unit Delay 0.05 s; initialize to baseline |
| settings (81x1) | `Settings vector` output through a 0.05 s Unit Delay, initialized to the snapshot |
| networkEnable | Boolean true for worker communication |

PAC status/limited/modelValid are delayed one tick to break output feedback;
raw q/dq are not delayed. Telemetry inputs are converted to double before Mux.
The final delayed tracking-limited sample is checked during first deceleration
sample as well. The original pressure-enable source feeds u(10); controller
pressureEnable replaces PAC input 8. Its adaptationEnable is ANDed with the
original operator mechanical-enable before PAC input 11.

Outputs 1/2 replace all old desired theta/phi destinations. Outputs 3/4/5 replace
the old three gamma source inputs of the tuning chart. Outputs 6/7 gate pressure
and mechanical adaptation. Outputs 8..12 are phase, candidate ID, done, reason,
and `[l2mm;l3mm;velocityRMS;roughness]`. These diagnostics/gains are logged.

`Settings vector` takes scale, rate, lower, upper, initial, adaptationOptions,
historyOptions, pressureOptions, supervisorOptions, lambda, Ks and the three
**pre-gate operator** switches (mechanical, common, freezeHysteresis). It encodes
positive infinite supervisor stop time as -1. Do not feed the intentionally
gated mechanical-enable output back as the operator setting.

Phase IDs: 0 startup, 1 hold, 2 gain acknowledgement/settle, 3 bend approach,
4 spin-up, 5 scored-motion batch, 6 spin-down, 7 bend withdrawal, 9 terminal
withdrawal, 10 pressure-off/done. Reason IDs: 1 clock, 2 operator disabled,
3 PAC fault, 4 invalid grace exceeded, 5 heartbeat, 6 settings mismatch,
7 actual gamma mismatch, 8 tracking limiting, 9 startup timeout, 10 rejected
proposal, 11 proposal/hold timeout, 12 movement during application;
22 target repeated, 23 budget exhausted, 24 worker rejected/failed.

R02 hardware smoke exposed a validity alignment issue: delayed PAC diagnostics
must not extend the continuous raw-measurement dropout timer. Runtime now times
raw current validity independently, retains delayed diagnostics conservatively
for score/hold eligibility, and applies a 1e-9 boundary tolerance. PAC's existing
fault grace remains unchanged. A transient invalid frame during gain hold waits;
actual fresh measured motion still aborts. See r02/review_v1 and
validatePacOnlineValidity; r03 is the fresh same-gain retry session.

## Configuration changes and limitations

Create a fresh session using `pacOnlinePrepareSession(newFolder,g,[g/2 g*2],T)`.
Inspect its saved snapshot; edit the preparation source for other fixed feedback
settings before generating. Rebuild the library and install into a **new** model
copy to bind that session. coder.load configuration is compiled into blocks:
changing MAT contents during a run does not update them. Do not reuse past
scale-1 results as scale-5 BO observations. Filter source identity must be
reviewed in the selected model; the settings vector checks numerical tuning,
not arbitrary changes to MATLAB Function source code.

The worker is independent of the simulation task. The runtime's heartbeat
watchdog can stop motion if the worker hangs, but it runs **inside Simulink**:
it cannot protect against a stalled Simulink process or PC failure. A physical
watchdog/kill mechanism and operator verification remain necessary for that
failure mode. Desktop timing/UDP operations have no hard real-time guarantee.
There is no automatic rollback after a fault, mass estimator, actuator model,
speed schedule, or proven safe exploration/stability guarantee. Preserving
learning makes BO's stationary objective assumption approximate.

Offline checks: `validatePacOnlineRuntime`, `validatePacOnlineLibrary`, and the
three-stage `validatePacOnlineTransport` synthetic fixture. Never run the
transport fixture concurrently with a real session because it uses the same
localhost ports. Earlier `validatePacTuneOnline` still tests the old host-only
prototype; it is not the new runtime integration test.
