# Open-loop pressure sources for real-arm dynamics data

2026-11-02 NDI diagnostics: the updated reader appends a32x2 base/tip status
matrix (BaseRotation port23, Position_only port20). With the model stopped,
refresh the System class and log the new last output as ndiDiagnostics.
See ndi_diagnostics/README.md for wiring/schema and offline tests. Exporter
adds named raw-status/rejection fields when present; existing logs are supported.
No SLX wiring or hardware run performed by the assistant.

Paste **one** of these standalone .m functions into a MATLAB Function block.
Inputs: sampled simulation time `t` in seconds (existing .05 s sampled clock).
Set explicit discrete sample time .05 s for predictable labels. No workspace
parameters or helper function dependencies are needed inside the source blocks.
Edit the constants near the top before a run; save the exact source/settings
with the trial. Time is simulation time, not wall time; retain NDI host timestamps.

Connect pressureRequest_bar [3x1] to the existing **Pneumatic Valves/Pressure input**
in place of the manual three-pressure Mux. Keep the downstream Kill and saturation
path. If comparing the forward model, feed it the same request separately.
The inverse-model prediction is not connected to hardware by these blocks.
There is no deadzone subtraction, inverse model or feedback in these sources.
The model prediction's pressureUsed_bar is not the dataset action.

## Choices and defaults

| Block | Default experiment | Duration |
|---|---|---|
| rlPressureStepsBlock | Sequential P1/P2/P3: peaks 1,2,3 bar; two repetitions, rise5/hold10/fall5/zero10 s; initialzero5 s | 545 s |
| rlPressureSineBlock | Sequential P1/P2/P3: 0..3 bar, omega=.2 rad/s, three cycles per trial, two trials per channel,10 s relaxation; initialzero5 s | 630.487 s |
| rlPressureCombinedBlock | Four modes: all-channel common pulse, symmetric s, antisymmetric a, rotating s/a; two cycles each at .2 rad/s; rise/fall5 s,relax10 s; initialzero5 s | 336.327 s |

For the first two blocks, channelSelection=1/2/3 selects P1/P2/P3. Set it to0
to run each muscle separately in one sequence. This triples the number of
trials while retaining one initial baseline. Defaults select0: all three sequentially. Select1 to reproduce the original P1-only schedule.
For sine, period=2*pi/omega=31.416 s; 30 s is not one complete cycle.
Combined defaults: common pulse0..2bar; differential bias1.5bar and amplitude.5bar.
The rotating differential mode is pressure excitation, not a commanded tip circle.
All values are raw pressure requests in bar within the model's existing0..3bar
command envelope. These are not newly validated bounds for every configuration;
begin with the already-tested single-P1 sine before combined-channel collection.

## Log the outputs

All blocks have the same output order:

1. pressureRequest_bar (double3x1): pre-Kill/pre-saturation requested pressures.
2. trialId (double scalar): 1-based trial;0 during initial baseline/completion.
3. phaseId (double scalar): -1 invalidtime,0 baseline,1 loading/preload,
   2 static hold,3 unloading,4 zero relaxation,5 periodic excitation,9 complete.
4. channelId (double scalar): 1/2/3 single channel;0 combined/baseline/completion.
5. cycleId (double scalar): sine/combined periodic cycle, or step repeat number;
   zero outside that labeled portion (step repeat also labels its relaxation).
6. done (Boolean): true after the programmed sequence; request remainszero.
7. testEnd_s (double scalar): total duration for the current settings.

Set Stop Time >= testEnd_s (optionally add5 s) to retain final relaxation.
done does not stop Simulink or replace the Kill. Do not use indefinite repeating
excitation for a trial that needs an identifiable ending.
Requested outputs are bounded by settings assertions, not silently clipped.
The physical command path remains responsible for its configured limits.

Record the final P1_NI_DAQ/P2_NI_DAQ/P3_NI_DAQ, not just this source. Also retain
measured orientation, XYZ plus corrected moduleTipXYZ_m, inferred lengths,
validity/observability, frames, host/poll timestamps, base orientation and Kill.
Pressure values are commands/proxies, not measured chamber pressure. Constant
signals may log sparsely; export must preserve their held-value semantics.
Do not label des_theta/des_phi predictions as physical targets or measurements.

No MAT files are saved automatically by these blocks. Save `out` plus trial/session,
source parameters, calibration and mounting/load metadata after each run.
`validateRlPressureBlocks` checks default waveform bounds, complete endings,
invalid-time handling, source slew and sine equivalence offline only.


## Angle-reference experiments (not pressure sources)

New standalone MATLAB Function blocks take the same sampled time `t`:

| Block | Default schedule | End time |
|---|---|---|
| rlCircleSweepReferenceBlock | phi .2,.8,1,1.2,1.5 sequentially; each rise5s, target wait3s, five20s circles, fall5s, zero hold10s; initialzero5s | 620 s |
| rlHelixReferenceBlock | Four trials; each phi0->1.5 over three20s turns, then1.5->0 over three20s turns, followedzero10s; initialzero5s | 525 s |

These output [theta,phi,trialId,phaseId,cycleId,direction,done,testEnd_s].
Theta/phi are desired angles in radians, not measured response or pressure.
Theta is unwrapped and linear during rotation, with no per-lap pauses. Theta
speed starts/stops at experiment phase boundaries; phi rise/fall uses smoothstep.
Circle direction defaults+1; helix alternates+1,-1,+1,-1 across trials atphi0.
The helix has continuous rotation at its peakphi; no midpoint stop. This is an
outward/inward, geometry-constrained tip spiral, not prescribed constant-pitch
Cartesian XYZ. Change turnsLoading/turnsUnloading/numTrials inside the block.

Circle phaseIds:0baseline,1rise,2targetwait,3fall,4zerohold,5circle,9done.
Helix adds6loading/7unloading. -1 indicates invalid time.

IMPORTANT: Do not connect theta/phi to the valve-pressure vector. To realize
these requested paths, a separate angle->desired length/rates->feedforward
inverse pressure path is required for open-loop collection. A closed-loop
tracking controller would produce controlled data, outside the current stated
open-loop plan. The saved validation model's simulated des_theta/des_phi are
forward predictions; substituting scope/comparison references alone will not
actuate the arm along these paths. Actual open-loop paths need not match them.
Log requested angles separately from forward-model predicted angles.

Set Stop Time to625s(circle sweep) or530s(helix), retaining the finalzerohold.
The existing actual pressure limits/Kill and timestamp logging still apply.
Defaults reachingphi1.5 describe requests, not a new experimentally verified
motion/pressure limit. validateRlReferences checks schedules, continuity and
rotation speed offline; it does not establish physical feasibility or tracking.

## CSV/report export

After saving a raw MAT containing `out`:

```matlab
addpath('RL_data_collection');
r = exportRlTrainingData('RL_data_collection/SineWave.mat');
% Or process all top-level MAT recordings in this folder:
% results = exportRlTrainingData('RL_data_collection');
```

An existing export is not overwritten; a timestamp suffix is used for another
export of the same file. Original MAT files are never modified. Nothing runs
Simulink or hardware. Outputs go to `RL_data_collection/exports/<recording>/`:

- REPORT.md: coverage, measured trial ranges, timing and interpretation.
- observations.csv: one row per native NDI orientation sample, with pressure
  commands aligned using previous-sample hold in simulation time. Explicit
  phi/theta/fresh-frame/transition validity flags; invalid angle targets are NaN.
- commands.csv: final NI commands on the union of their native log timestamps.
- model_predictions.csv: computed des_theta/des_phi, simulated lengths and
  inverse/effective-pressure predictions; keep separate from measured targets.
- experiment_intervals.csv and trial_summary.csv: labels from the recording,
  never schedules guessed solely from a filename. Absent if labels unavailable.
- acquisition_gaps.csv: host sample intervals above the configured threshold.
- signal_inventory.csv and column_dictionary.csv: exact source names, index,
  BlockPath, declared units and interpretation. Native multidimensional data
  flatten in MATLAB column-major order after the time dimension is moved first;
  rotation elements are R11,R21,R31,R12,R22,R32,R13,R23,R33.
- trial_manifest.json: source path, quality, assumptions, metadata and inventory.
- native_signals/: every numeric logged timeseries with ORIGINAL timestamps and
  values, including duplicates and invalid/held values. Indexed filenames avoid
  collisions. Original MAT remains the authority for unsupported signal types.

Duplicate aliases are resolved by known physical-reader/valve provenance. An
unresolved ambiguity raises an error rather than silently selecting a signal.
Constant-valued signals may be held after their final sparse log; every use is
noted. Nonconstant signals are not extrapolated beyond coverage. The default
physical transition-gap threshold is .10 s; set
`options.maxTransitionHostInterval_s` if a different criterion is justified.
That threshold is a quality filter, not a change to the hardware sample period.
Do not automatically use all rows as training transitions. Histories crossing
invalid samples, long gaps or trial/session boundaries require explicit handling.
The exporter does not assign train/test splits or invent missing metadata.

Supply experiment metadata explicitly when available:

```matlab
options.metadata = struct( ...
    'session_id','session_20261019_A', ...
    'trial_id','SineWave', ...
    'experiment_type','sequential_single_channel_sine', ...
    'module_id','ENTER_MODULE_ID', ...
    'calibration_file','ENTER_ACTUAL_CALIBRATION', ...
    'payload_g',0, ...
    'base_motion','fixed', ...
    'reset_notes','ENTER_ACTUAL_RESET_PROCEDURE');
% Also save the exact source code/parameters and mounting/frame convention.
% Do not copy metadata values to a new session without checking them.
r = exportRlTrainingData('RL_data_collection/SineWave.mat',options);
```

Set options.writeNative=false for a smaller summary-only package. Raw measured
pressures are unavailable: command values remain setpoint proxies, not sensors.
Sensor host elapsed times and simulation time are both retained, but the actual
NI-write host timestamps are absent; this export does not magically synchronize
physical action delivery to sensor acquisition or identify chamber-pressure lag.
Split independent evaluation by whole trials/sessions, not random adjacent rows.

The saved inverseModel input is [l2;l3;dl2;dl3;ddl2;ddl3], not [theta;phi]. For feedforward circle/helix excitation, compute desired lengths/rates and feed that vector to the inverse; use its pressure_bar output through the valve safety path. Existing angular references alone are not complete pressure generators. Ensure smooth angular start/stop and valid derivative generation. Actual open-loop tracking is not guaranteed.

Circle/helix exports now include theta_reference_rad and phi_reference_rad from theta_from_block/phi_from_block. If the source is feedforward inverse dynamics, pressureDesired_bar is retained as the pre-valve request; final P1/P2/P3_NI_DAQ remain the applied command proxies. Forward-model des_theta/des_phi stay in model_predictions.csv. All five saved experiment types were checked offline.
