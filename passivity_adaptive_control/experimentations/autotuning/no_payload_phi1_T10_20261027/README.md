# First operator-run framework trials

Prepared 2026-10-27. Model duplicate supplied by user: PAC_autotuning_manual.slx.
User-selected conditions: no payload, phi=1 rad, period=10 s/circle. No trial has
been executed by the assistant. This is a proposed sequence, not recorded results.

## Sequence

| Order | Filename | Circles | Tracking time | Purpose |
|---|---|---:|---:|---|
| 1 | smoke_00g_phi1_T10_3cycles_r01.mat | 3 | 30 s | Verify startup/logging/offline import |
| 2 | baseline_00g_phi1_T10_12cycles_r01.mat | 12 | 120 s | Observe learning and delayed roughness |
| 3 | baseline_00g_phi1_T10_12cycles_r02.mat | 12 | 120 s | Repeatability with fresh software state |

Run trial 1 first and have its log reviewed before trials 2/3. Review baseline r01
before repeating if it shows growing oscillation, faults or other abnormal behavior.
These are baseline characterization runs, not formal optimizer confirmation: gain
bounds, quantitative acceptance targets and physical reset/preconditioning are not
fully established yet. Do not start optimizer-selected trials at this stage.

## Keep fixed and record

- Use the current reproducibly smooth gains, all tuning vectors and switch values.
  Do not substitute the older disk defaults or change lambda/Ks/gamma mid-run.
- Keep adaptation/common-pressure/freeze switch settings identical across runs;
  log their actual values. The first run establishes the observed configuration.
- No payload must be present during calibration/trim and tracking. Record the
  calibration procedure and whether trim is redone for each trial.
- Start each run fresh so parameter estimates/history/supervisor state reinitialize.
  Do not resume learned state. Use the same settling/preconditioning procedure;
  software restart does not itself reset material or pneumatic history.
- Keep phi=1 rad, T=10 s, circle direction, startup and smooth return consistent.
  Change the trajectory's circle count, not just the model stop time while a
  longer reference is still circling. Retain the existing operator stop method.

If the verified circle start is 20 s, tracking windows are [20,50) for the smoke
run and [20,140) for each baseline. Let the existing smooth return/shutdown finish
after tracking; the total model stop time must include its actual duration.
These timing assumptions must be checked from logged reference, not filenames.
Evaluator defaults [20,320] are for old 30-cycle logs and must not be used here.

## Save the original out variable

Save once the manually run simulation/acquisition has finished, from the MATLAB
workspace containing that run's out (which must include out.logsout):

```matlab
dataDir = fullfile( ...
 'C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl', ...
 'passivity_adaptive_control','experimentations','autotuning', ...
 'no_payload_phi1_T10_20261027');
file = fullfile(dataDir,'smoke_00g_phi1_T10_3cycles_r01.mat');
assert(~isfile(file),'Choose a new repeat number; do not overwrite a recording.');
assert(exist('out','var')==1,'Save the completed run output as out first.');
save(file,'out','-v7.3');
```

For baseline runs change the filename to the matching table entry. Preserve
aborted/faulted runs under a distinct filename (for example append _aborted), and
record what happened; do not discard or overwrite them with a later successful run.

## Logging and notes

Keep existing full Dataset logging, including reference/measured lengths and
rates, angles/observability, validity/freshness, status, final pressures, nominal
residual, limiting, adaptive estimates/history, all ten tuning-source outputs,
lambda/Ks and adaptation/common-pressure/freeze enables. Retain s_mPerS,
tauBlended_N and forceBias_N for gate/allocator replay, force components for model
diagnostics, and compensated module-tip XYZ/validity for later spatial plots.
The evaluator checks exact names/provenance and reports missing/ambiguous data.
Do not silently rename or select arbitrary duplicate signals.

Add a same-stem notes.txt with: date/time, no-payload confirmation, calibration/
trim procedure, direction, initial/reset/preconditioning procedure, exact gains
if not logged, observations/interventions and whether the full return completed.
Record unknown details as unknown. Do not claim a verified reset or completed
trial if it was not checked. The formal campaign manifest can be prepared after
reviewing the first recording and agreeing the baseline/targets.

The assistant has not opened/modified/run the model, applied gains or created
a configured hardware campaign. Send the smoke MAT path for the first review.
