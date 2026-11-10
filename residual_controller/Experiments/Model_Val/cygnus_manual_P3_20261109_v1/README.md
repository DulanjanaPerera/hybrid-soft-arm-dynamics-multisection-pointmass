# Manual P3-direction motion reproduces tracker loss without valve actuation

2026-11-09. User provides a NEW Desktop/log.txt and reports moving the arm
manually in the P3 actuation direction with no pneumatic valve actuation.
This physical test condition is user-reported, not encoded in the tracker
protocol. Source3640632bytes,SHA256
a6a8cd6c8c30e24645b8519762691e560db93ff18dabe6494c9d1442fbe0d118.
HeaderMonNov09 17:50:08CST2026; lastwrite17:52:23. Preserve separately from
the previous manual and pressure-driven logs, all namedlog.txt originally.

## Result

**Yes: explicit tip-MISSING transforms recur during this no-actuation test.**
3035BXreplies,allheader/payloadCRCsandlengthchecks pass. Finalsession3033
replies over76.0373s; first2initial acquisition reports excluded.

| Runtime after startup | Base0A | Tip0B |
|---|---:|---:|
|Reports|3031|3031|
|Valid transforms|3031|3001|
|MISSING transforms|0|30|

The30tip losses(.990% of runtime reports) form28episodes:26single-update
and2two-update episodes. To next valid reported sample they span24.979ms
or49.959ms. All30have valid base transforms and matching base/tip frame
counters. No out-of-volume or partial-volume flags, systemstatus0throughout.
No repeated runtime frame counters or identical consecutive valid tip poses.
These are explicit handle status0x02missing, not just gaps in host queries.
They occur30.5498-69.9173seconds after the first final-session BX reply.
Detailedsource lines/frame numbers are in missing_tip_reports.csv and
missing_episodes.csv; no interpolation or held values fill invalid transforms.

The log's Track Frequency40.033Hz is distinct from Frame Frequency320.264Hz;
normal frame increment8is24.979ms. Final-session query increments are8x3023,
16x8and40x1. Query skips are distinct from explicit MISSING status. Missing
transform durations here refer to returned observations, not unknowable
validity of any unqueried intermediate frames.

## What this establishes

Under the user's stated no-valve-actuation condition, valve actuation is
**not necessary** to trigger the tip tracking loss. Therefore the earlier P3
losses cannot be explained solely by pneumatic valve switching. This does
not prove valves/electronics can never contribute in another configuration,
or that all nearby powered hardware was off. It does not identify a unique
physical cause. No evidence here establishes out-of-volume; the user's
inside-volume statement is consistent with all reported flags.

Three preserved captures now give:

| Capture | Runtime tip MISSING |
|---|---:|
|Pressure-driven standard sine, inferredP3phase|23|
|Earlier generic fast manual motion|0|
|New manual P3-direction motion without actuation|30|

This supports investigating direction/pose/orientation, sensor/cable motion
or measurement quality, rather than assuming a valve-switching trigger or a
universal motion-speed limit. Tip-only losses with valid base and cleanCRC
also argue against treating these records as entire-system serial corruption.
They do not rule out all acquisition or communication issues elsewhere.

## Pose comparison with earlier pressure-driven P3 losses

The new trajectory covers a broader range but some adjacent valid poses
overlap the earlier failure region. Distances compare base-frame raw sensor
separation, correcting for base translation/orientation, not calibrated
continuum-arm tip offsets.52unique valid poses adjacent to the new missing
reports versus30adjacent valid poses from the earlier sine capture give
nearest-previous-failure-neighbor distances min1.125mm,median12.197mm,
max88.284mm. This supports some overlap, not identical matched trajectories
or a unique small spatial failure region. Positions during the missing
transforms themselves are unavailable.

New adjacent-valid raw tipXYZ ranges(mm):
X[-29.33,106.45],Y[-164.96,-97.83],Z[-289.61,-184.91]. Earlier pressure-P3
adjacent-valid ranges:X[31.22,63.03],Y[-149.69,-111.66],Z[-220.29,-169.38].
These ranges must not be used as a claim that either path exits the volume.

A useful next manual distinction would be holding an implicated pose still
versus passing through it, with pneumatic actuation still disabled. Persistent
loss while held supports a different mechanism from loss only during motion;
the current logs do not separate them. No such additional test was run here.

## Implication for controller work

The original tracker failure can occur independently of pneumatic actuation.
The separate residual-command Switch bug/NaN feedback cycle and the observer's
long re-warm still need their own fixes; neither explains this independent
Cygnus loss. Command-handoff smoothing can limit a loss's pressure consequence
but cannot recover missing transforms or justify weakening observation gates.

## Artifacts and checks

source_log.txt,summary.json,tracking.csv,missing_tip_reports.csv,
missing_episodes.csv,freshness_checks.json,pose_comparison.json;
manual_tracking.PNG/PDF/FIG,source snapshots and SHA256manifests.
Generic manual plotter now derives status counts from each capture instead
of hardcoding the earlier zero-loss result. No source model, runtime NDI,
controller or SLX changed. No device/serial/hardware/model/training access,
run,commit orpush by the agent. Offline saved-data analysis only.
