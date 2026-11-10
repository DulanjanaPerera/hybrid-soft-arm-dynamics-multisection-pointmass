# New Cygnus log: fast manual motion without pressure

2026-11-09. User replaced Desktop/log.txt with a NEW capture of manual/random
fast arm motion without pneumatic actuation. Treat separately from the previous
sinewave log. Current source2973106bytes,SHA256
7253f8388b1db67f0818371db0850e5691deef8385abd3b4354504b70d314290.
Its header says MonNov09 17:11:05CST2026; file lastwritten17:13:06.
The earlier P3 capture is preserved in ../cygnus_sine_20261109_v1/source_log.txt.

## Direct answer

**No runtime MISSING sensor transforms appear in this manual-motion recording.**
2448total BX replies all pass independent header/payload CRC and length checks.
Final tracking session2446replies covers61.4493s. Its first2reports are initial
acquisition with both tools missing; every2444subsequent pair is valid.
Earlier setup session has2base-missing/tip-disabled replies, also not movement
failures. These startup/setup states must not be counted as runtime dropout.

After acquisition:

The visible fast-motion burst is around46-54s. All320returned pairs in that
window are valid, with exactly8frame-counter ticks(~25ms) between every pair.
There are no query gaps or MISSING statuses during that movement window.

| Check | Base0A | Tip0B |
|---|---:|---:|
|Valid transform reports|2444|2444|
|MISSING reports|0|0|
|Repeated frame counters|0|0|
|Exactly identical consecutive8-value transforms|0|0|
|Out-of-volume flags|0|0|
|Partial-volume flags|0|0|

Base/tip frame counters match in every decoded pair and advance; systemstatus
is0throughout. Position/quaternion data keep updating. This provides no
evidence of frozen raw transforms during the observed movement. It does not
independently prove every GUI number refreshes each frame or that every
computed transform is physically accurate.

The user states the sensor stays inside the tracking volume; the status flags
are consistent with that. Neither this log nor the prior log establishes
out-of-volume as the cause of the previous missing transforms.

## Query gaps are different from MISSING status

The tracker reports40.033Hz transform updates and320.264Hz frame-counter
frequency. Normal frame increment8means24.979ms. Time here is counter elapsed
from the first reply of final tracking, not an absolute timestamp or GUI clock.
After first valid acquisition, four reply intervals exceed the usual increment:

| Between elapsed seconds | Interval | Frame increment |
|---|---:|---:|
|.0999-.1499|49.96ms|16|
|41.4658-41.5157|49.96ms|16|
|61.0996-61.2744|174.86ms|56|
|61.2744-61.4493|174.86ms|56|

The two long gaps are the last two reply intervals, just before tracking stops.
They could relate to application/communication/stop scheduling; this log has
no host timestamps to distinguish the exact cause. They are not evidence of
mid-run tracker transform failure. Intermediate frames were not queried or
recorded, so their validity is unknown; do not call them explicit missing
sensor reports or assume the hardware stopped generating measurements.

## Relation to the earlier pressure-driven P3 capture

The earlier sinewave capture had23explicit runtime tip-MISSING reports,
17episodes,in the inferredP3motion group. This new manual capture has0.
Fast manual movement alone did not reproduce that failure here. This does
not prove pressure/valve current caused it: poses, orientation and loading
were not matched between experiments.

Raw tip position coverage also differs. Current tipZ is[-340.87,-222.60]mm;
valid samples immediately adjacent to earlier P3 misses hadZ approximately
[-220.29,-169.38]mm(median-188.34). Base median rawposition shifts only
[.181,.043,.097]mm between captures. This supports different pose coverage,
not an assertion that either trajectory is outside the tracking volume.
Position-dependent orientation quality, motion/cable effects and actuation
interference remain possible but unproven. Do not loosen measurement validity
or infer a universal speed limit from these recordings.

Seeing continuous motion in Cygnus is compatible with brief invalid reports
in the earlier capture: most transforms still arrived normally. The log alone
cannot verify how that GUI displays an invalid sample. In this new capture,
however, there simply is no runtime MISSING record to explain a visible
measurement freeze. A GUI refresh/rounding/held downstream value would need
its own corresponding timestamps or recording to locate; do not assume it
is a sensor failure solely from the display.

## Retained evidence

source_log.txt,summary.json,tracking.csv,header-only missing CSVs,
query_gaps.csv,freshness_checks.json,manual_tracking.PNG/PDF/FIG,source snapshots
and SHA256manifest. Raw invalid startup data remain NaN; no interpolation.
analyzeCygnusBxLog.py now handles zero-runtime-missing captures with empty
event lists/header-only CSVs. plotCygnusManualAudit.m is separate from the
sine alignment plotter; no sine timing template was imposed on manual motion.

No serial/device/model access, runtime/controller/NDI edits, SLX changes,
hardware/acquisition/training run,commit orpush. Offline saved-data audit only.
