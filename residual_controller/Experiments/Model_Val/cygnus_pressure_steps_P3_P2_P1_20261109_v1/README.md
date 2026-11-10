# Individual pressure-step capture: P3 first, then P2, then P1

2026-11-09. NEW Desktop/log.txt, user reports manually applying multiple
pressure steps to each PMA with significant pauses between PMAs, in order
**P3 ->P2 ->P1**. Source5986190bytes,SHA256
1976de55d612ba40bad58af96c84073e1eb87b7a6b466295563a1a46cbf10e87.
HeaderMonNov09 17:59:50CST2026; file lastwrite18:02:55. Archived source_log.txt
keeps it separate from all earlier captures that shared the Desktop filename.

## Findings

All4978BXreplies pass header/payload CRC and length/layout checks. Final
tracking session4976replies covers125.496s.3initial acquisition reports
precede the first jointly valid reading and are excluded from runtime counts.
Runtime4973:base0Avalid4973/missing0,tip0Bvalid4941/missing32.

The three response groups are clear from raw sensor motion and the long
inter-group pauses. Labels use the user-confirmed activation order. No actual
pressure command amplitudes or synchronized pressure-event timestamps exist
in this tracker log; windows below are observed response windows, not exact
pressure-on/off times.

| PMA/order | Approximate visible response window(s) | Runtime tip MISSING |
|---|---|---:|
|P3 first|8.97-37.32|32|
|P2 second|57.10-83.81|0|
|P1 last|96.25-121.22|0|

All32missing reports occur9.0925-37.2443s, within theP3group. They form
30consecutive-missing episodes:28single and2double reported updates, lasting
approximately25-50ms to the next valid reported sample. Every missing tip
report has a valid base transform and matching frame counters. No out-of-
volume/partial-volume bits are set; systemstatus0throughout. No flags here
establish out-of-volume as a cause. Device status0x02explicitly means the
tip transform is unavailable; these are not just gaps in host queries.

Tracker reports40.033Hz transform updates and320.264Hz frame counter rate.
Elapsed seconds=(frame-firstFinalSessionFrame)/320.264. Normal increments8;
larger query increments16/40/48are retained separately and do not constitute
known MISSING transforms for any unqueried intermediate slots.

## Transitions rather than held step plateaus

Losses repeatedly accompany P3 loading/unloading response transitions. Held
plateaus return valid transforms. Across each gap, the immediately preceding
and following valid poses change substantially; their position-difference
motion indicator spans181-741mm/s(median380). This is not ground-truth
instantaneous velocity and does not identify a speed threshold or reconstruct
the missing pose.

P2/P1 transitions also move quickly but have no explicit missing transform in
this capture. Their valid-pair raw-position rate maxima are approximately
610/610mm/s, overlapping the P3 missing-interval indicators. Sampling and
pose/derivative differences prevent identifying a universal speed limit.
The evidence supports a direction/pose/orientation-dependent issue, potentially
including motion/cable effects, but does not uniquely diagnose it.

## Interpretation across experiments

P3 still produces the losses when tested FIRST, while laterP2/P1 remainvalid.
Therefore the association is not explained merely by being the last channel
or by longer elapsed experiment time. Earlier standard sine also lost the
tip in inferredP3groups, and earlier manualP3-direction movement WITHOUT
valve actuation reproduced30missing reports. Together these support a P3-
direction tracking issue; valve actuation is not necessary. They do not prove
valve/electronics interference can never add an effect in other conditions.

The remaining useful distinction is an intermediate problem pose held still
versus moving through it. The current step plateaus are at the endpoints,
so valid endpoint holds do not prove the intermediate poses track when held.
No new manual/hardware experiment was performed by the agent.

## Group classification and artifacts

Classification uses valid raw-tip displacement>30mm from an initial1-5s
baseline. Active samples separated by at most8s are joined into a motion
group; long rest-gap midpoints47.211s/90.026s separate the3groups. No sensor
values were filled or interpolated. Counts per broad group including rests
are1868/1699/1406runtime reports and32/0/0missing, summing4973/32.
The PMA labels depend on user-reported order, not a pressure field in the log.

source_log.txt,summary.json,tracking.csv,missing_tip_reports.csv,
missing_episodes.csv,movement_groups.csv,movement_classification.json;
annotated manual_tracking.PNG/PDF/FIG,source snapshots and SHA256manifests.
The plotter retains a legacy output filename but the figure title explicitly
identifies this pressure-step experiment. Parser checks, event counts and
per-PMA totals are retained. No runtime/controller/NDI/SLX edit, serial/device
access,model/hardware/trainingrun,commit orpush. Offline log analysis only.
