# Important decisions

Preserve decisions when reversed; add the replacement and reason. Dates below
are activity dates or documented stages, not newly invented approval dates.

## D-01 — Fixed-backbone reduced coordinates (clarified 2026-10-05)

Use two independent changes l2/l3 and l1=-l2-l3. The physical backbone constrains
the length sum. Equal ideal forces are a bending-force nullspace, not a free
axial coordinate. This supersedes the preceding discussion's suggestion of
appreciable common-mode shortening; common pressure may instead change prestress,
PMA bulging and incremental bending response. See technical notes.

## D-02 — Causal filtered-difference baseline (2026-10-01)

Retain the simpler filtered difference for the first controller; retain the
constant-velocity Kalman filter as comparison. Neither proxy establishes true
velocity accuracy. Increase analysis default from 0.12 to 0.20 s based on reduced
jitter; 0.30 s remains a comparison and 0.50 s requires delay assessment.
Skip invalid measurement updates and use pollStart timestamps as practical timing.

## D-03 — Coupled stiffness and bound expansion (2026-10-01–02)

Keep K coupled as k*[2 1;1 2], with a separate diagonal smooth bound penalty.
Raise modeled reduced bounds from +/-30 to +/-35 mm at user request; this is
not a measured physical limit. Historical inverse analysis retains +/-30 mm;
Inv_Dyn_val was write-blocked and is still an unverified pending saved-model item.
The subsequent 650 N/m sweep selection is conditional on its recording and
fixed assumptions, not a universal parameter replacement for every model.

## D-04 — Empirical phi-only stiffness is conditional (2026-10-05 backfill)

Use PCHIP for the measured P1 mean curve; standalone evaluation returns NaN
outside the measured range, while the paste-ready block explicitly holds endpoint
values. Poor unloading agreement and large variability prevent treating the
curve as an identified general mechanical law. Refit coefficients do not
automatically update the hardcoded block.

## D-05 — All-active common pressure before broad tracking (2026-10-05–06)

Test equal common levels, then unequal neutral trim and separate symmetric/
antisymmetric perturbations. The minimum-pressure allocator leaves a channel
inactive and cannot identify an all-active common-bias response. Keep its source
and historical results for comparison. Choose a 1.6 bar mean local candidate
rather than assuming 2 bar is best; equal-pressure bend and headroom matter.
The extrapolated [1.766,1.766,1.268] candidate overshot in validation;
manual repeated trial 4 [1.6656,1.6492,1.4852] became the local feedback baseline.
This replacement does not establish global straightness or absence of buckling.

## D-06 — Working pressure assumption (explicit 2026-10-05)

At the user's direction, treat commanded pressure as chamber pressure in
analysis because sensors are unavailable. Prefer post-limiter/post-Kill logs
when present. Preserve the distinction between this working assumption and an
independent pressure measurement; do not repeatedly block progress on sensing.

## D-07 — Bounded incremental feedback (2026-10-06)

Replace the old allocator with trim-aware, all-active mean-1.6 bar allocation.
Anchor only static G+Keff*q at 50 s; subtract it from the absolute model force,
then allocate incremental corrections with +/-0.04 bar s/a limits and normal
0.8 bar/s slew. Do not add the old active-channel 0.8 bar offset again.
The discussion proposed initially targeting measured pose; actual supplied test
instead targets geometric zero with a five-second feedback blend. Preserve this
explicit implementation change. The invalid-input slew exception remains open.

## D-08 — Small-step test before escalation (2026-10-06)

Use 0.01 rad U/V reference steps near trial 4, retaining Kill and local bounds.
Proceed with existing NDI and visual observation because no camera system exists.
The visible local buckle and history dependence limit claims of shape stability.
The one-off a+ point remains a candidate to repeat, not an approved trim replacement.

## D-09 — Forward predictor before NMPC (2026-10-05–06)

Consider future forward-model NMPC because direct learned inverse performed
poorly and pressure/history redundancy matters. Do not infer recursive capability
from measurement-assisted one-step validation. First audit timing/preprocessing,
test held-out simultaneous-actuation rollouts with predicted future histories,
and enforce the backbone constraint. Older Bouc-Wen mechanical-force parameters
are candidates only and do not remove pressure-allocation deadzone discontinuities.

## D-10 — Collaboration and documentation (2026-10-01, 05, 06)

Discussion is default; code/model edits and hardware work require explicit
requests. October 5 authorized ongoing log updates for new code or important
technical discussion. October 6 requested this seven-file documentation system:
current state edited in place; history append-only after initial migration;
specialized files maintained for material work. Read state plus relevant specialty
and search historical entries only as needed. This restructuring does not change
the technical-edit collaboration preference.

## D-11 — Circle reference settings remain in the block (2026-10-07)

For the user-pasted circle MATLAB Function, retain one time input `t` and
place desired phi, seconds per revolution, and integer circle count as
editable variables inside the block code. This supersedes the previously
suggested external `des_phi` and `circle_time_s` Simulink input ports. Each
lap smoothly runs through 2*pi while phi holds its target. The code remains
paste-ready guidance until the user installs and logs it; the saved model
still has the earlier one-input reference hardcoded to 0.1 rad.
