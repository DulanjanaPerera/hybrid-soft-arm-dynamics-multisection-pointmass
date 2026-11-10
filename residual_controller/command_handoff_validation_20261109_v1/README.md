# Command handoff offline verification

2026-11-09. New residualCommandValid.m is a7input/4output stateful pressure
handoff. Full wiring and mode definitions: ../COMMAND_HANDOFF.md.

validateResidualCommandHandoff.m verified:

-Nominal startup, ramped residual entry, brief invalid hold, timed fallback
 fade and re-entry after recovery; no invalid candidate arithmetic.
-0..3bar rails and shared.30bar-per-.05s ordinary per-channel step bound.
-Repeated same-timestamp evaluations do not advance timers or state.
-Disable wins on the same timestamp; reset and bad inputs are handled.
-Long/fast call intervals restart at actual pressure without catch-up steps.
-Limiter anchors to actual applied feedback after Kill, not stale cache.

7501recorded test5 candidate/validity/nominal samples were then passed through
the new function with an explicit simulated one-tick command delay. The
candidate and validity inputs were frozen from the old physical recording;
no physical arm, sensor, GRU state or forward dynamics was simulated under
the changed commands. This tests pressure selection only, not new tracking.

Result: all finite outputs; maximum ordinary new command step.3000000000000003
bar,0steps above.3+numeric tolerance. The old recorded command had30such steps,
maximum1.513951bar.74grace-hold ticks,370fallback-fade ticks,339re-entry ticks.
Current source explicitly requires a fixed.05s tick and rejects inconsistent
call intervals. Disable/reset/input-fault zeroing is an intentional exception
to ordinary slew smoothing; physical Kill remains independent downstream.

Retained summary.json,verification.mat,frozen_input_commands.csv,
command_handoff_example.PNG/PDF/FIG,source snapshots and SHA256manifest.
No SLX edited/loaded/compiled/run, no hardware/acquisition/pressure command,
no GRU/model/training change, no commit or push. Actual operator model build,
real-time timing and physical behavior remain unverified for this new source.
