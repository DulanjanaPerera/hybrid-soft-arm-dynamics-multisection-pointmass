# Continuous tuning integration (offline prototype)

**Update 2026-10-28:** The new fixed-rate runtime, separate worker, generated
library and wired model copy are documented in `ONLINE_SIMULINK_WORKFLOW.md`.
The text below preserves the earlier host-only prototype contract; its missing
adapters statement does not describe the new implementation. The new model
copy still requires operator-run hardware integration verification.

`pacTuneOnlineConfig` and `pacTuneOnlineStep` implement a host-side batch state
machine around the existing BO proposer. They do not run Simulink, apply gains,
collect scores or operate devices. Do not paste the host function into a MATLAB
Function block: it uses tables and a host-side optimizer. The SLX integration
adapter and independently operating hardware watchdog are not implemented.

The user requested unattended tuning on 2026-10-28. Automatic operation here
means no approval between candidates after the operator starts a configured
session, not autonomous assistant operation of the arm.

## Contract

Call `[state,request]=pacTuneOnlineStep(state,telemetry,config)` with empty state
only to create a new session. Persist state and append raw batch logs outside
the real-time task. An unexpected host restart must stop the tuning session;
resumption needs reconciled controller/session IDs, not an automatic replay.
Bounds and baseline gamma are deliberately NaN until configured. Do not reuse
scale-1 campaign observations with scale-5 settings: the gradient multiplier is
25 times larger at equal gamma. Archive the complete current tuning snapshot,
filters, period and physical conditions before starting a separate campaign.
Current user-supplied hysteresis scale/rate/upper settings are 5/10/50; actual
gammaK/D/H and approved new search bounds remain to be selected.

States: hold -> apply -> settle -> run -> score -> hold, or terminal stopped.
The request includes mode, candidate ID, gamma, apply flag and explicit
preserveAdaptiveState=true. The adapter must atomically apply all three gamma
groups in a verified stationary hold and acknowledge the exact ID and values.
Do not reset PAC weights, history or filters. `holdReady` must mean the actual
trajectory ramp has completed and measured motion meets the chosen stationary
criterion; a requested hold alone is insufficient. Run requests must be handled
once per candidate ID, with deterministic trajectory start timing. Timing of
the host run request and actual reference start must be checked by the adapter;
this prototype's timer alone is not sufficient evidence of a completed batch.

Telemetry fields are documented in pacTuneOnlineStep. `healthy` must combine
latched fault, measurement validity/freshness, controller status, command limits,
transport heartbeat and external enable. Its host-side handling is supplementary:
a separate real-time watchdog must enforce stop/hold if the host hangs. A stop
request must invoke the existing supervised stop/ramp policy; do not abruptly
zero pressures solely because optimization ended. Emergency kill remains
independent. Optimizer elapsed-time checks cannot interrupt a hung optimizer.

`scoreReady` is an audited completed-batch result, not an arbitrary recent
window. `metrics` order is l2/l3 RMSE in mm, vector velocity-error RMS in m/s,
and fixed-raw velocity-error-rate RMS in m/s^2. Use the existing evaluator with
actual start/stop times, final scoreCircles window, full-batch fault/limit checks
and coverage checks. Missing, malformed, mismatched tuning or incomplete logs
must never be marked ready with invented scores. No evaluator adapter is yet
provided. Healthy must latch transient hazards across the full batch.

## Search behavior and limitations

The first candidate is the configured baseline. Infeasible accuracy results
enter BO; pressure limiting, unhealthy telemetry, invalid score and timeouts
stop the session. Local proposals are restricted to a multiplicative factor
(default 1.5) around the current gains within global bounds. This is a step-size
restriction, not a proven safe set. The BO feasibility score is heuristic.
Observations outside the local region are excluded from that proposal's GP;
the full session ledger retains them and determines the global budget.

Two consecutive passing batches at identical gains stop successfully at the
operational target; these are not independent reset confirmations or a
stability guarantee. Preserved learning makes the objective history-dependent.
The prototype uses stationary BO heuristically, without a contextual GP or
causal attribution. It does not automatically schedule speeds or expand bounds.
No automatic rollback is implemented: a failure terminates to supervised stop.

## Verification before hardware integration

Run `validatePacTuneOnline` offline. Next implement and test the model adapter,
nonblocking transport, heartbeat loss handling, full raw-log evaluator bridge,
gain transition and candidate-ID acknowledgement, actual trajectory timing and
restart behavior. Use operator-run tests before enabling unattended batches.
No SLX has been edited, model run or gains applied by these files.
