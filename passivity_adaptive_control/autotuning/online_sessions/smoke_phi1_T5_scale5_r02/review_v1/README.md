# r02 smoke review: aborted before completing a circle

Recorded data: online_smoke_00g_phi1_T5_3cycles_r02.mat plus worker_log.mat.
Review created offline on 2026-10-28; no hardware/model execution.

Startup/trim/blend normal. Stationary hold20.55s, gain acknowledgement20.60s,
bend approach22.60s, angular spin-up27.60s, constant-speed circles29.60s.
At34.00s online reason4 caused terminal withdrawal; this is not expected
budget reason23. Only4.40s/0.88 circles of tracking, zero completed batches.
No completed scoring window and no optimization occurred.

All gamma values constant: .0515/.0005999/1.49999; hysteresis scale/rate/upper
5/10/50 verified. No tracking limiting and no latched PAC status5 faults.
PAC status4 occurred for isolated50ms samples and recovered to3.

Root cause is supervisor time alignment plus strict floating-point boundary:
raw measurement invalid33.90s -> valid33.95s -> invalid34.00s. At33.95s the
delayed PAC status/modelValid still reflect the preceding invalid frame. The
old controlled predicate therefore stayed false across all three samples,
joining separate dropouts. 34.00-33.90 is slightly above.10 in binary floating
point and triggered the strict >.10 comparison. Actual PAC recovered between
the raw dropouts and did not latch a fault.

Shutdown verified: phase10 at39.00s, command pressure below.01bar at39.30s,
simulation ended40.30s after1s zero-command dwell. Final commands and all three
DAQ outputs zero. Peaks[2.13823,2.42442,2.60675]bar. Physical pressure not sensed.

Early incomplete-run length RMSE[5.2206,8.9064]mm is exploratory only, not a
last-two-circle score, settled performance or evidence of instability. Adaptation
active90.91% in the short run. Nonbaseline gain changes and between-batch
state preservation were not exercised because no batch completed.

Corrected pacOnlineRuntimeStep: raw current validity controls dropout grace;
delayed diagnostics remain conservative score/hold eligibility and PAC's own
fault logic remains active. Added1e-9 comparison tolerance without increasing
configured grace. Invalid hold frames wait rather than being called movement;
fresh measured movement still aborts. Recorded prefix replay through34.20s
no longer falsely aborts; sustained invalidity and actual movement tests still
abort. This replay is not proof of counterfactual later tracking performance.
Fresh r03 session/model prepared with identical fixed gains and no search.

Artifacts: review.mat, phase_transitions.csv, status_transitions.csv,
abort_neighborhood.csv. Reproduce with reviewPacOnlineSmoke and a new output
folder; validatePacOnlineValidity tests the saved-prefix regression.
