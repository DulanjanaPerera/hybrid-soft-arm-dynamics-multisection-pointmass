# R033: motion completed; score rejected for coverage

Source online_smoke_00g_phi1_T5_3cycles_r033.mat and worker_log.mat. Offline
review only. Three5s circles29.60-44.60s; withdrawal complete51.60s. Gamma
constant .0515/.0005999/1.49999, H scale/rate/upper5/10/50. No tracking limiting
or latched PAC faults. Peak pressures[2.6315,2.6191,2.6050]bar.

Last-two-circle exploratory metrics34.60-44.60s:

| Metric | Recorded | Limit |
|---|---:|---:|
| l2 RMSE | 4.681777 mm | 2 mm |
| l3 RMSE | 6.780124 mm | 2 mm |
| Vector velocity-error RMS | .01210224 m/s | .01 m/s |
| Fixed-raw roughness | .03081964 m/s^2 | .025 m/s^2 |
| Minimum sample/pair coverage | .8844221 | .95 |

All four metrics independently reconstructed from worker packets to<1e-8.
Inadequate coverage excludes them from GP training. '0 completed batches'
in old worker text meant zero accepted observations, not zero physical motion
batches: one motion batch completed. New worker distinguishes these counts.

Score window200samples:192rawvalid(96%),184controlled-valid(92%) after delayed
diagnostic masks,176/199consecutive pairs(88.442%). Eight tip-validity dropouts
at34.9,35.4,36.45,36.85,39.85,40.3,43.8,43.9s; basevalid1 andoverrun0. Exact
missing/stale/tracking cause not established. Raw-only pairs183/199=91.96%,
still below95%, so this is not solely delayed-mask conservatism. Coverage
policy was not loosened and invalid samples were not filled. No UDP packet
gaps in this saved stream. Inspect tip tracking/frame quality before search.

Per-lap controlled-valid length RMSE(mm):[5.3172,8.3457],[4.4699,7.1445],
[4.8845,6.3950]. l3 improves; l2 is not monotonic. Three circles do not prove
settled tracking or adaptation benefit.

Found rejection handshake bug: runtime requested2 but rejected score response
still carried1. Model waited30s and timed out reason11 at81.60. Phase10 at86.60,
commandszero86.90, stop87.90; all finalDAQzero. Corrected reply ID, explicit
rejection reasons and terminal handling without stationarity wait. Repeated
packets cannot mutate a terminal error. Rejected attempts are now retained
separately and checkpointed. Saved-data regression and synthetic tests pass;
no counterfactual physical performance claim and no new hardware run.

Artifacts: review.mat, phase/status transitions CSV, abort_neighborhood.csv,
lap_metrics.csv. Preserve this completed session. Further physical trial needs
a fresh session; no new experiment or gain changes queued by this review.
