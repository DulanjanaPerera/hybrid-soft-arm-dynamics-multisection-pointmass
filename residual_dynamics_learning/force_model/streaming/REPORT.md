# Streaming generalized-force residual: completed experiment sequence

2026-11-03. Eleven completed offline fits: six architecture/horizon screen fits,
one burn-in comparison and four additional finalist-seed fits. October19/27
only for training and selection. All finalist weights and the representative
research checkpoint were frozen before this November2 replay. The previously
inspected November2 data are evaluation only, not blind confirmation.

The selected configuration is **B: 64 units, 
0.5s maximum physical training rollout**. Selection uses
the lowest mean older-holdout score across three seeds; seed113 was declared
in advance as the representative, rather than picking the luckiest seed.
Research checkpoint: `release_20261103_v1/force_model.mat` in this folder.
SHA256: `8150352dd47230a2a4094b240e1e886b641b8c22ebad0cbfe300f63e9f526c6a`. Previous parent release is unchanged.

**Main limitation:** the force GRU substantially improves the nominal model,
but constant velocity remains better over all tested horizons through1s.
The numerical inverse works on self-generated targets; this does not yet
justify physical pressure correction or controller deployment.

## Contract and checks

The GRU updates each layer once from current measured q, causal estimated v,
current pressure command and previous physical dt. Its outputs are two
generalized forces in N on the mechanical RHS; lengths follow from the same
fixed nonlinear nominal dynamics and SDIRK integrator. Chronological .5s
training chunks carry measured-input hidden states, detached between updates.
Auxiliary .5/1s free prediction starts with a separate pre-current state copy,
uses integrated q/v and recorded future commands, and never commits its state.
It is a known-command prediction test, not unknown-command forecasting.

Original whole-trial splits and TRAIN-only normalizer preserved exactly. 5533
train chunks, 2087 selection chunks in 241 valid mature segments;
128 common selection anchors across eight eligible recordings. Day2 StepWave
has two train chunks and no holdout coverage; day2 combined has neither.
All anchors obey actual poll-start time caps and original validity/gap guards.

Numerical/state checks pass for 1x32,1x64,2x32: stable gates, finite gradients,
gradient versus central differences, chronological state replay, masked carry,
future-target independence and auxiliary-state isolation. Full anchor time/
segment/partition audit passes. Saved one-step metrics independently recomputed
within1e-8mm; old release, normalizer, source/data/model hashes unchanged.

## Six-model screen

Old-only score = .7 macro one-step vector RMSE + .3 macro .5s free-path vector
RMSE. Vector means sqrt(e2^2+e3^2) before time RMS, then equal recording means.
This is NOT a single-coordinate or closed-loop tracking error.

| ID | GRU units by layer | Horizon (s) | Parameters | Old score (mm) |
|---|---|---:|---:|---:|
| A | 32 | 0.5 | 4018 | 1.1661 |
| B | 64 | 0.5 | 14162 | 0.5441 |
| C | [32 32] | 0.5 | 10258 | 0.6568 |
| D | 32 | 1 | 4018 | 2.0105 |
| E | 64 | 1 | 14162 | 0.6183 |
| F | [32 32] | 1 | 10258 | 0.7854 |

All fits used seed113, batch32, Adam .0003, gradient clip1 and maximum1200
updates with the same early-stop policy, checkpoint cadence, data, normalizer,
force bound and loss weights. Length terms separately average one-step and
free-path errors (.5/.5), plus .001 force-rate and .0001 magnitude penalties.
Two layers did not win this screen. This is a practical capacity comparison,
not parameter-matched proof that depth is inferior. Learning curves were still
improving near the budget cutoff; rankings are budget-specific.

Leader burn-in .5s: 0.54411mm; 1s: 
0.56325mm. The worse 1s variant was not retained.

| Finalist | Seed113 old score | Seed211 | Seed307 | Mean +/- SD (mm) |
|---|---:|---:|---:|---:|
| B | 0.5441 | 0.5299 | 0.5830 | 0.5523 +/- 0.0275 |
| E | 0.6183 | 0.6199 | 0.6314 | 0.6232 +/- 0.0071 |

## November2 frozen evaluation

47,953 common one-step anchors per model/seed; 60 matched free windows
per model/seed, each scored within .1/.25/.5/1s physical caps. Primary replay
refreshes q/v each sample; free windows have no measured-state resets.
Finalist numerical guard failures: 0. No failed paths silently removed.

| Predictor | One-step vector RMSE (mm) | .1s path | .25s path | .5s path | 1s path |
|---|---:|---:|---:|---:|---:|
| B, three seeds mean +/- SD | 0.5629 +/- 0.0684 | 0.5026 +/- 0.0559 | 0.3777 +/- 0.0209 | 0.3613 +/- 0.0148 | 0.4272 +/- 0.0183 |
| E, three seeds mean +/- SD | 0.6391 +/- 0.0675 | 0.5628 +/- 0.0540 | 0.4688 +/- 0.0263 | 0.4357 +/- 0.0191 | 0.4777 +/- 0.0292 |
| nominal | 4.7638 | 4.5842 | 6.5245 | 6.9113 | 7.0445 |
| constant_velocity | 0.0621 | 0.0212 | 0.0522 | 0.1138 | 0.3315 |

The released seed113 representative has November2 one-step macro RMSE
0.5359mm and .5/1s path RMSE
0.3722/
0.4379mm. The table above reports
the full seed distribution, not only that representative.

Mean/spread are across three fixed seeds of equal-run macro scores, not
confidence intervals or independent hardware repeats. Per-coordinate RMSE,
p95 absolute error, waveform phase and measured-speed groups are in CSVs.
Speed groups are norm(v)<=.005, .005-.015 and >.015m/s; these are measured
motion groups, not new cycle-speed experiments.

Max residual magnitude across finalist replays 19.2616N; maximum
per-run occupancy within1% of the80N head bound 0.0000%;
largest per-run force-rate p95 5.6670N/s. Force quality and all spikes
are retained without post-hoc filtering. Rates use adjacent samples within
valid segments; reset boundaries remain explicit in saved segment labels.
These are effective compensation
forces, not measurements or uniquely identified physical parameters.

The parent release used a different warm-up/free-prediction contract and
anchor set. Its previous 2.151mm one-step result cannot isolate a width/depth
effect against this streaming experiment. No post-November2 tuning occurred.

## Offline pressure inversion

Best mean OLD-holdout configuration, fixed seed113. A known logged pressure
generates its own target through B*p+r(q,v,p,h_previous). The solve holds
common pressure at the logged mean, with0..3bar bounds and0.3bar step cap
(6bar/s at intended.05s). The0..3bar range is the offline RL command domain,
not a changed PAC setting. A deliberately displaced starting guess avoids
trivial solves on holds. Trial optimizer calls reuse immutable pre-current
hidden states; only the selected pressure commits one new state.

128/128 logged targets satisfy those prescribed bounds;
excluded targets are explicitly retained as infeasible, not counted successes.
Max force-equation residual 4.7431e-05N; max bound violation 2.2204e-16bar.
Minimum observed differential pressure Jacobian singular value 9.7578N/bar.
Solve+selected-GRU latency median 5.190ms; p95
13.927ms; 2/128 exceed50ms.
Maximum latency 405.837ms includes the first optimizer call;
optimizer cold-start cost was retained, not excluded to improve the timing.
Exact MATLAB summary and per-case timings/exit flags/continuity gains are
retained. Latency excludes acquisition, nominal-control computation and
real-time scheduling; this is NOT qualification for the50ms hardware cycle.

Self-consistent inversion cannot establish true pressure correction or
tracking improvement. Commands are pressure proxies; q is inferred from
orientation. Missing measured force/chamber pressure and current-command
dependence/actuator lag remain material limits. No PAC/PDFB/SLX integration
or arm execution was performed.

## Retained artifacts

- campaign_20261103_v2: predeclared protocol, all checkpoints, source hashes,
  split coverage, screen/history/seed CSVs and all six frozen finalist weights.
- release_20261103_v1: selected research checkpoint, training source and hash.
- evaluation_20261103_v1: all predictions/forces/truth/commands/time in MAT,
  per-run/phase/speed and baseline CSVs, source hashes, independent verification.
- inverse_20261103_v1: targets/pressures/states, cases CSV, timing and protocol.
- Three FIG/vectorPDF/300dpiPNG sets: streaming_comparison,
  circle_streaming_example and circle_short_free_example. Example windows are
  first chronologically selected, with prediction mode explicit in captions.
- Failed campaign_v1 retained for reporting-error provenance; no model released.

New November3 recordings were not pooled. File migration, hardware testing,
code/model deployment, Git commit and push remain unperformed.
