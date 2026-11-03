# November2 evaluation of frozen generalized-force GRU

2026-11-03. Canonical model: ../release_20261103_v1/force_model.mat.
Full model/training/numerical description: ../release_20261103_v1/REPORT.md.
Release SHA256 b5cfd7d0cd13b2b9067a35d810b8ae166e032e96421b6a30c55d703039af268f
was checked unchanged during evaluation. No fitting, normalization refit,
candidate selection, Simulink or hardware execution in this evaluator.

Five November2 recordings supply50,078 guarded one-step anchors. The same
exported data were already audited; their earlier inspection is acknowledged.
All targets are inferred lengths l2/l3, not independent strain or force truth.

| Predictor | One-step macro vector RMSE,mm | .5s path | 2s path | 5s path |
|---|---:|---:|---:|---:|
| Corrected nonlinear nominal, zero residual | 4.7352 | 6.9508 | 7.0438 | 6.6372 |
| Constant velocity | .06235 | .18575 | .89271 | 3.4062 |
| Force-residual GRU | 2.1510 | 1.9373 | 1.8375 | 2.1267 |

Sixty matched initial windows per model/horizon,12/run. Each begins from
observed history and then uses predicted q/integrated v only. Actual recorded
pressure/time sequences are known exogenous inputs; command alignment uses
the declared NDI-row ZOH approximation, not measured chamber pressure.
No future measured response is fed back. Across horizons, eligible start
pools can differ. Run path metric is sqrt(mean window pathMSE); macro gives
each of five experiments equal weight. Windows overlap and are not independent
experiments. Individual scores, durations and indices are retained in CSV.

All540 model/window propagations complete without the .05m independent-
coordinate numerical guard; zero one-step failures. Guard status is not a
physical safety, passivity or control-stability guarantee. Successful-only
metric columns preserve explicit failure counts even though all arezero here.

Conclusion: the force GRU improves the matched nominal physics predictor but
has weak short-horizon accuracy. Its5s macro error is37.6% below constant
velocity; its one-step/.5s/2s errors are worse. No controller/inverse readiness,
unique physical-force identification or best-architecture claim. This model
was NOT refit after observing these results. New force outputs remain in N,
inside the mechanical equations, unlike earlier position-correction predictors.

Retained: one_step_by_run.csv, per-recording MAT predictions/forces,
rollout_by_window.csv, rollout_windows.csv, rollout_summary.csv,
rollout_predictions.mat with every q/truth/p/residual/time path,
residual_output_quality.csv, heldout_sources.csv and protocol.mat.
force_model_evaluation.fig/pdf/png contains the old-only selection curve and
held-out one-step/2s/5s comparisons. Source evaluateForceResidualModel.m retained.
Paper caption must state prediction in inferred lengths, old-day training,
prior-inspected held-out date, approximate command timing and no controller test.
