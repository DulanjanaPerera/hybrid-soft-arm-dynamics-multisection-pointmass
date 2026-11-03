# November2 circle, sine and helix: frozen force-model prediction plots

Frozen GRU64 seed113 selected on October19/27, .5s training horizon. Model hash
8150352dd47230a2a4094b240e1e886b641b8c22ebad0cbfe300f63e9f526c6a.
No weights, normalizer, controller or recorded measurement values were changed.
November2 is previously inspected evaluation; these are additional long-horizon
diagnostics, not a fresh blinded confirmation or new model-selection campaign.

## Plots

- circle_free_prediction: one full requested-phi1 circle,259-279s simulation.
- sine_free_prediction: one period for each pressure channel, rows1/2/3.
- helix_free_prediction: first loading/unloading episode after initial warmup,
 6.6-124.95s simulation (119.85physical seconds); first1.6s of the reference
 episode is observed-history initialization, not part of the free prediction.
- circle_one_prediction, sine_one_prediction, helix_one_prediction: EXACT saved
 one-step predictions from the frozen streaming evaluation, shown separately.
- nominal_free_prediction: poor nominal predictions retained in their own figure.
Each figure has editable FIG, vector-requested PDF and300dpi PNG. Underlying
predictions/observations/pressures/physical time/masks are in predictions.mat;
scores.csv and selected_windows.csv give definitions and actual selections.

Free predictions receive measured q/v and the observed recurrent prefix ONCE,
followed only by integrated q/v, recorded commands and timestamps. No hidden-
state or length resets occur during these windows. Future pressure inputs are
known recorded commands; this is not unknown-command forecasting.
The long windows exceed the .5s training and previous1s primary test horizons.
All15 model/window free predictions completed without the inherited numerical
guard stopping them. This guard is not a controller stability guarantee.

## Full-window free prediction scores

Length-vector RMSE=sqrt(mean(e_l2^2+e_l3^2)), millimetres, over real valid
measurements. Windows selected chronologically by declared input/time criteria,
not minimum model error. These five windows are examples, not independent repeats.

| Window | Force GRU64 (mm) | Constant velocity (mm) | Nominal (mm) |
|---|---:|---:|---:|
| circleSweep | 0.9440 | 2.2530 | 11.4814 |
| SineWave4 channel1 | 2.3355 | 9.4811 | 6.0698 |
| SineWave4 channel2 | 1.7784 | 17.4802 | 9.6563 |
| SineWave4 channel3 | 1.7958 | 21.7312 | 7.5594 |
| helix | 0.7627 | 4.7396 | 8.9723 |

The force GRU is closer than CV and nominal on all five selected long windows.
This coexists with CV winning the previous .1-1s tests and the one-step panels:
long recorded-input propagation and measurement-conditioned one-step prediction
are different questions. Do not claim universal superiority from these examples.

## Measurement gaps and geometry

Future NDI loss does not by itself require terminating a free rollout: known
pressures/time still drive the model. Three missing labels in sine channel3
and three in the helix stay NaN, unscored and visibly marked by gray strips.
No sensor value was interpolated or held as truth. Time/pressure/trial guards
are enforced throughout. The one-step panels retain the original gap reset
and recovery warmup, so they have additional unavailable prediction intervals.
Gray strips identify missing observations, not every warmup omission.

PCC X/Y/Z are reconstructed from inferred/predicted lengths with r=.013m and
L=.17411m. They are not independent sensor XYZ measurements. Requested theta/
phi are gray dotted references, separate from measured arm response. The
circle requested phi1rad but actual phi is much smaller, as in earlier reports.
These open-loop recorded-pressure experiments are model-response comparisons,
not a newly run PAC/PDFB tracking test. The reference named helix is the recorded
outward/inward theta/phi spiral, not an exact constant-pitch Cartesian helix.
Theta plots unwrap finite observable spans and lift by integer2pi only; near-
straight/unobservable direction is blank. Theta error is circular, evaluated
only where observed and predicted directions are meaningful.

## Verification and November3 companion audit

All15 first free transitions reproduce the exact frozen saved one-step outputs
within1e-10m. Measured phi reconstructed from q agrees with logged orientation.
Raw/prepared/export/model hashes remain unchanged; previous force release kept.
The plotting source and numerical propagation dependencies are archived/hashed.
No training, output smoothing, raw-data filling, SLX/hardware, commit or push.
November3 exports and the separate53-sample interpolation PREVIEW are under
../motion_predictions_20261103_v1/short_gap_audit/REPORT.md. They have not
entered any training or frozen-evaluation pipeline.
