# Length-policy baseline diagnostic replay

Existing user-acquired recordings reevaluated with pac-tune-2. Raw MAT files are
unchanged. Configuration required limits remain NaN; all results are diagnostic,
not accepted optimizer candidates. See ../analysis_review_v1/README.md for the
full experiment assessment, settings and initialization caveats.

Vector length RMSE remains the objective, with separate l2/l3 acceptance limits.
Angular limits are optional. Required smoothness uses vector velocity-error RMS
and RMS of velocity-error slopes over consecutive valid ungapped samples.

Late baseline 1/2 velocity RMS: .0030014/.0029658 m/s. Error-rate RMS:
.019338/.019177 m/s^2. Pair coverage:97.33/98.33%. These values characterize the
baseline; no numeric threshold has been selected from them. Length RMSE and
tracking results match the earlier review. Reproduce using reviewFirstTrials
with a fresh output folder; keep the same filtering and sample rate.
