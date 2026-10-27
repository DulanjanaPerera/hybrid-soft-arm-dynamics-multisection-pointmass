# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\adaptCommonP\1rad_10s_gainTesting_v9.mat`

SHA-256: `d3b7e0d0ebe6896535728762ad33707af2b041caddf3b7664ebafceeb32c8756`

Classification: **invalid_data**. Success: 0. Schema: pac-tune-1.

- Targets are unset; historical diagnostics cannot certify success.
- Coverage below configured minimum of 0.95.
- Live tuning changes detected; exploratory trial.

Reference period: 2 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 74
end_s | 80 | 80
samples | 1200 | 120
coverage | 0.9483333 | 0.9583333
thetaCoverage | 0.9483333 | 0.9583333
q2Rmse_mm | 7.498189 | 8.082374
q3Rmse_mm | 9.838568 | 9.648521
qVectorRmse_mm | 12.37014 | 12.58645
thetaMae_rad | 1.081054 | 1.076262
phiMae_rad | 0.293333 | 0.2815402
phiMax_rad | 1.01844 | 0.9704639
velocityRmse_mPerS | 0.04133608 | 0.04182448
thetaRateStd_radPerS | 2.010039 | 1.38885
nearStopFraction | 0.01751152 | 0
slowOrReverseFraction | 0.02119816 | 0
activeFraction | 0.05416667 | 0.03333333
limitedFraction | 0.8925 | 0.925
faultSamples | 0 | 0

Tuning changes: 2. Replay: 1138 samples, 0 allocation limits, 1071 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.94833.

Observed parameter rate-cap fractions: `[0 0.015385 0.10769 0.015385 0 0 0 0 0.015385 0]` (active transitions: 65).

Minimum estimated K eigenvalue: 65.3461 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[2 2 2]`.

Tracking parameter change: `[856.2 -585.653 -649 12.3402 2.34559 0 0.918252 0 -0.00142857 0.0317072]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
