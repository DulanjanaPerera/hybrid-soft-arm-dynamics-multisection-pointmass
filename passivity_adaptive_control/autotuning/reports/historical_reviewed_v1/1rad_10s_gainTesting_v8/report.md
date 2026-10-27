# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\adaptCommonP\1rad_10s_gainTesting_v8.mat`

SHA-256: `a7fa378a622bfcd5e6047a5667d0a26df3fad0899902c58c3bfd13ef0bdea344`

Classification: **protocol_deviation**. Success: 0. Schema: pac-tune-1.

- Targets are unset; historical diagnostics cannot certify success.
- Live tuning changes detected; exploratory trial.

Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 290
end_s | 320 | 320
samples | 6000 | 600
coverage | 0.9896667 | 0.985
thetaCoverage | 0.9896667 | 0.985
q2Rmse_mm | 2.024831 | 0.8625407
q3Rmse_mm | 2.30449 | 1.396839
qVectorRmse_mm | 3.067672 | 1.641686
thetaMae_rad | 0.1786755 | 0.09942834
phiMae_rad | 0.07635119 | 0.06006832
phiMax_rad | 1.26505 | 1.15769
velocityRmse_mPerS | 0.003387314 | 0.003065404
thetaRateStd_radPerS | 0.3464339 | 0.3746036
nearStopFraction | 0.04 | 0.0395189
slowOrReverseFraction | 0.04034043 | 0.0395189
activeFraction | 0.9615 | 0.9216667
limitedFraction | 0 | 0
faultSamples | 0 | 0

Tuning changes: 2. Replay: 5938 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.98967.

Observed parameter rate-cap fractions: `[0 0 0.0043335 0 0 0 0 0 0 0]` (active transitions: 5769).

Minimum estimated K eigenvalue: -347.0622 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.52329 1.91353 1.36078]`.

Tracking parameter change: `[-391.844 -576.317 -401.16 10.8649 14.6492 11.9587 11.3969 11.5066 -1.11969 0.496887]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
