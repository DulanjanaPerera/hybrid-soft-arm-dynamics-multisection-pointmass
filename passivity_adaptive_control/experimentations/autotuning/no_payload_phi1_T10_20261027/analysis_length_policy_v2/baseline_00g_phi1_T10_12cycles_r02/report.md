# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\baseline_00g_phi1_T10_12cycles_r02.mat`

SHA-256: `422f88a9ad89757bab0a58a21c18d5c481376154e9d4ea7dff89726785b7e326`

Classification: **unconfigured_targets**. Success: 0. Schema: pac-tune-2.

- Targets are unset; historical diagnostics cannot certify success.

Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.9883333 | 0.9916667
thetaCoverage | 0.9883333 | 0.9916667
phiCoverage | 0.9883333 | 0.9916667
velocityCoverage | 0.9883333 | 0.9916667
velocityPairCoverage | 0.9766667 | 0.9833333
q2Rmse_mm | 2.969287 | 1.528962
q3Rmse_mm | 3.413898 | 1.856881
qVectorRmse_mm | 4.524529 | 2.405355
thetaMae_rad | 0.3007121 | 0.1483596
phiMae_rad | 0.1066997 | 0.06234168
phiMax_rad | 1.246135 | 1.104598
velocityRmse_mPerS | 0.003792854 | 0.00296576
velocityErrorRateRmse_mPerS2 | 0.01828605 | 0.01917664
thetaRateStd_radPerS | 0.3182145 | 0.3268788
nearStopFraction | 0.03540956 | 0.04745763
slowOrReverseFraction | 0.03583618 | 0.04745763
activeFraction | 0.9766667 | 0.9833333
limitedFraction | 0 | 0
faultSamples | 0 | 0

Tuning changes: 0. Replay: 2372 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.98833.

Observed parameter rate-cap fractions: `[0 0 0.0017065 0 0 0 0 0 0 0]` (active transitions: 2344).

Minimum estimated K eigenvalue: -323.2088 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.03917 1.0645 1.11397]`.

Tracking parameter change: `[-251.617 -707.212 -240.661 7.31034 8.04778 10.0542 9.20486 9.48719 -1.11846 -0.167142]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
