# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\candidate04_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `ff1090d4abf5232a02ee1dcd9c3d6112b5845bb30c862f7bc7d7660626bd02e0`

Classification: **valid_infeasible**. Success: 0. Schema: pac-tune-2.

Campaign classification: **valid_infeasible**. Eligible for candidate selection: **0**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.99125 | 0.9883333
thetaCoverage | 0.99125 | 0.9883333
phiCoverage | 0.99125 | 0.9883333
velocityCoverage | 0.99125 | 0.9883333
velocityPairCoverage | 0.9820833 | 0.9766667
q2Rmse_mm | 2.98754 | 1.692167
q3Rmse_mm | 3.422139 | 1.940342
qVectorRmse_mm | 4.542734 | 2.574559
thetaMae_rad | 0.3090862 | 0.1669633
phiMae_rad | 0.1071276 | 0.06272322
phiMax_rad | 1.282075 | 1.106293
velocityRmse_mPerS | 0.003819606 | 0.003027893
velocityErrorRateRmse_mPerS2 | 0.019019 | 0.02072852
thetaRateStd_radPerS | 0.3249554 | 0.3444063
nearStopFraction | 0.04412389 | 0.05290102
slowOrReverseFraction | 0.04412389 | 0.05290102
activeFraction | 0.9825 | 0.9766667
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 1 | 1.692167 | exceeded/unavailable
q3Rmse_mm | 1 | 1 | 1.940342 | exceeded/unavailable
velocityRmse_mPerS | 1 | 0.01 | 0.003027893 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.02072852 | within limit
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 2.574559 | disabled
thetaMae_rad | 0 | NaN | 0.1669633 | disabled
phiMae_rad | 0 | NaN | 0.06272322 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.3444063 | disabled

Tuning changes: 0. Replay: 2379 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.99125.

Observed parameter rate-cap fractions: `[0 0 0.0059372 0 0 0 0 0 0 0]` (active transitions: 2358).

Minimum estimated K eigenvalue: -412.7168 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[0.994376 1.16532 1.13897]`.

Tracking parameter change: `[-204.875 -723.009 -214.236 13.8338 15.1942 9.52449 8.47983 8.77226 -1.00737 -0.394539]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
