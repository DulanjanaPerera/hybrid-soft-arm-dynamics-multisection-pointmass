# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\candidate07_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `24753e2f2c5590bc24759db599827d77705a55b46536cfbb073f222978a2812f`

Classification: **valid_feasible**. Success: 1. Schema: pac-tune-2.

Campaign classification: **valid_feasible**. Eligible for candidate selection: **1**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Objective: velocityErrorRateRmse_mPerS2 / 0.025 = 0.9050694.

Motion metric source: fixedRawLength.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.98625 | 0.985
thetaCoverage | 0.98625 | 0.985
phiCoverage | 0.98625 | 0.985
velocityCoverage | 0.98625 | 0.985
velocityPairCoverage | 0.9720833 | 0.97
q2Rmse_mm | 2.635934 | 1.236682
q3Rmse_mm | 3.191741 | 1.66263
qVectorRmse_mm | 4.139488 | 2.07213
thetaMae_rad | 0.2605379 | 0.1214906
phiMae_rad | 0.09434197 | 0.06289083
phiMax_rad | 1.253243 | 1.12301
velocityRmse_mPerS | 0.003678133 | 0.003012494
velocityErrorRateRmse_mPerS2 | 0.01935712 | 0.02262673
thetaRateStd_radPerS | 0.3438207 | 0.3826328
nearStopFraction | 0.03814831 | 0.04810997
slowOrReverseFraction | 0.04157737 | 0.05841924
activeFraction | 0.9704167 | 0.965
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 2 | 1.236682 | within limit
q3Rmse_mm | 1 | 2 | 1.66263 | within limit
velocityRmse_mPerS | 1 | 0.01 | 0.003012494 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.02262673 | within limit
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 2.07213 | disabled
thetaMae_rad | 0 | NaN | 0.1214906 | disabled
phiMae_rad | 0 | NaN | 0.06289083 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.3826328 | disabled

Tuning changes: 0. Replay: 2367 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.98625.

Observed parameter rate-cap fractions: `[0 0 0 0 0 0 0 0 0 0]` (active transitions: 2329).

Minimum estimated K eigenvalue: -206.9406 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.14376 1.23326 1.272]`.

Tracking parameter change: `[-258.566 -724.205 -298.92 8.9137 10.8452 10.5808 10.0009 10.6856 -0.738041 -0.00856582]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
