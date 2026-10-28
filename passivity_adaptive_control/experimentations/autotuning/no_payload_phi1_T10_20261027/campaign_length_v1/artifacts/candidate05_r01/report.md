# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\candidate05_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `37bd0ec3b39cf656033a2d4aa75230c7c0a37393e35491faf25d0d91f75cb808`

Classification: **valid_infeasible**. Success: 0. Schema: pac-tune-2.

Campaign classification: **valid_infeasible**. Eligible for candidate selection: **0**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.9916667 | 0.9933333
thetaCoverage | 0.9916667 | 0.9933333
phiCoverage | 0.9916667 | 0.9933333
velocityCoverage | 0.9916667 | 0.9933333
velocityPairCoverage | 0.9829167 | 0.9866667
q2Rmse_mm | 3.008784 | 1.836021
q3Rmse_mm | 3.571382 | 2.068622
qVectorRmse_mm | 4.669855 | 2.765894
thetaMae_rad | 0.3233841 | 0.1840056
phiMae_rad | 0.1025942 | 0.06195131
phiMax_rad | 1.236063 | 1.195993
velocityRmse_mPerS | 0.004034063 | 0.003724713
velocityErrorRateRmse_mPerS2 | 0.02528716 | 0.03358536
thetaRateStd_radPerS | 0.4063129 | 0.4858277
nearStopFraction | 0.05595591 | 0.1030405
slowOrReverseFraction | 0.06782535 | 0.1182432
activeFraction | 0.9833333 | 0.9866667
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 1 | 1.836021 | exceeded/unavailable
q3Rmse_mm | 1 | 1 | 2.068622 | exceeded/unavailable
velocityRmse_mPerS | 1 | 0.01 | 0.003724713 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.03358536 | exceeded/unavailable
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 2.765894 | disabled
thetaMae_rad | 0 | NaN | 0.1840056 | disabled
phiMae_rad | 0 | NaN | 0.06195131 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.4858277 | disabled

Tuning changes: 0. Replay: 2380 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.99167.

Observed parameter rate-cap fractions: `[0 0 0 0 0 0 0 0 0 0]` (active transitions: 2360).

Minimum estimated K eigenvalue: -0.4628236 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.00632 2.72223 3.71942]`.

Tracking parameter change: `[-28.7872 -701.82 -188.725 29.5047 34.5795 7.5611 6.96829 7.42652 -0.471032 0.0233324]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
