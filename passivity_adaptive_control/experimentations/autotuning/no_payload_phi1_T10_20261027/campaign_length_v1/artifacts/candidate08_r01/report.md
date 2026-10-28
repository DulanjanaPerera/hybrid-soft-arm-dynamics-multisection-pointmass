# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\candidate08_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `5c4fc2183d529392ee88fa04f07ab48937f9eef6f16a34afda3fd1e7ea98c261`

Classification: **valid_infeasible**. Success: 0. Schema: pac-tune-2.

Campaign classification: **valid_infeasible**. Eligible for candidate selection: **0**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.99 | 0.995
thetaCoverage | 0.99 | 0.995
phiCoverage | 0.99 | 0.995
velocityCoverage | 0.99 | 0.995
velocityPairCoverage | 0.9795833 | 0.99
q2Rmse_mm | 2.545854 | 1.172517
q3Rmse_mm | 2.743959 | 1.302891
qVectorRmse_mm | 3.743085 | 1.752804
thetaMae_rad | 0.2355296 | 0.1097406
phiMae_rad | 0.09625459 | 0.07187077
phiMax_rad | 1.296659 | 1.167303
velocityRmse_mPerS | 0.00366216 | 0.00299378
velocityErrorRateRmse_mPerS2 | 0.02339203 | 0.02391251
thetaRateStd_radPerS | 0.381802 | 0.3957877
nearStopFraction | 0.06975755 | 0.08080808
slowOrReverseFraction | 0.07868992 | 0.09090909
activeFraction | 0.9733333 | 0.9733333
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 1 | 1.172517 | exceeded/unavailable
q3Rmse_mm | 1 | 1 | 1.302891 | exceeded/unavailable
velocityRmse_mPerS | 1 | 0.01 | 0.00299378 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.02391251 | within limit
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 1.752804 | disabled
thetaMae_rad | 0 | NaN | 0.1097406 | disabled
phiMae_rad | 0 | NaN | 0.07187077 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.3957877 | disabled

Tuning changes: 0. Replay: 2376 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.99000.

Observed parameter rate-cap fractions: `[0 0 0.0094178 0 0 0 0 0 0 0]` (active transitions: 2336).

Minimum estimated K eigenvalue: -471.9256 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.02776 1.54975 1.93927]`.

Tracking parameter change: `[-443.876 -650.846 -303.077 21.6378 21.2447 11.1697 9.93564 10.4169 -1.39938 -0.718804]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
