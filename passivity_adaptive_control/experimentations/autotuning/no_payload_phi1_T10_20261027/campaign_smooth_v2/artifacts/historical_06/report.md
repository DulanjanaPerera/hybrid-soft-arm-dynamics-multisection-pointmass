# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\candidate06_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `ed66339a1f8574af8fec5c90806c092ebd1bdadde98e7731e75e6add00ed2044`

Classification: **valid_infeasible**. Success: 0. Schema: pac-tune-2.

Campaign classification: **valid_infeasible**. Eligible for candidate selection: **0**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Objective: velocityErrorRateRmse_mPerS2 / 0.025 = 0.7750229.

Motion metric source: fixedRawLength.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.9895833 | 0.9883333
thetaCoverage | 0.9895833 | 0.9883333
phiCoverage | 0.9895833 | 0.9883333
velocityCoverage | 0.9895833 | 0.9883333
velocityPairCoverage | 0.97875 | 0.9766667
q2Rmse_mm | 3.695611 | 2.633751
q3Rmse_mm | 4.208709 | 2.903274
qVectorRmse_mm | 5.600962 | 3.919903
thetaMae_rad | 0.4096935 | 0.2843771
phiMae_rad | 0.1332664 | 0.07111649
phiMax_rad | 1.267785 | 1.045224
velocityRmse_mPerS | 0.004266461 | 0.003455605
velocityErrorRateRmse_mPerS2 | 0.01782514 | 0.01937557
thetaRateStd_radPerS | 0.3192082 | 0.3199823
nearStopFraction | 0.0391656 | 0.03754266
slowOrReverseFraction | 0.04044274 | 0.03924915
activeFraction | 0.9791667 | 0.9766667
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 2 | 2.633751 | exceeded/unavailable
q3Rmse_mm | 1 | 2 | 2.903274 | exceeded/unavailable
velocityRmse_mPerS | 1 | 0.01 | 0.003455605 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.01937557 | within limit
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 3.919903 | disabled
thetaMae_rad | 0 | NaN | 0.2843771 | disabled
phiMae_rad | 0 | NaN | 0.07111649 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.3199823 | disabled

Tuning changes: 0. Replay: 2375 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.98958.

Observed parameter rate-cap fractions: `[0 0 0 0 0 0 0 0 0 0]` (active transitions: 2350).

Minimum estimated K eigenvalue: -321.0461 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.12304 1.51758 1.49994]`.

Tracking parameter change: `[5.1413 -793.873 -124.265 11.3876 12.3317 7.00572 6.34075 6.61367 -0.556745 -0.338725]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
