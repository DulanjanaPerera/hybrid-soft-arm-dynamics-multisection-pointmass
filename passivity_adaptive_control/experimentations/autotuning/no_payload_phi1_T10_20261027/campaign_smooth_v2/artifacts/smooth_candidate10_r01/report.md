# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\smooth_candidate10_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `4dd97af7e46fcb588f4af0fb4823efb079670cfccb7b1f9698fccaa5cde86b13`

Classification: **valid_infeasible**. Success: 0. Schema: pac-tune-2.

Campaign classification: **valid_infeasible**. Eligible for candidate selection: **0**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Objective: velocityErrorRateRmse_mPerS2 / 0.025 = 0.6431046.

Motion metric source: fixedRawLength.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.9908333 | 0.995
thetaCoverage | 0.9908333 | 0.995
phiCoverage | 0.9908333 | 0.995
velocityCoverage | 0.9908333 | 0.995
velocityPairCoverage | 0.98125 | 0.99
q2Rmse_mm | 3.9436 | 2.880204
q3Rmse_mm | 4.486152 | 3.166379
qVectorRmse_mm | 5.973068 | 4.280366
thetaMae_rad | 0.4434736 | 0.3106698
phiMae_rad | 0.1425524 | 0.07809105
phiMax_rad | 1.186989 | 1.050749
velocityRmse_mPerS | 0.004427766 | 0.003572203
velocityErrorRateRmse_mPerS2 | 0.01448023 | 0.01607762
thetaRateStd_radPerS | 0.275688 | 0.2778131
nearStopFraction | 0.0118896 | 0.01346801
slowOrReverseFraction | 0.0118896 | 0.01346801
activeFraction | 0.9816667 | 0.99
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 2 | 2.880204 | exceeded/unavailable
q3Rmse_mm | 1 | 2 | 3.166379 | exceeded/unavailable
velocityRmse_mPerS | 1 | 0.01 | 0.003572203 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.01607762 | within limit
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 4.280366 | disabled
thetaMae_rad | 0 | NaN | 0.3106698 | disabled
phiMae_rad | 0 | NaN | 0.07809105 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.2778131 | disabled

Tuning changes: 0. Replay: 2378 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.99083.

Observed parameter rate-cap fractions: `[0 0 0 0 0 0 0 0 0 0]` (active transitions: 2356).

Minimum estimated K eigenvalue: 111.4992 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[0.90536 0.970629 0.866588]`.

Tracking parameter change: `[-2.78266 -710.472 -213.583 4.44266 4.79063 7.34437 6.73563 7.08599 -0.602556 -0.213305]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
