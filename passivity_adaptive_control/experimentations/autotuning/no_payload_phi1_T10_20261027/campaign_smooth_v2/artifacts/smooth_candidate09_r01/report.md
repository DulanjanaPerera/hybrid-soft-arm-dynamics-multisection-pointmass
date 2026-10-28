# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\smooth_candidate09_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `0e48a30a4fab7c95ffc1e2c083871532dd0b26d0a9d2fa242dd5ebe82ca760a9`

Classification: **valid_infeasible**. Success: 0. Schema: pac-tune-2.

Campaign classification: **valid_infeasible**. Eligible for candidate selection: **0**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Objective: velocityErrorRateRmse_mPerS2 / 0.025 = 0.8830536.

Motion metric source: fixedRawLength.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.99 | 0.99
thetaCoverage | 0.99 | 0.99
phiCoverage | 0.99 | 0.99
velocityCoverage | 0.99 | 0.99
velocityPairCoverage | 0.9795833 | 0.98
q2Rmse_mm | 3.731408 | 2.527911
q3Rmse_mm | 4.043302 | 2.752004
qVectorRmse_mm | 5.501972 | 3.736825
thetaMae_rad | 0.4091306 | 0.2708057
phiMae_rad | 0.1249676 | 0.07273854
phiMax_rad | 1.333286 | 1.073127
velocityRmse_mPerS | 0.004397457 | 0.003515224
velocityErrorRateRmse_mPerS2 | 0.02071302 | 0.02207634
thetaRateStd_radPerS | 0.3417199 | 0.3587269
nearStopFraction | 0.05869843 | 0.06802721
slowOrReverseFraction | 0.06465334 | 0.07653061
activeFraction | 0.98 | 0.98
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 2 | 2.527911 | exceeded/unavailable
q3Rmse_mm | 1 | 2 | 2.752004 | exceeded/unavailable
velocityRmse_mPerS | 1 | 0.01 | 0.003515224 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.02207634 | within limit
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 3.736825 | disabled
thetaMae_rad | 0 | NaN | 0.2708057 | disabled
phiMae_rad | 0 | NaN | 0.07273854 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.3587269 | disabled

Tuning changes: 0. Replay: 2376 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.99000.

Observed parameter rate-cap fractions: `[0 0 0.0097789 0 0 0 0 0 0 0]` (active transitions: 2352).

Minimum estimated K eigenvalue: -507.1893 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.12632 1.5685 1.719]`.

Tracking parameter change: `[-116.872 -938.12 -82.3988 12.8309 12.6245 7.87985 6.91132 6.66512 -0.819711 -0.321337]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
