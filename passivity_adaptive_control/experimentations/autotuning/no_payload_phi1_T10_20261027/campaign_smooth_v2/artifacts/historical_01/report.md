# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\baseline_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `fc508a484416369b2e4786bd3a5581b27779a9a3a8b15747d3c2507370c69794`

Classification: **valid_feasible**. Success: 1. Schema: pac-tune-2.

Campaign classification: **valid_feasible**. Eligible for candidate selection: **1**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Objective: velocityErrorRateRmse_mPerS2 / 0.025 = 0.7748384.

Motion metric source: fixedRawLength.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.9916667 | 0.9866667
thetaCoverage | 0.9916667 | 0.9866667
phiCoverage | 0.9916667 | 0.9866667
velocityCoverage | 0.9916667 | 0.9866667
velocityPairCoverage | 0.9829167 | 0.9733333
q2Rmse_mm | 2.943851 | 1.568437
q3Rmse_mm | 3.415848 | 1.888731
qVectorRmse_mm | 4.509355 | 2.455056
thetaMae_rad | 0.3010761 | 0.1520176
phiMae_rad | 0.1049582 | 0.06247963
phiMax_rad | 1.2718 | 1.111328
velocityRmse_mPerS | 0.003805289 | 0.003014625
velocityErrorRateRmse_mPerS2 | 0.0184934 | 0.01937096
thetaRateStd_radPerS | 0.324156 | 0.3344155
nearStopFraction | 0.04239084 | 0.05308219
slowOrReverseFraction | 0.04408648 | 0.05308219
activeFraction | 0.9833333 | 0.9733333
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 2 | 1.568437 | within limit
q3Rmse_mm | 1 | 2 | 1.888731 | within limit
velocityRmse_mPerS | 1 | 0.01 | 0.003014625 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.01937096 | within limit
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 2.455056 | disabled
thetaMae_rad | 0 | NaN | 0.1520176 | disabled
phiMae_rad | 0 | NaN | 0.06247963 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.3344155 | disabled

Tuning changes: 0. Replay: 2380 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.99167.

Observed parameter rate-cap fractions: `[0 0 0.0038136 0 0 0 0 0 0 0]` (active transitions: 2360).

Minimum estimated K eigenvalue: -344.698 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.11271 1.10716 1.09687]`.

Tracking parameter change: `[-239.976 -719.402 -266.077 7.27136 8.16864 10.1056 9.24897 9.64365 -0.880395 -0.243344]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
