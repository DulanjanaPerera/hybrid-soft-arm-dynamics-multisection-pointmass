# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\filter_smoke_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `f9c82bcbb02ef4f1e6400f694436d298ab56de59cfaa463b50dcadea3939f2f9`

Classification: **valid_infeasible**. Success: 0. Schema: pac-tune-2.


Objective: velocityErrorRateRmse_mPerS2 / 0.025 = 1.363557.

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
q2Rmse_mm | 2.60587 | 1.189358
q3Rmse_mm | 2.899842 | 1.416503
qVectorRmse_mm | 3.898672 | 1.849609
thetaMae_rad | 0.2479151 | 0.1071551
phiMae_rad | 0.08785485 | 0.06118963
phiMax_rad | 1.295139 | 1.211117
velocityRmse_mPerS | 0.003968155 | 0.00357242
velocityErrorRateRmse_mPerS2 | 0.02891259 | 0.03408892
thetaRateStd_radPerS | 0.4500196 | 0.5226175
nearStopFraction | 0.08535032 | 0.1010101
slowOrReverseFraction | 0.104034 | 0.1346801
activeFraction | 0.9770833 | 0.9716667
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 2 | 1.189358 | within limit
q3Rmse_mm | 1 | 2 | 1.416503 | within limit
velocityRmse_mPerS | 1 | 0.01 | 0.00357242 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.03408892 | exceeded/unavailable
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 1.849609 | disabled
thetaMae_rad | 0 | NaN | 0.1071551 | disabled
phiMae_rad | 0 | NaN | 0.06118963 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.5226175 | disabled

Tuning changes: 0. Replay: 2378 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.99083.

Observed parameter rate-cap fractions: `[0 0 0.0081023 0 0 0 0 0 0 0]` (active transitions: 2345).

Minimum estimated K eigenvalue: -452.9043 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.14394 1.94221 2.5648]`.

Tracking parameter change: `[-392.335 -730.932 -362.454 22.868 23.7126 11.8044 10.0108 10.4578 -0.732885 -0.25217]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
