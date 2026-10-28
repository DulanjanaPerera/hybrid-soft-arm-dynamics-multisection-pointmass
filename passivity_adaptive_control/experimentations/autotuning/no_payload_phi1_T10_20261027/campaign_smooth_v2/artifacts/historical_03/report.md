# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\candidate03_00g_phi1_T10_12cycles_r01.mat`

SHA-256: `93f1c6298709cde0f694639628e975b01b75a3a98d42d2797daa30e27ac9ff71`

Classification: **valid_feasible**. Success: 1. Schema: pac-tune-2.

Campaign classification: **valid_feasible**. Eligible for candidate selection: **1**.


The success flag above describes evaluator criteria only; campaign eligibility also checks metadata/configuration.


Objective: velocityErrorRateRmse_mPerS2 / 0.025 = 0.8319052.

Motion metric source: fixedRawLength.


Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 110
end_s | 140 | 140
samples | 2400 | 600
coverage | 0.98875 | 0.9883333
thetaCoverage | 0.98875 | 0.9883333
phiCoverage | 0.98875 | 0.9883333
velocityCoverage | 0.98875 | 0.9883333
velocityPairCoverage | 0.9775 | 0.9766667
q2Rmse_mm | 2.900247 | 1.425487
q3Rmse_mm | 3.267694 | 1.690009
qVectorRmse_mm | 4.369126 | 2.210915
thetaMae_rad | 0.2866861 | 0.1361176
phiMae_rad | 0.1046765 | 0.06543739
phiMax_rad | 1.298536 | 1.131219
velocityRmse_mPerS | 0.003861928 | 0.003163922
velocityErrorRateRmse_mPerS2 | 0.01831814 | 0.02079763
thetaRateStd_radPerS | 0.322484 | 0.3662152
nearStopFraction | 0.04433078 | 0.06313993
slowOrReverseFraction | 0.04475703 | 0.06484642
activeFraction | 0.9775 | 0.9766667
limitedFraction | 0 | 0
faultSamples | 0 | 0

Acceptance limits (scoring window)

Metric | Required | Limit | Observed | Status
--- | --- | ---: | ---: | ---
q2Rmse_mm | 1 | 2 | 1.425487 | within limit
q3Rmse_mm | 1 | 2 | 1.690009 | within limit
velocityRmse_mPerS | 1 | 0.01 | 0.003163922 | within limit
velocityErrorRateRmse_mPerS2 | 1 | 0.025 | 0.02079763 | within limit
limitedFraction | 1 | 0 | 0 | within limit
qVectorRmse_mm | 0 | NaN | 2.210915 | disabled
thetaMae_rad | 0 | NaN | 0.1361176 | disabled
phiMae_rad | 0 | NaN | 0.06543739 | disabled
thetaRateStd_radPerS | 0 | NaN | 0.3662152 | disabled

Tuning changes: 0. Replay: 2373 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.98875.

Observed parameter rate-cap fractions: `[0 0 0.0072464 0 0 0 0 0 0 0]` (active transitions: 2346).

Minimum estimated K eigenvalue: -444.9931 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[1.12069 1.24541 1.31288]`.

Tracking parameter change: `[-319.002 -746.89 -275.963 3.51305 3.73485 11.0513 10.0329 10.2528 -1.08872 -0.259753]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
