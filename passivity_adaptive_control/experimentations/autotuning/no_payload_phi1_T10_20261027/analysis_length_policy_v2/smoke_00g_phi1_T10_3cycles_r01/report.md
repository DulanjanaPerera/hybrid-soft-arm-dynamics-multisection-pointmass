# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\autotuning\no_payload_phi1_T10_20261027\smoke_00g_phi1_T10_3cycles_r01.mat`

SHA-256: `f24964411fd146d0ad8ebddd535458e8d60c8952d50ff7eb7e73369a5e369b5e`

Classification: **unconfigured_targets**. Success: 0. Schema: pac-tune-2.

- Targets are unset; historical diagnostics cannot certify success.

Reference period: 10 s; direction: 1. First fault: NaN s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 40
end_s | 50 | 50
samples | 600 | 200
coverage | 0.9866667 | 0.985
thetaCoverage | 0.9866667 | 0.985
phiCoverage | 0.9866667 | 0.985
velocityCoverage | 0.9866667 | 0.985
velocityPairCoverage | 0.9716667 | 0.97
q2Rmse_mm | 4.190861 | 4.252805
q3Rmse_mm | 5.075567 | 4.337554
qVectorRmse_mm | 6.58215 | 6.074597
thetaMae_rad | 0.4858468 | 0.4668396
phiMae_rad | 0.2050988 | 0.1504833
phiMax_rad | 1.265127 | 1.035594
velocityRmse_mPerS | 0.004954952 | 0.004420634
velocityErrorRateRmse_mPerS2 | 0.01746066 | 0.01744866
thetaRateStd_radPerS | 0.3414267 | 0.2938683
nearStopFraction | 0.04459691 | 0.06701031
slowOrReverseFraction | 0.04459691 | 0.06701031
activeFraction | 0.9733333 | 0.97
limitedFraction | 0 | 0
faultSamples | 0 | 0

Tuning changes: 0. Replay: 592 samples, 0 allocation limits, 0 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.98667.

Observed parameter rate-cap fractions: `[0 0 0.013699 0 0 0 0 0 0 0]` (active transitions: 584).

Minimum estimated K eigenvalue: -338.1877 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[0.964327 0.917095 0.810634]`.

Tracking parameter change: `[117.619 -909.602 46.3822 3.52913 4.31943 4.39536 3.6062 3.79406 -0.6219 -0.437485]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
