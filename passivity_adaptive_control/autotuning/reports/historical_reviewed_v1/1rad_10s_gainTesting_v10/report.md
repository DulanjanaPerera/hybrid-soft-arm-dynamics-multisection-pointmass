# PAC offline evaluation

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\passivity_adaptive_control\experimentations\adaptCommonP\1rad_10s_gainTesting_v10.mat`

SHA-256: `7a91a002dcb9f9c666450f593e6723e3f5403e282efe183e1afd54c1af89ffb5`

Classification: **aborted_or_faulted**. Success: 0. Schema: pac-tune-1.

- Targets are unset; historical diagnostics cannot certify success.
- Recording ends before planned tracking completes.
- Coverage below configured minimum of 0.95.
- Live tuning changes detected; exploratory trial.

Reference period: 4 s; direction: 1. First fault: 104.85 s.

Metric | Tracking | Scoring
--- | ---: | ---:
start_s | 20 | 94
end_s | 140 | 102
samples | 1793 | 160
coverage | 0.6829167 | 0.9625
thetaCoverage | 0.6829167 | 0.9625
q2Rmse_mm | 3.590972 | 2.710495
q3Rmse_mm | 4.037489 | 2.946797
qVectorRmse_mm | 5.403369 | 4.003798
thetaMae_rad | 0.3697526 | 0.2442787
phiMae_rad | 0.152838 | 0.1192631
phiMax_rad | 1.401232 | 1.401232
velocityRmse_mPerS | 0.01284579 | 0.01184815
thetaRateStd_radPerS | 1.234153 | 1.226767
nearStopFraction | 0.1403398 | 0.1756757
slowOrReverseFraction | 0.1736941 | 0.2162162
activeFraction | 0.8449526 | 0.8625
limitedFraction | 0.05744562 | 0.0875
faultSamples | 96 | 0

Tuning changes: 2. Replay: 1639 samples, 6 allocation limits, 92 target/command mismatches, 0 logged flag disagreements.

Explicit constant-input extensions: enableAdaptation, enableCommonPressureAdaptation, freezeHysteresisWeights. Aligned finite coverage: 0.91411.

Observed parameter rate-cap fractions: `[0 0.0039604 0.0046205 0 0 0 0 0 0 0]` (active transitions: 1515).

Minimum estimated K eigenvalue: -312.6277 N/m (diagnostic only).

Peak commanded pressure rates (bar/s): `[6 6 6]`.

Tracking parameter change: `[11.7641 -1104.82 -423.073 65.4748 70.6958 6.14877 5.32356 6.00374 -0.509505 -0.0437621]`.

Tracking masks require normal status, validity and freshness. Rates are sampled angle increments. Nominal residual is not measured force error. No hardware was run.
