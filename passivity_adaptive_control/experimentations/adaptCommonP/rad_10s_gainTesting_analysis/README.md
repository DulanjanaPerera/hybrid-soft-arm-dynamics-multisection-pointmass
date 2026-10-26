# Gain-testing recording audit — 2026-10-26

Source: `../rad_10s_gainTesting.mat`; 0-330 s, thirty 10 s circles, phi=1 rad.
Payload condition is not established by this filename. All ten tuning vectors
are logged (history options use the recorded spelling `hystoryOptions`).

`audit.txt` contains actual settings, all feedback-gain changes, disjoint update
gate counts, rate-cap fractions, stiffness eigenvalue checks, per-lap tracking,
and metrics per setting interval. The gate reconstruction matches the logged
adaptationActive flag at every sample. Categories use the first blocking gate;
for example, an excessive nominal residual during limiting is counted as limiting.

Circle samples: 6000; updates active5534 (92.233%); inactive466:
358 allocation/slew limiting,67 invalid data/model,41 unacceptable update
interval; no independent residual-only, deadband, manual-disable or unexplained
cases. Invalid gaps can also block the first returning update because the elapsed
measurement interval exceeds .075 s. This is not evidence of convergence.

At 243.75-250 s, second Ks component100: 123 of125 samples are limiting;
adaptation active0.8%,thetaMAE.314676,phiMAE.152063. Adjacent intervals have
lower errors. Gain intervals are partial cycles with evolving adaptive states,
so they do not establish an optimum or controlled comparison. Final lap:
thetaMAE.187799,phiMAE.060273,length RMSE1.967/1.904 mm.

Actual history-weight gamma=1 and rate caps=2 N/s, not .7/.5 from earlier chat.
Spring caps450 reached only .40/.14/1.62% of all active updates; no other
parameter caps reached. Raising caps is not the main opportunity demonstrated.

Signed stiffness becomes indefinite during valid circle samples23.95-25.35 s
(29samples,minimum eigenvalue -87.667N/m) and during final return/hold
(whole-run minimum -170.793N/m). No faults were latched; modelValid does not
test estimated stiffness SPD. Observed bounded tracking does not reinstate the
previous SPD/passivity argument. k3 reaches its1N/m lower bound; k2 never reaches
-500 (minimum -318.435). k2 changes periodically rather than settling at a single
constant; final-after-return k2=-171.791 is not its circle mean.

Figures: adaptation_gate_diagnostics; all ten retained parameter estimates;
feedback gains; adaptation switches. Additional gain/parameter/switch exports
include vector PDF/SVG and editable FIG. lap_metrics.csv columns:
lap,thetaMAE_rad,phiMAE_rad,l2RMSE_mm,l3RMSE_mm,activePercent,deadbandPercent.

Reproduce offline: `run('tmp/auditPacGainTesting.m')` from the project root.
No controller changes, model execution, hardware use, commit or push.
