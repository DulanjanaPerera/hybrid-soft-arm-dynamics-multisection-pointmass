# v3 vs v2: history adaptation gain comparison

Source: `../1rad_10s_gainTesting_v3.mat`, compared with v2. Thirty10s circles,
phi1rad,lambda[28;25.5],Ks[25;.1]. Among logged tuning inputs, only history
weight adaptation gamma changes1->1.5. Payload/preconditioning not established.

Final lap thetaMAE .118949 vs .168024rad;phi .052894 vs .054390rad;lengthRMSE
1.200/1.543 vs1.643/1.881mm. Theta firstbelow.2rad onlap17 vs24. Learning is
faster, but continued adaptation/physical parameter identification not proven.

Late290-320s l3 velocity peaks15.917 vs16.067mm/s: bursts persist. Apparent
theta-rate SD.514 vs.520rad/s,slow-motion fraction14.56% vs14.73%;largest
apparenttheta-rate peakdecreases2.59 vs3.34rad/s. These are sampled angle
differences, not independent velocity measurements. Phi maxincreases1.220 vs
1.168rad;timeabove1.05rad9.46% vs6.25%. Fasterlearning doesnotresolve roughness.

Adaptationactive96.02%;circleoff118limiting,64invalid,54dt,3deadband. History
rate capsnotreached;historyweightsapproach9.25N of10Nbounds. No latchedfaults.
EstimatedK indefiniteearlycyclesandreturn/hold;noSPD/passivity claim restored.

Files: audit.txt(actualsettings,gatecounts,perlapmetrics),lap_metrics.csv,
v2_v3_roughness.txt(comparisondefinitions/results),v2_v3_last3cycles.png/fig,
gate diagnostics and parameter/gain/switch plots. No controller/model changes.
Reproduce from root: `run('tmp/auditPacGainTestingV3.m')`, then
`run('tmp/comparePacV2V3.m')`. No hardware/model execution.
