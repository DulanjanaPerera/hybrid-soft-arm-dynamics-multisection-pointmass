# v4 vs v3: second feedback damping component 0.1 to 1

Source: `../1rad_10s_gainTesting_v4.mat`;30x10s circles,phi1. Loggedlambda
[28;25.5],Ks[25;1] constant;other logged tuning matchesv3(historygamma1.5).
Payload/preconditioning equivalence and repeatability not established.

Final3 mean-per-lap results (v3 -> v4):thetaMAE .121886 -> .115522rad,
phiMAE .052813 -> .052724rad,lengthRMSE1.241/1.541 ->1.211/1.453mm.
This is a small improvement, not a demonstrated cure for within-cycle roughness.

Matching290-320s valid-data velocity/roughness metrics:
l3velocity-errorRMS2.699->2.564mm/s;peak15.917->15.732mm/s.
Apparenttheta-rateSD.514->.505rad/s;nearstopfraction14.555->14.576%;
max2.591->3.285rad/s. Rate is sampled angle difference,not independent speed.
Phi maximum1.220->1.196rad,buttimeabove1.05rad9.459->11.765%.
Limiting3->3.333%. Thus improvements are mixed;stop-burst motion remains.

Circleadaptationactive96%;firstblockedgates137limited,54invalid,42dt,7deadband.
No latchedfaults. Historyweightsapproach9.264Nof10N;ratecapsnotlimitinghistory.
EstimatedK remainsindefiniteearlycycles/return;noSPD/passivity claim restored.

audit.txt:actualsettings,gates,perlapmetrics;lap_metrics.csv;
v3_v4_roughness.txt;v3_v4_last3cycles.png/fig;parameter/gain/switchplots.
Reproducefromroot:`run('tmp/auditPacGainTestingV4.m')`,
then`run('tmp/comparePacV3V4.m')`. No controller/model/hardware execution.
