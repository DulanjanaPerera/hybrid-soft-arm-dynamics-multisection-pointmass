# v6: spring adaptation gamma 0.8

Recorded springgammas .8/.8/.8,historygamma1.5;lambda18/16.5,Ks15/10fixed
fromstartup. Historyupperbounds20vs10inv5,nonbindinginboth. v5hadearlyfeedback
changes,so v5vsv6doesnotcleanlyisolate springgamma. No sourcetuningcodechanged.

Final-lapthetaMAE .192586vs v5.188067rad,phi.060439vs.064654rad;lengthRMSE
1.867/2.148vs1.884/2.201mm. No clearoverallperformancegain;learningcontinues.
Springratecaps450hit15.8/15.5/10.1%activeupdatesvsunder.6%inv5.

Late290-320s thetaRateSD.320vs.317rad/s,nearstop4.10vs3.24%,peak2.23vs
1.54rad/s;sampledanglevelocitynotindependentmeasurement. l3velocityerrorRMS
2.302vs2.368mm/s;phi peak1.085vs1.107rad. Mixedsmoothnesschanges.

Adaptactive98.17%,no circlelimiting/deadband,55invalid+55dtgates,nofaults.
Estimatedstiffnessindefinite774circlesamples(12.9%,761valid),vs44(.7%)inv5.
Thiscontinueswellbeyondstartup(until~239.65s);currentmodelValiddoesnotcheckKSPD.
CoupledSPDprojection/stabilityanalysisremainsanunimplementeddesignissue.

Nextcomparisonproposed:repeatfixedv6feedbackwithspringgamma.4,historygamma1.5,
andallotherinputsunchanged;donotraisecapsbasedonthisrun. No hardware/model
execution. audit.txt,lap_metrics.csv,v5_v6_roughness.txt,v5_v6_last3cycles.png/fig
andgain/parameter/switchplotsprovideevidence. Reproducefromroot:
`run('tmp/auditPacGainTestingV6.m')`,then`run('tmp/comparePacV5V6.m')`.
