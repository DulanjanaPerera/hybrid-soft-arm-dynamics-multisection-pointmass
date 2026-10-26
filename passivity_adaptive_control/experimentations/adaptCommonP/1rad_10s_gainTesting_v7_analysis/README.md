# v7 dynamic adaptation-gain audit

Feedbacklambda[18;16.5],Ks[15;10]fixed;gamma k/h/d changes12times.
Springratecaps650instead of450. Historyupperbounds20. See audit.txt for exact
timeline,allsettings,perlaptracking,gates and stiffness eigenvalue checks.

Historygamma5.8from81scoincideswithfastererrorreduction,butk/dgammasalsochange
andallstatespersist:theseare not independent trials. Final-lapthetaMAE.063644,
phi.051364,lengthRMSE.787/1.013mm. Thislowerrorcomeswithworselate motion.

At254.95s,gamma_d.02->2;d estimates49.95/52.40->106.87/117.88by272.6s,
then~112/124late. Loweringgammaagaindoesnotundoexistingd estimates.
Final3circlesnearstopfraction27.4%,phi max1.322rad,pressurelimiting20.2%.
The large learned damping channel is model compensation Dh*v,not equivalent
to raising explicit feedback damping Ks. Nooptimalgamma provedbythissweep.

Proposednexttest: fresh initialization withfixedfeedback18/16.5,15/10 and
gamma_k.5,gamma_h5.8,gamma_d.002,biasgamma1;other v7inputsunchanged.
Thisisacandidatecomparison,notvalidatedoptimum. No automaticmodel/hardwarerun.
SignedKisindefiniteon227circlesamples;noSPD/passivityclaimrestored.

Plots:gamma_and_damping_history.png/pdf/fig;v6_v7_last3cycles.png/fig;
parameter/gain/switchfigures. gamma_stage_audit.txt giveswithin-runroughness.
Reproducefromroot:run tmp/auditPacGainTestingV7.m,tmp/auditPacV7Stages.m,
andtmp/comparePacV6V7.m. No controller/model changes,commit orpush.
