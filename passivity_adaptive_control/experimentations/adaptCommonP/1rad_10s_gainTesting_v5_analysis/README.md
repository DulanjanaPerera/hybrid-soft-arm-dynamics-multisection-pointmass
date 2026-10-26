# v5: smoother feedback,with a lag tradeoff

Source `../1rad_10s_gainTesting_v5.mat`.30x10s circles,phi1;payloadunverified.
Gamma is unchanged fromv4: historyweights1.5. Lambda changes28/25.5->28/16.5
at40s->18/16.5at59.5s;Ks25/1->25/10at70.25s->15/10at88.2s.
After88.2s feedbacksettingsremainfixed. Otherloggedtuningmatchesv4.

Last3meanlaperrors(v4->v5):theta.115522->.195919rad,phi.052724->.064290rad,
length1.211/1.453->1.927/2.263mm. Tracking is less accurate,but much smoother.
Afterfinalgainchange,thetaMAE falls.488633lap8->.286515lap20->.188067lap30.

Last3roughness(v4->v5): apparenttheta-rateSD.505->.317rad/s,nearstopfraction
14.576->3.242%,max3.285->1.544rad/s. These are sampledangleincrements,not
independentmechanicalvelocity. l3filteredvelocitypeak15.732->11.786mm/s;
velocityerrorRMS2.564->2.368mm/s. Phi max1.196->1.107rad;no pressurelimiting.

Adaptactive97.933%circle samples;off62invalid+62dt,nonebydeadbandorlimiting.
Historyweightsremainbelow10N,butgrowtoward~8.5-8.9N. EstimatedK indefinite
earlyandduringreturn;noSPDclaimrestored. No latchedfaults.

This supports a smoother feedback operating point and continued adaptation,
but does not isolate gamma:it wasunchanged. Firstrepeatwithfinalfeedbackgains
fixedfromstartupbeforecomparinggammas. No controller/model/hardwarechanges.

audit.txt:settings/gates/perlapmetrics;lap_metrics.csv;v4_v5_roughness.txt;
v4_v5_last3cycles.png/fig;parameter/gain/switchplots. Reproducefromroot:
`run('tmp/auditPacGainTestingV5.m')`,then`run('tmp/comparePacV4V5.m')`.
