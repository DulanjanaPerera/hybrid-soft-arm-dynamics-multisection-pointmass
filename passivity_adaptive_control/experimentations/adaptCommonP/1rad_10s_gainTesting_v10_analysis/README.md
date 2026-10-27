# v10: four-second circle onset and fault audit

Reference4s/circle,phi1rad,from20s. Recordingends109.6s. Lambda18/16.5,
Ks15/10fixed;slew6bar/s. Initialgamma_k.5,d.02,h5;d->.002at78.7,k->.1at93.5.
No gain change at55s. Parametercaps850/20/5/.2 bygroup,lowerk2-500,hupper20.

The log shows growing within-cycle roughness whilemeantrackingerrorcontinues
falling. At55s,status3,modelvalid,no commandlimit or sensorhold. Estimatedd
hasgrown40/40->86/90N.s/m;~104/109by78.7. Lowergammaafterwarddoesnotresetd.
Kisindefiniteintermittentlyearlier,andcontinuouslyfrom55.2. Thisdoesnotalone
proveinstabilityoridentifyphysicalcause. Unmodeledactuator/filterdynamicsand
evolvingmodelcompensationareplausible;oneclosedlooplogcannotisolatethem.
k2firsthits-500at62s,afterreportedonset;nootherparameterbound/rateeventat55.

Nearstopfractiongrows6.6%(36-44s)->16.8%(52-60)->24.3%(60-68);phi peaks
~1.28then1.40rad. Laterpressurelimitingincreases. Initial20-44sallocatorreplay
showsnolimits,unlikev9'spersistent2bar/sslewblocking.

Separateending:tipinvalid104.75/.80/.85;hold4thenlatchedfault5at104.85after
the.10sgrace. Sensorrecovers104.90,butfaultretains. Themainauditlap22
(104-108s)includesfault/unloading;do notuseitaspuretrackingerror. Use
onset_audit.txtnormalstatus3-onlywindowsforperformance. modelValidalone
doesnotexcludeafaultedpressurecommandpath.

onset_timeline.png/fig alignsphi,velocityerror,d estimates,K eigenvalue and
pressure;verticalmarkers55s,78.7s,andfault104.85. Exactsettings,gatesand
fullperlapvaluesareinaudit.txt. Suggestedfreshgamma_d.002comparisonisdiagnostic,
notvalidatedcure. No controller/model/hardwarechanges. Reproducefromroot:
tmp/auditPacGainTestingV10.m andtmp/auditPacV10Onset.m.
