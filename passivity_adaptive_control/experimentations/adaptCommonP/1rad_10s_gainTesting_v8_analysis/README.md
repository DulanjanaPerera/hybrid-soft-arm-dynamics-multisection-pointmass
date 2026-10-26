# v8 saturation and speed-headroom audit

30x10s circles,phi1. Lambda18/16.5,Ks15/10fixed. Gamma k.7->.5at221.65s,
d.002->.02at271.45s,h5.8fixed;caps650/10/2/.2 bygroup;hupper20.

Circleboundcontacts: k2=-500on24samples;k3=1on93;historyweightszeroearly
on16/5/17samples. No upperboundcontacts,no damping/biascontacts. Onlyk3
ratecapreached(25updatesearly). NONEof10parametersatvalueor ratecapsinthe
final3circles290-320s. Finalreturnholdlowerboundcontactsareseparatefromtracking.

Adaptactive96.15%;off62invalid,62dt,107deadband;no pressure/slewlimiting.
Lap8thetaMAE.1983,lap20.0975,lap30.0994rad;phi settles~.062rad.
Final-lengthRMSE.911/1.359mm. Remaininglateerrordoesnotappearbound/rate-limited.
Latevsv7muchlessrough,butnotassmallpositionerror. No latchedfaults;estimated
stiffnessremainsindefiniteon7.25%circlesamples,currentmodelValiddoesnotcheckit.

2sgoal:40PACsamples/cycle;5xvelocityand25xaccelerationrelative10s.
For13mmamplitude,referencepeak40.84mm/s,128.30mm/s^2. Compressingrecorded
pressurewaveform5x(notpredictingactual2sclosedloop)requirespeak7.62/9.57/6.80
bar/svs2bar/slimit;~35-43%lateincrementswouldexceedthelimit. Largeradaptation
gainsalonecannotestablish2stracking;actuatorbandwidth/lagandsignedKneedattention.

audit.txt:actualsettings,gates,perlapmetrics;bound_and_rate_audit.txt:exact
bound/ratecountsandpressureheadroom;v7_v8_roughness.txtandcomparisonfigure;
parameter/gain/switchfigures. Reproducefromroot:tmp/auditPacGainTestingV8.m,
tmp/auditPacV8Limits.m,tmp/comparePacV7V8.m. No controller/modelchanges or
model/hardwareexecution. No commit/push.
