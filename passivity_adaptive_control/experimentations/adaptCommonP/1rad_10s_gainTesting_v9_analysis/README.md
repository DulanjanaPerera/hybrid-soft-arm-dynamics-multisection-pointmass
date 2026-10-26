# v9: actual two-second circle audit

Filename contains10s,but loggedreference is30two-secondcirclesatphi1rad,
20-80s;recordingends89.55s. Lambda18/16.5,Ks15/10constant;gamma_k.5->.05
at63.4s,gamma_h5.8->.1at69.7s,gamma_d.02unchanged.

Mainresult:2bar/s pressurecommandslewlimitblocksnominalrequestedforce.
1071of1200circlesampleslimited(89.25%);adaptationactive65samples(5.42%).
Exactallocatorreplaymatchesalllimits:zero nominalforceinfeasibleallocations,
1071commandsnotyetatallocator targets. Peakcommandresidual9.385N;allocator
targetresidualroundoff. Commandsstaybelow3barupperrails. Thisdoesnotmeasure
physicalchamberpressureorforce. See allocation_and_bounds.txt andfigure.

Final-lapthetaMAE1.067rad,phi.275rad,lengthRMSE8.07/9.76mm. Somecomplete
lapscontainnoadaptiveupdates. k3/h1/h3frequentlyatlowerbounds,butmoregamma
cannotbypassmissingactuationauthority/updategates. Historyratecapsnotreached.
EstimatedKSPDthroughcircles;nofaults. 62invalidcirclesamplesareasecondaryissue.

Do not disablegatingorraiseratecapsasaparameter-learningfix. Physicalactuator
responseandallowablecommandslewneedcharacterizationbeforechangingthe2bar/s
limit;otherwise reducecircleamplitude/speed. No controllerchanges/hardware/model
executionperformed. Reproducefromroot:tmp/auditPacGainTestingV9.m and
tmp/auditPacV9Allocation.m. audit.txt uses2sperlapwindows,not10s.
