# Fixed-gain v2 audit

Source: `../1rad_10s_gainTesting_v2.mat`. 30 circles at phi=1 rad and10s/circle.
Lambda [28;25.5] and Ks [25;.1] are constant throughout; all logged adaptation
settings are also constant. Payload condition not established by filename.

Theta MAE decreases from .435323rad onlap2 to .168024 onlap30;phi from .138773
to .054390. Length RMSE improves3.989/4.191 to1.643/1.881mm. See lap_metrics.csv
for each circle (columns lap,thetaMAE,phiMAE,l2RMSEmm,l3RMSEmm,active%,deadband%).

Adaptation is active96.183% of circle samples. First-blocking gate counts:
105limiting,69invalid,55updateinterval,0deadband. Gate reconstruction exactly
matches the logged flag. History/damping rate caps are never reached; spring
caps reached only .687/.285/2.044% of active updates. Gamma could change update
magnitudes below caps, but faster/better physical tracking is not guaranteed.

History weights continue growing; their mean values onlap30 are7.275/6.778/7.307N,
below10N bounds. The update requests shrink as the error correlations change.
This is ongoing compensation learning, not established physical identification.

Estimated stiffness is indefinite early in the run:81circle samples,79modelvalid,
23.4-25.9s and33.9-35.35s;minimum eigenvalue -257.412N/m. No latched faults;
modelValid does not test estimatedK SPD. k2 lowerbound -500 is reached early.

audit.txt contains full settings, gate/rate/eigenvalue checks and per-lap metrics.
PNG/FIG diagnostics and parameter/gain/switch plots are included. No controller
or model changes/hardware runs. Reproduce: `run('tmp/auditPacGainTestingV2.m')`.
