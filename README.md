# Single-module feedback-control baseline

This branch contains the starting pieces for a new feedback-linearization model. It does not contain a controller or the prior dynamics-comparison experiments.

## Models

- NDI_readingOnly.slx is the original XYZ-based two-sensor reader, with the confirmed arm/NDI axis rotation and offsets. It uses COM12, 921600 baud, one BX poll per paired sample, sensor 0A as base, sensor 0B as tip, and a nominal 0.05 s sample period. Its small MATLAB dependencies are AuroraReader_s2_Position_only.m, AuroraDriver_2.m, PortHandle.m, ndiAcquisitionSettings.m, ndiValidSample.m, ndiConfiguration.m, ndiSensorGeometry.m, and f20260219_2_task2config_withLenExt.m.
- Pressure_actuationOnly.slx is the original NI/SMC pressure-regulator model. Its NI Analog Output blocks address physical hardware when the model runs.

## Dynamics

armS_single_entry.m is the single-section RHS with q = [dl2; dl3]. It calls armS_single_dynamics.m and armS_single_core.m, which use the three included math helpers. The RHS accepts a full 2-by-2 stiffness matrix. A nominal coupled 1350 N/m axial stiffness corresponds to K = 1350*[2 1; 1 2]. The confirmed module length is 0.17411 m, actuator radius is 0.013 m, and the arm hangs downward; supply gravity in the arm frame when using the RHS.

The baseline NDI model estimates angles from XYZ. It does not use the later orientation calibration or feedback-controller classes. Secure the sensors and verify straight-pose geometry before using their measurements for control. Keep the pressure model stopped while checking the NDI measurement chain.

No experiment recordings, MEX binaries, old three-section models, or prior feedback files are included. MATLAB and the NI/SMC hardware support installed on the experiment computer are still required for the physical Simulink models.
