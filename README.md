# Single-module feedback-control baseline

This branch contains the starting pieces for a new feedback-linearization model. It does not contain a controller or the prior dynamics-comparison experiments.

## Models

- NDI_readingOnly.slx reads two-sensor XYZ and orientation, with the confirmed arm/NDI axis rotation and offsets. It uses COM12, 921600 baud, one BX poll per paired sample, sensor 0A as base, sensor 0B as tip, and a nominal 0.05 s sample period. Its MATLAB dependencies are AuroraReader_s2_Position_only.m, AuroraDriver_2.m, PortHandle.m, ndiSupportFnc.m, and f20260219_2_task2config_withLenExt.m. Geometry, acquisition defaults, sample validation, and XYZ comparison are methods in ndiSupportFnc.m.
- Pressure_actuationOnly.slx is the original NI/SMC pressure-regulator model. Its NI Analog Output blocks address physical hardware when the model runs.

## Dynamics

armS_single_entry.m is the single-section RHS with q = [dl2; dl3]. It calls armS_single_dynamics.m and armS_single_core.m, which use the three included math helpers. The RHS accepts a full 2-by-2 stiffness matrix. A nominal coupled 1350 N/m axial stiffness corresponds to K = 1350*[2 1; 1 2]. The confirmed module length is 0.17411 m, actuator radius is 0.013 m, and the arm hangs downward; supply gravity in the arm frame when using the RHS.

The NDI model reports calibrated, base-relative sensor-orientation angles at its main `theta` and `phi` ports. It retains XYZ inverse-kinematics estimates at `thetaXYZ` and `phiXYZ` for comparison, and exposes `lengthChange_m`, `orientationValid`, `thetaObservable`, and `nonBendingRotation_rad`. The reader obtains XYZ and quaternions from the same BX poll; `ndiSupportFnc.m` contains orientation, calibration, preflight, and angle-to-length methods.

On a missing paired reading, the Scope traces hold their most recent valid position and angles instead of leaving a gap. The raw XYZ path into inverse kinematics remains NaN on an invalid sample, and `sensorValid`, `configurationValid`, and `orientationValid` remain false; held outputs are stale display values, not new measurements. Before the first valid reading (or before theta becomes observable), a trace can still be NaN.

With the pressure model stopped and the arm held straight, untwisted, and still, run `[cal, calibrationFile] = ndiSupportFnc.runCalibration` in MATLAB. It records about five seconds at the configured 0.05 s interval and saves the raw poses and fitted `cal` together in one MAT-file under `ndi_orientation_recordings`. Fixed lateral sensor offsets are recorded as the reference rather than forced to zero. A genuinely bent reference pose would be mistaken for a mounting rotation, so confirm the pose physically and validate the resulting angles with independent bends. The September 28 calibration files should be treated as stale if either sensor mount moved.

Set the MATLAB System block's `CalibrationFile` property to the full path of that verified MAT-file. Each reader startup then takes a five-second preflight while the arm is held at the same straight reference pose. It checks paired readings, actual poll timing, pose stability, relative orientation, and displacement against the saved reference. Failure stops setup and closes COM12; it does not silently recalibrate. The default is controlled by the block's `PreflightSeconds` property. For a separate preflight with the reader stopped, run `[report, raw] = ndiSupportFnc.preflight(calibrationFile)`. Only one process may own COM12 at a time. While `CalibrationFile` is blank, orientation angles and lengths remain NaN and `orientationValid` is false; the XYZ path still operates without preflight. Keep the pressure model stopped while checking the NDI measurement chain.

The preflight's initial engineering limits are at least 90% valid pairs, median poll interval at most 1.5 times the configured period, median bend and extra rotation at most 5 degrees, relative XYZ change at most 15 mm from the saved pose, and 95th-percentile motion within 2 degrees and 3 mm. These are editable gates in `ndiSupportFnc.evaluatePreflight`, not NDI accuracy specifications. Review the returned `report` and raw sample timing if a gate fails.

The Aurora driver starts by trying `INIT` at its host default of 9600 baud. If the tracker does not reply because a previous session left it at the configured 921600 baud, the driver switches only the host serial settings and retries `INIT`. It also retries one `ERROR02` response after clearing stale input, since a wrong-baud probe may leave a delayed malformed-command reply. Other errors stop startup. A successful 9600-baud startup still changes the tracker to 921600 for acquisition. The legacy MATLAB `serial` deprecation warning does not by itself indicate an acquisition failure.

No experiment recordings, MEX binaries, old three-section models, or prior feedback files are tracked on this branch. MATLAB and the NI/SMC hardware support installed on the experiment computer are still required for the physical Simulink models.

## Offline single-module animation

Run `runDynamicSimulation_armS_single_custom` in MATLAB to solve and animate the single-module dynamics without connecting to NDI or the pressure valves. Edit the parameter block at the top of the script for initial length changes, signed pressure input, base orientation, and recording options. The defaults use the confirmed 0.17411 m module length, a downward-pointing arm, and the nominal coupled stiffness `1350*[2 1; 1 2]` N/m. The script leaves `t`, `X`, `X0`, and `params` in the workspace. Set `recordVideo=true` to save an MP4 and MAT file under `single_module_recordings`.

The animation uses `drawingArm_single_stationary_base.m` and `HTM_nume.m`; the simulation uses the existing `armS_single_entry.m` dynamics. No compiled MEX is needed.
