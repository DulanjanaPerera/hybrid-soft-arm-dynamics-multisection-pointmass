# RL data report: synthetic_diagnostics

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\ndi_diagnostics\validation_20261102_v1\synthetic_diagnostics.mat`

## Coverage and quality

- Sensor rows: 3, simulation interval 0.000..0.100 s.
- Valid orientation: 3/3; fresh paired frames: 3/3.
- Observable/valid theta: 3/3.
- Final command range: [0.5 0.5 0.5] to [0.5 0.5 0.5] bar.

## Experiment intervals



Transitions crossing host acquisition gaps > 0.100 s are flagged, not silently used as uniform steps. The threshold is an export-quality setting, not a measured controller period. First/invalid transitions are also rejected.

## Signal interpretation

Commands are pressure proxies, not chamber measurements. des_theta/des_phi are forward-model predictions. theta_reference_rad/phi_reference_rad, when present, are requested circle/helix angles, not observed targets. pressureDesired_bar is the inverse-model calculation and may drive the circle/helix feedforward path; final NI commands remain authoritative actions. Inferred lengths are not independent strain measurements. Corrected tip XYZ uses measured sensor positions plus orientation and assumed offsets. Native CSVs retain invalid/held source values; observations masks invalid theta/phi/inferred lengths. Use explicit per-feature validity; do not assume the whole row is valid for every target.

## Export and timing

observations.csv uses NDI simulation timestamps; commands.csv preserves the union of NI command log times. Native files preserve original signal times and values. Never treat the two clocks as interchangeable. Previous-sample alignment avoids interpolation through dropouts. See column_dictionary.csv for source BlockPath, signal index and constant-tail assumptions.

## Metadata and limitations

A recording is not a complete training dataset. Complete session/calibration/source settings/payload metadata in the supplied manifest. Split training/validation/test by whole trials or sessions. Do not estimate actual pressure lag solely from command-to-motion data. Source MAT was not modified.

## Warnings

- Optional signal absent: sensor_host_s
- Optional signal absent: sensor_poll_start_host_s
- Optional signal absent: sensor_valid
- Optional signal absent: sensor_overrun
- Optional signal absent: base_orientation_valid
- Optional signal absent: length_inferred_m
- Optional signal absent: sensor_relative_xyz_m
- Optional signal absent: module_tip_xyz_m
- Optional signal absent: module_tip_valid
- Optional signal absent: R_NDI_from_arm
- Optional signal absent: gravity_arm_mps2
- Optional signal absent: xyz_configuration_valid
- Optional signal absent: xyz_fit_residual_m
- Optional signal absent: nonbending_rotation_rad
- Optional signal absent: length_velocity_filtered_mps
- Optional signal absent: length_acceleration_filtered_mps2
- Optional signal absent: ndi_raw1
- Optional signal absent: ndi_raw2
- Optional signal absent: ndi_raw3
- Optional signal absent: ndi_raw4
- Optional signal absent: ndi_raw5
- Optional signal absent: ndi_raw6
- Optional signal absent: ndi_raw7
- Optional signal absent: ndi_raw8
- Optional signal absent: p1_voltage_V
- Optional signal absent: p2_voltage_V
- Optional signal absent: p3_voltage_V
- Pre-valve pressure request not identified; final NI commands remain the action.
- Optional signal absent: kill_switch
- Optional signal absent: trialId
- Optional signal absent: phaseId
- Optional signal absent: channelId
- Optional signal absent: cycleId
- Optional signal absent: direction
- Optional signal absent: done
- Optional signal absent: testEnd_s
- Optional signal absent: measurement_time_sim_s
- No host acquisition timestamps: physical transition timing cannot be validated.
- Trial/phase labels absent: experiment type and settings must be supplied manually.
