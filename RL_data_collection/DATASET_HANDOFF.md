# Real-arm dataset handoff

Read the recording's exported REPORT.md, trial_manifest.json and
column_dictionary.csv before using its CSVs. Keep the source MAT as authority.

- Complete dynamics learning uses final NI pressure COMMANDS as actions and
  measured orientation/position as observations. Pressure is not independently
  measured; actuator dynamics cannot be separately identified without assumptions.
- des_theta/des_phi, comp_lengths, pressureDesired_bar and pressureUsed_bar are
  model predictions/transforms in these validation recordings, not real targets.
- observations.csv is sensor-aligned in simulation time. Actual host intervals
  vary; use sensor_dt_host_s and transition_from_previous_valid and respect trial
  boundaries. Actual NI-write wall timestamps are absent.
- Native CSVs preserve original logged timestamps and invalid/held samples.
  Duplicate signal names are indexed and have block-path provenance.
- Lengths are inferred from NDI orientation; module-tip XYZ combines measured
  sensor XYZ, orientation and assumed axial offsets (.055/.050 m). It is relative
  to the current module-base frame. It is not independent strain sensing.
- Verify trial/session/source parameters/calibration/mounting/payload metadata.
  Export does not invent absent metadata or establish physical safe bounds.
- Evaluate on held-out trials/sessions; do not randomly mix neighboring samples.

Initial exported full recording: exports/SineWave_reviewed/.
It contains six sequential trials (P1,P1,P2,P2,P3,P3), three cycles each,
636 simulation seconds,12721sensorrows,12710validorientationrows. Hostelapsed
648.3145s,12gaps>.10s,maxgap4.499s,17overruns. Inspect P3's larger measured
phi excursions (~2.49-2.51rad) rather than assuming equal channels have equal
motion. These are measured excursions, not commanded safe angle limits.
