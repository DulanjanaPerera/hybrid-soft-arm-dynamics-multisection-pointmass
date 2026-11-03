# Force-model preparation: two older days only

October19+October27, ten raw sources hash-checked against the interface audit.
training_data.mat stores copied measured-state/command/clock/quality sequences,
new whole-trial partitions, normalizer and corrected-nominal model configuration.
Legacy inverse-pressure and acceleration residual labels are not copied/used.
The original two_day_v2 files and normalizer remain unchanged.

Within each day/channel stratum, middle listed trial is internal selection;
other trials train.59,630 mature TRAIN rows fit the eight-feature normalizer.
2,883 deterministic training windows (24steps),160 selection windows(40steps),
10frame initial history. Existing validity/gap segments, physical timestamps,
causal derived velocities and warmup remain intact; no window crosses a trial,
partition or acquisition gap. Window/time/trial provenance is in CSV.

October27 StepWave contributes3 training windows/no selection; combined pressure
contributeszero windows due to acquisition gaps. Other eight recordings each
contribute360 training/20 selection windows. Do not claim all data were usable.

November2 is absent from this preparation, model-fitting and normalization.
It is evaluated only after release freeze. Inputs are q2,q3,v2,v3,p1,p2,p3,dt;
future q forms the integrated-output-error target only. Force labels are not
assumed measured. Reported dt/pressure alignment retains the earlier limitations.

The training-probe MAT verifies neural/solver gradient consistency and target-
free propagation. hybrid_probe.mat records the practical CPU runtime. The
abandoned full-autodiff compilation probe performed no parameter fitting.
