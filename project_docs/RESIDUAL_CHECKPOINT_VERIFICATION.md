# Residual checkpoint verification

2026-11-03. Reran current-source tests with two MATLAB computation threads:
validateArmSBounds, validateForceResidualPhysics and validateNdiDiagnostics.
NDI tests use an injected byte-stream mock and synthetic logged signals; no
serial/DAQ connection or Simulink model execution. Full local transcript:
tmp/git_checkpoint_20261103_v1/verification.txt.

- PASS: legacy/corrected nominal RHS equivalence, inverse roundtrip and allocation identity, including implied-l1 bound.
- PASS: zero-lag original nonlinear RHS, no future-state leakage, delayed-input causality.
- PASS: three-actuator bounds, potential gradient, cyclic symmetry, legacy RHS and inverse/PAC consistency.
- PASS: exact nominal coefficients (1.55431e-15), force RHS sign/units, integration gradients (2.00496e-09), nonlinear ODE convergence.
- PASS:32x2 matrix ordering, named base/tip columns and legacy export unchanged.
- PASS:32 reader replay cases, legacy outputs, same-BX status, incomplete/recovery and named export. No hardware/model run.

Latest capacity/fine-tuning checks remain in their versioned VERIFICATION.md
reports. All448 capacity artifacts matched their frozen hashes before this
checkpoint; current64-unit reference and128-unit candidate hashes unchanged.
No training was repeated for this Git operation.
