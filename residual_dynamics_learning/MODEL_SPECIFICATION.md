# Step 1: residual model interfaces and learning targets

Step3 implementation note2026-11-01: initial forward benchmarking uses a
position-transition residual with frozen-coefficient nominal steps to avoid
treating derived accelerations as ground truth. It is distinct from the
acceleration-residual template below; output units/injection are not
interchangeable. Full nonlinear nominal RHS rollouts are also benchmarked.
See benchmarks/screen_20261101_v2/REPORT.md for trained models and limitations.

Completed2026-11-01 as an offline interface specification, not trained-model
selection. See audits/interfaces_20261101 for ten-run inventory and provenance.

## Two separate acquisition cohorts

| Cohort | Raw directory | Acquisition dates in saved metadata | Runs |
|---|---|---|---|
| day1_root | RL_data_collection/ | 2026-10-19 | StepWave, SineWave, combinePressure, circleSweep, helix |
| day2_20261027 | RL_data_collection/20261027/ | 2026-10-27 | Same five filenames, different recordings |

Dates confirmed by SimulationMetadata.TimingInfo.WallClockTimestampStart, not
file modification times. Preserve the recorded timestamp without inventing a
timezone absent from its value. Day1 model name RL_dataCollection; day2
RL_dataCollection_v2. Pressure runs have77signals onday1 and82onday2; circle/
helix83bothdays. Run duration and signal sampling also differ. Each run ID
contains cohort plus stem; do not join or overwrite by filename alone.
cohort_provenance.csv contains hashes. No MAT data modified or merged.

## Physical system interface

- System input u=[P1_NI_DAQ;P2_NI_DAQ;P3_NI_DAQ], bar, from Pneumatic Valves
  provenance. All30 scalar channels located across ten runs. Commands are
  setpoint proxies, not measured chamber pressure or delivered force.
  Current saved V2 valve wiring confirms these taps are after Kill/saturation
  and before the10/9 pressure-to-voltage gain; P*_voltage is a separate role.
- Observed response: NDI-derived q=[l2;l3], metres; l1=-l2-l3. Primary measured
  signal lengthChnage_m at the reader is three channels; select2/3. Retain
  phiOrientation, thetaOrientation/observability and corrected XYZ as diagnostics.
- State x=[q;dq], four coordinates. dq is derived/filter-dependent, not a
  separate sensor. Retain physical poll/host and simulation/command timestamps.
- Required quality: reader validity, paired fresh frames, base orientation,
  finite values, elapsed times, Kill/limits. Do not drop straight configurations
  solely because theta is unobservable; inferred q can remain useful there.
- Exclude nominal des_theta/des_phi, comp_lengths, pressureUsed_bar and
  pressureDesired_bar from physical response/action labels. Their provenance
  identifies simulated forward outputs or inverse predictions. Reference
  theta_from_block/phi_from_block is a different role again.

Commands in step/sine/combined logs are typically sampled at.05s; circle/helix
NI commands are sampled at.01s. Never assume equal rates or align by array row.
Exact indices, native dimensions, times, units and ranges are in signals.csv.
Historical metadata/model identity is not proof of current saved SLX byte
identity. Extracted SLX scripts and saved-model hashes are current snapshots.

## Nominal mechanics and input map found in source

armS_single_entry -> armS_single_dynamics -> armS_single_core implements

    qdot = dq
    ddq_nom = M(q)^(-1) [tau_nom - (C(q,dq)+D)dq - G(q) - K_eff(q)q]

Current saved collection model forward blocks use L=.17411m, r=.013m,
mass=.1kg, K=650*[2 1;1 2], D=40*eye(2), g=[0;0;-9.81], mu=2000,
lengthBounds=[-.035;.035;1e6]. These are model constants, not identified truth
or verified payload metadata. The stiff input is unused in the extracted blocks.

    tau_nom = alpha * A * clamp(p_command,0,3)
    alpha = pi*(.013/2)^2*1e5 N/bar
    A = [-1 1 0; -1 0 1]

Active source has zero deadzone subtraction; the .8bar subtraction is commented.
No actuator memory model exists here. The saved inverse computes
tau_nom=M*ddq+(C+D)dq+G+K_eff*q, then chooses the smallest nonnegative pressure
offset and scales effective pressures if peak>3bar. V1 invalid inverse inputs
return zero; V2 holds last successful outputs with modelValid=false. Therefore
held or zero inverse outputs cannot become valid residual labels by themselves.
Step2 must verify the selected nominal implementation against recorded model
predictions before treating it as the historical baseline.
Also define whether inverse targets use unconstrained demanded pressure or
bounded nominal pressure; do not absorb allocator clipping as missing physics.

## Forward residual: missing dynamics prediction

Recommended contract:

    ddq_hat = ddq_nom(q,dq,p_command) + r_a_psi(H_k)

Output r_a has two entries in m/s^2. Integrate alongside qdot=dq. Use measured
q sequences for output-error fitting/validation when acceleration labels are
unreliable. This is an effective dynamics residual, not measured force.
Inputs H_k include available state/motion history, current/past commanded
pressures, intervals and known gravity/conditions. No future measured motion
or reference-tracking error is a forward prediction feature. A one-step label
can use future observed motion as TARGET, never as an available forward input.

## Inverse residual: three-pressure interface with explicit allocation

User-facing output remains three pressure commands:

    p_hat_k = p_nom_k + delta_p_psi(H_k, T_k, allocation_k)

H_k uses current/past measured q/dq and pressures only through k-1. Pressure
derivatives must also stop at k-1. T_k is the required motion/preview. During
inverse identification the recorded observed motion can define that motion
requirement; at execution preview must come from the planned desired trajectory,
not unknown future measurements. Validation must check that substitution.

Raw target p_command - p_nom is NOT yet ready for extraction: nominal inverse
must be evaluated on the same observed motion, time convention and allocation
policy, with valid derivatives/state history. Do not subtract a nominal pressure
computed from a different desired path, or blindly use logged inverse predictions.

Crucial allocation ambiguity: A*[1;1;1]=0. The nominal inverse's minimum-offset
choice differs from arbitrary recorded common pressure, so a naive three-output
residual absorbs allocation choice as well as missing dynamics.

Initial recommendation is an explicit common-pressure input p0 and two learned
differential residuals, then reconstruct the three-output inverse model:

    d=[p2-p1;p3-p1]
    p=p0*[1;1;1]+B*d
    B=[-1/3 -1/3; 2/3 -1/3; -1/3 2/3]

Nominal d comes from nominal inverse; learn delta_d. At deployment p0 is an
externally specified allocation choice. Using recorded p0 during training is
valid only for that conditional task; don't claim it predicts three independent
unknown commands without that input. Keep final bounds/rate handling explicit.
Alternatively a three-output residual needs a specified common-pressure policy
or additional allocation objective. Common pressure can affect the real arm's
stiffness/activation; the nominal nullspace does not prove physical irrelevance.

Forward and inverse residuals are distinct models/targets. A pressure-output
network alone does not supply an independent forward dynamics predictor.

## Architecture benchmark and feature policy

No winner selected before training. Benchmark on the same day-aware splits:

1. Nominal only / zero residual.
2. Simple linear history or low-order actuator-memory extension.
3. MLP with explicit delayed state/command windows.
4. Small GRU; consider LSTM only if held-out tests justify extra complexity.

Select history/preview duration, derivatives and capacity on development data.
Do not assume a memory duration from a filename. Start with q/dq/pressure and
timing; add pressure rates or angle features only by ablation. Theta features
use sin/cos or documented unwrap; q/angles are geometrically redundant.

Reserve cross-day testing: propose day1 development, day2 locked external
evaluation after checking protocol comparability. Hyperparameters need a
separate whole-trial/run validation split inside development, with sufficient
history/preview guard intervals. Do not tune repeatedly against day2 or treat
similar waveform repeats as independent adjacent random samples. Final fresh
experiments provide additional confirmation. Calibration/load/mounting equality
between days is not established by timestamps or signal names.

Evaluate forward rollouts and inverse required-command behavior, not training
loss alone. Physically equivalent pressure allocations need forward/motion
validation as well as pressure MAE. Freeze learning before matched nominal/
augmented PDFB and PAC comparisons. Residual learning is model augmentation,
not a separate error-to-command feedback policy.

References motivating benchmark structure (not evidence of an arm-specific winner):
MathWorks nonlinear ARX delayed regressors:
https://www.mathworks.com/help/ident/ug/what-are-nonlinear-arx-models.html
Lin et al., physics-guided GRU residual inverse identification (2025 preprint):
https://arxiv.org/abs/2507.14052

Next step2: prepare quality-controlled, aligned, day-aware training datasets
and verify nominal forward/inverse baselines. No training or hardware run yet.
