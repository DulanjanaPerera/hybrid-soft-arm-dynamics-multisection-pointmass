# Fixed-parameter passivity controller and isolated pressure command path

Isolated implementation for q=[l2;l3] in metres; l1=-l2-l3. No payload,
parameter adaptation, integral action or hardware I/O. Stage 1 is force only;
stage 2 adds an isolated common-pressure allocator and command supervisor.
Existing models and dynamics functions are unchanged.

## Controller

With e=q-qd, dqr=dqd-Lambda*e, ddqr=ddqd-Lambda*(dq-dqd), s=dq-dqr:

    tau = M*ddqr + C*dqr + D*dq + G + Keff*q - Ks*s

`armS_single_core` supplies the existing nominal distributed-mass M/C/G;
the elastic force includes the same smooth bound penalty as
`armS_single_dynamics`. For an exact model and delivered requested force,
M*sdot+C*s+Ks*s=0. For V=s'*M*s/2 and skew-symmetric Mdot-2C,
Vdot=-s'*Ks*s. This is the ideal mechanical-model result, not a stability
claim for the sampled, filtered, pressure-limited physical system.

The core preserves the existing dynamics coordinate implementation (including
its internal l=[0,q'] representation); it does not rederive its kinematics.
Gravity uses the existing model's sign convention and must be checked against
physical orientation before hardware integration. Large deformations outside
the tested grid are not validated by the supplied checks.

## Files and use

- `pacNominalParameters.m`: geometry, nominal K=650*[2 1;1 2], D=40*I,
  length-bound penalty, Lambda=[2;2] 1/s and Ks=[10;10] N*s/m. Gains are
  offline starting values, not experimentally tuned replacements for PID gains.
- `pacFixedController.m`: stateless force controller, parameterized for checks.
- `pacFixedBlock.m`: paste-ready wrapper for a NEW MATLAB Function block.
- `validatePacFixed.m`: dynamics identities, invalid-data handling, gravity
  direction and hardware-free matched-model integration. Writes no files.

From the project root:

```matlab
addpath('passivity_adaptive_control');
report = validatePacFixed();
```

Block inputs are measured (4x1: q,dq), desired (6x1: qd,dqd,ddqd),
gravity_arm (3x1), measurementValid and baseValid (scalar logicals).
Connect fresh measured/desired rates from the existing filters. Measured
acceleration is not used. The wrapper has ten outputs: requested force,
s, dqr, ddqr, inertia, Coriolis, damping, static, feedback force, and validity.
All vectors except gravity are 2x1. Requested/feedback/component forces are N.

Do not feed force outputs directly to NI or the pressure-command input. Invalid
inputs return zero and valid=false, which must trigger the downstream
supervisor; zero generalized force is not zero pressure. The stateless wrapper
does not detect stale held samples by itself. It must run at an explicitly
chosen discrete controller sample time when later integrated. MATLAB Function
compilation in a Simulink model and physical pressure-path integration are not
validated by a successful MATLAB numerical check. Stage-2 wrapper C source
generation has separately passed MATLAB Coder; binary compilation was not run.

## Stage 2: common-pressure wrapper

`pacCommonPressureBlock.m` is the combined paste-ready wrapper for a NEW
MATLAB Function block, with `pacPressureConfig.m` holding explicit settings.
`pacAllocatePressure.m` and `pacPressureStep.m` separate the algebra from
stateful supervision; `pacPressureInitialState.m` creates reproducible state.
`validatePacPressure.m` tests these paths without Simulink or hardware:

```matlab
addpath('passivity_adaptive_control');
report = validatePacPressure();
```

The nominal trim is p0=[1.665599;1.649182;1.485219] bar. Capture the
nominal static model force tau0 at the settled measured trim configuration.
Let Fhat(p)=f*[max(p2-0.8,0)-max(p1-0.8,0);
max(p3-0.8,0)-max(p1-0.8,0)], f=13.2732 N/bar. Define the constant bias

    b = Fhat(p0) - tau0

The allocator requests pressure force tau_req+b, or equivalently adds the
increment tau_req-tau0 to trim. Thus tau_req=tau0 reproduces p0, while
changes in current gravity remain in tau_req. This is a local calibrated
offset, not a globally identified force map. The existing nominal model and
trim need not be physically correct away from calibration.

With ds=(deltaTau1+deltaTau2)/(2*f), da=(deltaTau2-deltaTau1)/(2*f),
deltaTau=tau_req-tau0, pressures are computed from s0+ds,a0+da and mean c.
The algebra is exact only in the equal-gain all-active region before limits.
The nominal threshold map exposes the force residual outside that region.
Log tauEquivalent=Fhat(p_command)-b, and
allocationResidual=tauEquivalent-tauBlended. During engagement tauBlended
interpolates from tau0 to tauRequested; after engagement they are equal.
Residual/equivalent outputs are NaN before calibration; residual is NaN in
hold, fault and disabled states. These outputs are not chamber measurements.

Even under the calibrated nominal force model, the mechanical energy balance
after full engagement becomes Vdot=-s'*Ks*s+s'*allocationResidual. A nonzero
allocation residual therefore cannot be ignored in a passivity claim. Actual
force-map mismatch adds another residual. The blend also changes the request.

### Inputs and outputs

Configure the block for discrete 0.05 s execution, not continuous/inherited
minor-step updates. Inputs, in order:

1. measured: 4x1 [l2;l3;dl2;dl3], metres and metres/second.
2. desired: 6x1 [qd;dqd;ddqd], consistent filtered reference derivatives.
3. gravity_arm: 3x1, existing calibrated convention, m/s^2.
4. measurementValid: scalar logical, actual paired-pose validity.
5. baseValid: scalar logical.
6. measurementTime_s: timestamp of last successful paired measurement, not
   Clock assigned to held data. Verify the sensor's actual timing signal.
7. t: Clock seconds.
8. enablePressure: scalar logical.

All seventeen outputs are named in the wrapper signature. Log pressureCommand,
tauRequested, tauEquivalent, allocationResidual, s, the five force components,
modelValid, limited, status, allAboveActivation, anchorReady, forceBias and
tauBlended. Force/s vectors are 2x1; pressure is 3x1; flags/status are scalars.

Keep the existing physical Kill and downstream NI bounds. This wrapper
replaces the old PID virtual acceleration/inverse pressure block path in a
future duplicated model; do not cascade it through the old inverse block.
There is no integral block in this path. Parameter helpers contain compile-time
settings; changing them requires recompiling/updating the model.

### Timing and limitations

Defaults preserve saved base-rotation values: mean 1.6 bar, s/a corrections
bounded to +/-2.5 bar, rails 0.2-3 bar, slew 2 bar/s. The 0.2 bar rail is below
the nominal 0.8 bar activation point; allAboveActivation explicitly reports
this. Keeping the old rail is not a claim of guaranteed deadzone avoidance.

- Zero command until 10 s; then slew toward the trim.
- Anchor capture at/after 15 s requires actual command at trim, valid/fresh
  input, and both measured speeds <=0.5 mm/s. Otherwise remain in trim.
- Feedback blends over five seconds from the actual anchor capture time.
  Keep the desired trajectory in its startup hold until status=3. A fixed
  20 s start is only appropriate if calibration actually completes at 15 s.
- Samples older than 0.075 s are rejected. Brief invalidity holds command up
  to 0.10 s after the last usable update; sustained invalidity latches and
  slews pressure to zero. Recovery of validity cannot clear the latch.
- Disable slews to zero and clears calibration. Re-enable repeats trim.
  Stop time defaults to Inf: end by disabling or the external experiment gate.
- Invalid Clock immediately commands zero; reversed Clock latches a fault.
  Clear persistent state only when starting a new run with output disabled.

`modelValid` describes force-model inputs only; it does not imply the
supervisor is engaged or fault-free. Check status and anchorReady as well.
Status codes: 0 off, 1 trim/waiting, 2 blend, 3 control, 4 invalid hold, 5 fault.

## Remaining integration gate

The force/trim relationship is now explicit, but its constant-offset accuracy
needs experimental evaluation. No existing SLX was edited. Check the sensor
timestamp wiring, discrete rate, gravity sign, external Kill and post-Kill
logs in a future duplicated model before calling this hardware-ready.

Reference: Azizkhani et al., Dynamic Control of Soft Robotic Arm: An Experimental
Study, DOI 10.1109/LRA.2023.3243802, https://arxiv.org/abs/2209.11332.


## Adaptive implementation (2026-10-11)

The earlier sections describe the preserved nominal stage. The new adaptive
path is pacAdaptiveCommonPressureBlock.m; see [ADAPTIVE_INTEGRATION.md](ADAPTIVE_INTEGRATION.md)
for all inputs/outputs, settings, scope, tests and limitations. K/D/history-force
weights/bias adapt; no mass or internal Bouc-Wen coefficient adaptation.

Workspace use: analyzePacAdaptive() or analyzePacAdaptive(out). Workspace calls print/plot without writing report files; MAT-path input retains report saving.

For switchable common-pressure allocation and external tuning inputs, see TUNABLE_BLOCK.md and pacAdaptiveTunableBlock.m. Run pacTunableDefaults to populate Constant-block values; validatePacTunable is offline only.
