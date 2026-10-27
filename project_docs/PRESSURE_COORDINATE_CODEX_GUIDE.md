# Handoff guide: pressure-coordinate model validation, identification, and control

Prepared 2026-10-20 for the Codex agent working in `Dyn_Pressure/Matlab/dynamics`.

## 1. Purpose and evidence boundary

Use the length-coordinate project as a software and experiment-design reference, not as a source of pressure-coordinate equations or transferable gains. The intended sequence is: commission logging; build forward/inverse validation blocks; identify effective command-to-motion parameters; validate on held-out trials; implement PID computed-effort control and then PAC; write the paper from verified results.

This guide is based on the pressure project's seven `project_docs` files, its single-section RHS `arm_dynamics_pressure_2DoF_MSF.m`, `ROBOSOFT_2027_DATA_LOG_SPEC.md`, and `ROBOSOFT_2027_SINGLE_SECTION_PAPER_PLAN.md`. Reference implementation inspection covered the saved MATLAB Function scripts inside `RL_dataCollection.slx`, `armS_single_dynamics.m`, `commonPressureFeedbackInverseBlockV2.m`, and the fixed/adaptive PAC sources. No model, hardware, identification fit, or controller was executed for this guide. The blocks below are interface recommendations, not newly implemented pressure controllers.

The pressure documents currently establish numerical verification of selected mechanics/geometry/MEX cases. They do not establish physical command coupling, identified hardware parameters, or hardware controller validation. Their current state/decisions take precedence over older historical proposals. Chirp/least-squares identification below is a proposed workflow; this guide does not claim the length project completed that exact workflow.

Pressure-project root:
`C:/Users/dperera/OneDrive - Texas A&M University/Lab/Research/Controlling/Dyn_Pressure/Matlab/dynamics`

Reference-project root:
`C:/Users/dperera/OneDrive - Texas A&M University/Lab/Research/fDyn_ExpCom/hyDyn_Ctrl`

## 2. Resolve the coordinate/input contract first

The user's distinction is correct: pressure-valued coordinates do not make pressure the input to the mechanical equations. In this project they are **pressure-equivalent shape coordinates**, not chamber-pressure measurements.

| Quantity | Meaning | Size / unit |
| --- | --- | --- |
| `q_eq` | Independent mechanical shape coordinates | 2x1, Pa |
| `dq_eq`, `ddq_eq` | Shape-coordinate rate and acceleration | 2x1, Pa/s; Pa/s^2 |
| `Q` | Work-conjugate generalized effort | 2x1, J/Pa |
| `p_cmd` | Three physical regulator pressure setpoints | 3x1, Pa in analysis; optionally bar at hardware boundary |
| `p_ch` | Actual chamber pressures, currently unmeasured | 3x1, Pa |
| `theta`, `phi` | Bending-plane direction and bend magnitude | scalars, rad |

Mechanical equations:

`M(q_eq)*ddq_eq + C(q_eq,dq_eq)*dq_eq + G(q_eq) + F_el(q_eq) + D*dq_eq = Q`.

`Q'*dq_eq` has units W. Thus `M` has units J*s^2/Pa^2, `C` and `D` have units J*s/Pa^2, and a linear elastic coefficient has units J/Pa^2. A command-coupling matrix mapping Pa-valued setpoints to Q has units J/Pa^2. Length-coordinate gains and forces cannot be copied numerically.

Use the physical chain explicitly:

`p_cmd -> effective actuation/lag model -> Q -> mechanical forward dynamics -> q_eq -> measured-pose prediction`.

Inverse control reverses the mechanical stage, then separately allocates effort:

`desired motion + measured state -> Q_requested -> feasible command allocation -> p_cmd`.

Never connect `Q_requested` directly to the regulators. A two-coordinate chart also does not imply two physical chambers or require three actual pressures to sum to zero. Negative chart coordinates may be valid while all regulator commands remain nonnegative.

### Coordinate transformation and the unresolved K issue

The pressure project's documented length-change chart is

`ell = c*[-q2-q3; q2; q3]`, `c=3*A*r^2/(2*K)`.

For independent lengths `ell_red=[ell2;ell3]`, `ell_red=J*q_eq` with constant `J=c*I` only when the chart parameters are fixed. Then

`M_q=J'*M_ell*J`, `C_q=J'*C_ell*J`, `G_q=J'*G_ell`, `D_q=J'*D_ell*J`, and `Q=J'*tau_ell`.

These relations describe a consistent constant coordinate change; configuration-dependent charts require the extra Jacobian time-derivative terms. Reconcile gravity and constitutive conventions before claiming the existing two implementations are dynamically equivalent. The pressure documents flag opposite historical gravity signs.

The pressure geometry formula needs K in N*m, while the inherited expression `(3/20)*k*L*r^2` has N*m^2 if k is N/m. Resolve this discrepancy before reporting physical stiffness or freezing the chart. If K is an optimization variable, recompute the chart/output map for each candidate; do not transform measured angles once with K0 and then optimize K independently. Prefer freezing a dimensionally justified chart and fitting separate effective elastic parameters.

Do not adapt the coordinate-chart K online as though it were merely an elastic gain. Changing the chart changes measured coordinates, rates, inertia, and reference conversion. Online elastic adaptation should operate in a fixed chart unless a full time-varying coordinate derivation is provided.

## 3. What the length-coordinate blocks actually do

The saved `RL_dataCollection.slx` embeds MATLAB Function code; standalone files with similar names are not automatically the active saved code. Inspect and compare the chart scripts before asserting that a model is updated.

### Forward block

Its current `forwardModel(pressureCommand_bar,X,stiff)` accepts pressure commands and the externally integrated state `X=[ell2;ell3;dell2;dell3]`. It clamps commands to 0..3 bar, maps them to two forces with

`tau = area*1e5*[-1 1 0;-1 0 1]*pressureUsed_bar`,

and calls `armS_single_entry` to return the state derivative. Simulink's external integrator owns the mechanical state. It also outputs predicted lengths, theta/phi, nominal force, and theta observability. In the inspected saved model the 0.8-bar subtraction is commented out, and the `stiff` port is unused: K is hardcoded to `650*[2 1;1 2]`, D to `40*I`, and gravity to `[0;0;-9.81]`.

### Inverse block

`inverseModel(l,orientationValid,stiff)` receives six desired or supplied motion quantities `[ell2;ell3;dell2;dell3;ddell2;ddell3]`, evaluates the same mechanical model and smooth length-bound penalties, and computes

`tau_requested=M*ddell+(C+D)*dell+G+Keff*ell`.

It then applies a separate minimal-common nonnegative pressure allocation, a 3-bar ceiling, and logs requested versus nominal available force. Its current deadzone is 0.0 bar and `stiff` is also unused. The inverse code derives effective pressure area from `params.r/2`; the forward code has a separately named pressure diameter. Their numerical equality is specific to this hardware. Keep backbone-to-PMA radius and pressure-area radius separate in the new project.

An inverse block fed only desired motion is feedforward. It becomes computed-torque feedback only when measured state and a feedback virtual acceleration enter the mechanical effort calculation. Logging a measured validity flag alone does not turn feedforward into closed-loop motion control.

### Reuse architecture, not constants

Reuse: shared dynamics, external integration, fixed-size interfaces, explicit SI conversion, causal estimation, validity gates, force diagnostics, separate allocation/supervision, and offline replay. Do not reuse: 650/40 mechanical values, 0.8-bar activation assumption, pressure trims, area assumptions, gravity signs, or length-coordinate feedback/adaptation scales without derivation and data.

## 4. Recommended code/block organization

Keep the numerical core pure: no serial calls, filesystem access, plots, workspace lookups, or hidden hardware commands. Use a thin MATLAB Function wrapper with `%#codegen`, fixed-size outputs and typed Boolean flags. Put tunable values in one configuration/tuning block or a defined bus; avoid dozens of scattered constants and duplicated parameter tables.

Suggested functions:

| Function | Responsibility |
| --- | --- |
| `pressureMechanicalTerms(q,dq,params)` | Shared M, C, gravity, elastic effort, damping, validity |
| `pressureForwardBlock(q,dq,Q,params)` | Return ddq and term diagnostics; external integrators own q,dq |
| `pressureInverseEffortBlock(q,dq,a,params)` | Return Q=M*a+C*dq+G+F_el+D*dq and validity |
| `commandEffortBlock(p_cmd,q,dq,actuatorState,params)` | Identified effective command-to-effort model; optional separate lag state |
| `effortAllocationBlock(Qrequested,q,dq,params)` | Feasible 3-command allocation, equivalent effort, residual, limits |
| `pressureStateEstimatorBlock(...)` | Pose-to-fixed-chart conversion and fresh-sample causal rates |
| `pressureReferenceBlock(t,...)` | Smooth q reference with analytic first/second derivatives |
| `pressureSupervisorBlock(...)` | Enable, blend, sample age, fault policy, command bounds/slew |
| `pressureTuningBlock()` | Named fixed-size gains/options and before-run structural parameters |

Use unambiguous names such as `QRequested_JPerPa`, `QEquivalent_JPerPa`, `allocationResidual_JPerPa`, `qEquivalent_Pa`, `thetaMeasured_rad`, `thetaReference_rad`, `thetaPredicted_rad`, and `pressureCommand_bar`. The pressure project's older `tau` and `G_and_K` names need wrapper clarification: inspect whether a helper already includes elasticity before adding it again.

Mechanical core pseudocode:

```matlab
% Forward: q,dq,Q have the declared shape-coordinate units.
[M,C,G,Fel,D,valid] = pressureMechanicalTerms(q,dq,params);
rhs = Q - C*dq - G - Fel - D*dq;
ddq = M\rhs;

% Inverse uses exactly the same terms and parameter snapshot.
Qrequested = M*a + C*dq + G + Fel + D*dq;
```

This is a contract, not ready-to-paste code: implement finite-input, mass-matrix, size and failure checks. Never use `inv(M)`; check symmetry/positive definiteness with scale-aware criteria. Do not silently turn NaN/Inf acceleration into zero and report success. The existing single-section pressure RHS does zero nonfinite acceleration; revise that policy explicitly for validation/control wrappers, exposing failure and terminating offline fits when appropriate. Any regularization must be logged and its effect checked.

### Simulink layouts

1. **Mechanical forward/inverse harness:** synthetic q,dq,ddq -> inverse effort -> forward acceleration; compare recovered ddq. No hardware blocks. This tests algebra, not command physics.
2. **Open-loop validation model:** final three commands -> effective actuator model -> mechanical forward block + integrators -> pose prediction; NDI -> measured pose/logging in parallel. Log both without feeding prediction into measurement.
3. **Inverse feedforward validation:** reference -> desired q,dq,ddq -> inverse effort -> allocator/supervisor -> final commands; measured motion is an evaluation output. Include a parallel forward replay branch if useful.
4. **Controller model:** desired reference and measured state -> PID/PAC effort -> allocator/supervisor -> final commands. Keep Kill and device interfaces downstream and explicit.

These can share one acquisition/logging subsystem. Their difference is more than which input wire changes: state ownership, feedback route, actuation map and validity semantics must remain explicit. Avoid an algebraic loop through a memoryless inverse, direct-feedthrough actuator model and measured-state calculation.

## 5. Sampling, derivative estimation and NDI handling

The reference reader is `AuroraReader_s2_BaseRotation.m`, supported by `AuroraDriver_2.m`, `PortHandle.m`, `ndiSupportFnc.m` and calibration data. It provides paired XYZ, orientation-based theta/phi, inferred lengths, frame/host-time diagnostics, base attitude/gravity and validity flags. Inspect the pressure project's actual reader before assuming identical ports or mounting/calibration.

- Verify compiled controller sample time. Changing Constant block sample times alone does not prove the MATLAB Function executes at that period. The reference system has used 0.05-s sensing/controller updates and faster NI output tasks; these are historical settings, not requirements for the new arm.
- Separate simulation scheduling time, fresh-sensor sample time, host acquisition time and actual command-dispatch time. Host polling overruns must not be hidden by uniform simulation timestamps.
- Update estimation/adaptation once per new valid frame, using a bounded positive dt. Repeated values at faster downstream task rates are not new measurements. Handle clock reset/new run and long gaps explicitly.
- Desired trajectories should supply analytic rates/accelerations when possible. The inspected length reference filter uses backward differences plus first-order velocity/acceleration filters with 0.10/0.15-s constants; these cause lag and are not universal settings.
- Measured derivative filters must be causal for control, with their own documented initialization, delay and invalid-frame policy. Offline acausal smoothing is allowed for analysis copies but must not masquerade as deployable estimation.
- Near straight shape theta is poorly observable. Use smooth bend components `u=phi*cos(theta)`, `v=phi*sin(theta)` for estimation/fitting. Gate plane-angle metrics using measured bend/observability; do not punish arbitrary theta at zero bend.
- Wrapped theta error is `atan2(sin(theta_des-theta_meas),cos(theta_des-theta_meas))`. Phi is the bend magnitude, normally compared directly on the supported branch.
- Distinguish raw sensor-relative XYZ from offset-corrected module-tip XYZ and identify the frame. The reference `moduleTipXYZBlock.m` removes axial 55-mm base and 50-mm tip offsets using current base attitude and measured tip direction. Those mounting offsets are hardware-specific.
- Held last-valid displays are stale data. Carry valid/fresh/age flags separately; a stale value must not silently enter fitting or adaptation as a fresh observation.

Keep NDI setup/preflight independent from control profiling. In this project serial polling slowed severely until MATLAB was restarted; inherited controller rate changes did not explain a preflight polling fault. Profile acquisition separately from pure model computation.

## 6. Identification: single-PMA and combined-PMA chirps

First inspect one short commissioning recording. Confirm pressure channel order, polarity, units, frame convention, final post-saturation commands, clocks, validity, and pose geometry. Do not start a long fit campaign using assumed signal names.

Plan repeated loading/unloading static holds and low-frequency tests before chirps. For each PMA, excite around an established feasible bias, keeping the other two commands explicitly specified. Then excite common mode, differential modes, and independent/multisine combinations. Under a simplified difference-input law, equal commands are a nominal null direction; physical common pressure may still change preload, extension, stiffness or hysteresis. Record it and test its influence rather than discarding it.

A smooth chirp may use `p_i(t)=bias_i+amp_i*w(t)*sin(psi(t))`, with `psi=2*pi*(f0*t+0.5*(f1-f0)*t^2/T)` for a linear-frequency chirp. Distinguish Hz from rad/s. Choose bias, amplitude, duration and frequency endpoints from hardware limits and commissioning evidence; this guide supplies no approved excitation limits. Check the whole combination against rails and slew, including chirp phase derivatives and window ramps. Use labels for channel, mode, repeat, direction and loading/unloading phase. Repeated trials must preserve initial/rest conditions and session/calibration identity.

### Bounded least-squares pipeline

1. Freeze whole-trial train/validation/test splits before fitting. Do not randomly split adjacent samples of the same chirp.
2. Start with geometry, mass and coordinate chart fixed where independently justified. Select a small effective parameter vector: elastic terms, positive damping, identifiable input gains, and only evidence-supported deadband/lag/asymmetry. Avoid simultaneously freeing every scale.
3. Forward-simulate each recorded final command sequence using zero-order hold and actual timing assumptions. Initialize from measured shape and defensible rate estimates; estimate initial states as constrained nuisance variables if needed.
4. Construct fixed-length residuals from valid measured angles or bend components. For plane angle use wrapped residuals and exclude unobservable theta. Normalize residual units/noise scales. Add XYZ only with declared weighting/correlation, since it may share NDI sensing with angles.
5. Use `lsqnonlin` (or an equivalent bounded solver) with parameter scaling and multiple starts. Positive damping and positive-definite elasticity should be enforced structurally, e.g. Cholesky parameters, rather than repaired after the fit.
6. Failed integration or invalid model states must yield a deterministic finite penalty vector of the same dimension; never silently remove difficult samples or return zero residuals.
7. Inspect residuals versus channel, common command, motion direction, amplitude, speed and time. Inspect sensitivity/Jacobian rank, parameter correlation and bound hits. Static data usually identify input-gain/stiffness ratios; a good fit alone does not prove unique physical parameters.
8. Compare held-out free-run prediction against a static command-to-shape map and a simple dynamic baseline. Report transient/steady error, phase lag, XYZ error, repeat variability and simulation failures. One-step prediction that repeatedly resets to measured state is a different task and must be labeled.

Do not differentiate noisy measured pose twice to manufacture the sole least-squares target. Forward output-error fitting avoids making acceleration noise the main regression signal. A separate inverse-effort residual can be a diagnostic if derivative uncertainty is quantified.

Without chamber sensors, label identified dynamics as **effective command-to-motion** behavior. An identified lag cannot uniquely be called regulator lag, and command setpoints must not be renamed measured chamber pressure.

## 7. PID feedback linearization / computed-effort control

For fixed pressure-equivalent shape coordinates, use

`e=qd-q`, `de=dqd-dq`, `zdot=e`,

`a=ddqd+Kd*de+Kp*e+Ki*z`,

`Qrequested=M(q)*a+C(q,dq)*dq+G(q)+Fel(q)+D*dq`.

This evaluates cancellation terms at measured state. Reference-only inverse dynamics is feedforward, a different baseline. Gain units are Kp: 1/s^2, Kd: 1/s, Ki: 1/s^3 in this virtual-acceleration form. An outer effort-PID formulation has different units; do not mix them. For ideal diagonal independent axes, positive gains and `Kd_i*Kp_i>Ki_i` give the cubic stability condition; this does not certify the delayed, saturated physical loop.

Start without integral action and with conservative bandwidth relative to sensing/actuation delay. Add integral only after feasible tracking is established. Bound integral states and use conditional integration or justified back-calculation when allocation/slew/fault prevents requested effort. Keep a reference/enable blend and reset policy; resetting the integral arbitrarily on every reference change can create steps.

Allocate using the same identified command map as the forward model. Log Qrequested, Qequivalent under final commands, residual, active constraints and common pressure. If rate limiting changes the command, recompute the equivalent effort/residual after that stage; allocator feasibility before slew does not mean final effort was delivered. With unmeasured chamber pressure, Qequivalent remains a model prediction.

## 8. PAC and adaptive extension

Use PAC only after the fixed-parameter effort loop works. For `e=q-qd`, define

`v=dqd-Lambda*e`, `a=ddqd-Lambda*(dq-dqd)`, `s=dq-v`.

Two consistent damping treatments are possible:

- Fixed-controller reference implementation: `Q=M*a+C*v+G+Fel+D*dq-Ks*s`.
- Adaptive damping form used in the adaptive reference: `Q=M*a+C*v+G+Fel_hat+Dhat*v+Hhat+bhat-Ks*s`.

Do not accidentally replace `D*dq` by `D*v` while assuming identical force or proof. Derive the selected closed-loop s dynamics and adaptation regressor together. For pressure coordinates Lambda has units 1/s and Ks has units J*s/Pa^2. It is distinct from physical/model damping D. Neither Lambda nor Ks is automatically adapted by the parameter-update law.

For linear-in-parameter terms `Y*theta_hat`, a baseline update is `theta_hat_dot=-Gamma*Y'*s`, with signs matched to the declared s convention. The reference implements normalization, projection bounds, rate caps, deadband and bounded fresh-sample dt; its normalized rate is `-(scale.^2).*gamma.*(Y.'*s)`. Rebuild scales, rates and units for the new chart. Keep mass/inertia and geometry fixed initially; do not force nonlinear uncertain terms into an invalid linear regressor.

Expose separately: pressure enable, parameter-adaptation enable, common-pressure policy enable, and optional hysteresis-weight freeze. Disabling parameter updates should hold estimates; it does not necessarily remove their force contributions or freeze internal history states. Log which behavior is actually implemented.

Gate updates on valid/fresh measured state, valid dynamics, active control phase and acceptably small final allocation residual. Frozen parameters can still yield an unstable delayed/saturated loop; turning adaptation off does not isolate every controller mechanism. Any projection/rate/deadband/constraint modification also needs its own proof qualifications rather than inheriting an ideal adaptive theorem unchanged.

The reference's switchable common-pressure logic is primarily **constrained allocation**, not an identified adaptive pneumatic constitutive model: it projects the requested two-force vector onto feasible pressure differences and chooses a feasible mean near a nominal target. Changing physical common pressure may change the plant. Name and report this distinction.

A k3 bound of zero in the structured length stiffness is not by itself proof of the observed oscillation or permission to allow negative stiffness. In a pressure-coordinate implementation enforce the actual elastic matrix/energy constraints. Hysteresis weights and bias can absorb modeling errors; their numerical drift does not prove physical material identification.

Practical evidence from this arm: reducing external PAC Ks from 40 to 20 greatly reduced limiting/oscillation in one static protocol; reducing Lambda also helped, with slower settling. This is evidence about gain/delay/allocation interaction in that experiment, not a recommended numerical pressure-coordinate gain. The earlier high-gain runs showed large allocation residuals and invalid sensor fractions during oscillations. Analyze time-aligned requested/delivered effort, velocity, limits, sample age and adaptation gates before attributing cause solely to hysteresis.

## 9. Logging and validation gates

Retain original MAT/logsout plus non-overwriting exports. Resolve duplicate names by physical block-path provenance, not by silently selecting the first match. The reference exporter `RL_data_collection/exportRlTrainingData.m` distinguishes final NI commands, requested feedforward pressure, measured NDI angles, requested circle/helix angles, and forward-model predictions. Adapt those semantics, not just column spelling.

Minimum records: original times; final three commands; pre-limit requests; measured theta/phi; corrected and raw XYZ with frames/validity; inferred lengths; q_eq and its fixed chart parameters; estimated rates; frame freshness; acquisition timestamps/overruns; reference and experiment labels; base orientation/gravity; model/policy/configuration versions; fit parameters; calibration identity; enables/faults; predicted motion; requested/equivalent generalized efforts and residuals; all clipping/slew/activation indicators. Record inferred quantities as inferred. No pressure-sensor channel exists; keep that absence explicit.

Acceptance sequence:

1. Unit/frame/geometry checks, straight-pose limits, finite dimensions and sensor-offset checks.
2. M symmetry/positive definiteness, potential-gradient gravity sign, elastic energy and passive dissipation checks.
3. Shared-core forward/inverse round trip over the actual operating envelope; direct and compiled parity where used.
4. Allocation round trip for feasible requests, meaningful residuals for infeasible requests, rails/slew and mode transitions.
5. Fresh-frame/repeated-time/invalid/gap/reset tests; Boolean port typing and persistent-state reset.
6. Offline replay at intended sample time and measured delay; held-out prediction and matched baseline comparisons.
7. Hardware commissioning by the researcher, followed by controlled comparisons only after the earlier gates pass.

Do not run hardware-bearing Simulink models merely to compile-check a helper. Loading, updating or simulating a model can trigger callbacks or device setup; inspect these first and use hardware-free harnesses for automated tests. Do not modify the user's SLX wiring, deploy a controller, commit or push unless authorized. Keep a source snapshot/hash of embedded chart code so the paper can identify what actually ran.

## 10. RoboSoft manuscript workflow

Follow the pressure project's current single-section plan. This guide does not revalidate submission deadlines, page limits, anonymity or AI policy; check current official venue instructions when preparing submission.

Draft hardware, notation, assumptions and reproducible methods early. Add numeric claims only after frozen analysis. Suggested evidence chain:

1. Explain why Pa-valued shape coordinates and generalized effort are useful, with work/units and the physical input bridge made explicit.
2. Describe hardware channels, absent chamber sensing, coordinate chart, sensor offsets, command/pose clocks and operating envelope.
3. Report excitation protocols, repeated trials, whole-trial split, parameter bounds/scales, objective, solver and identifiability limits.
4. Show held-out dynamic predictions and simple baselines, with command traces, repeat variability and error metrics.
5. If justified, compare matched static feedforward, feedback-only/PID and PAC runs using the same references, constraints, initialization and load. Separate adaptation ON/OFF and common-pressure policy effects.
6. Discuss hysteresis, common-pressure dependence, unmeasured pressure, estimator delay, chart/K ambiguity and numerical approximation limits.

Numerical energy/MEX checks are implementation evidence, not experimental model validation. XYZ and angles from one tracker are correlated; do not present them as independent sensing validation. A bounded adaptive estimate is not automatically an identified physical parameter. An ideal effort-level passivity argument is not automatically a regulator-command-level stability guarantee.

Maintain a result manifest mapping each figure/table to raw trial files, code/model version, parameters, processing settings and metric windows. Keep exploratory fitting plots separate from held-out results. Compile and visually inspect the final IEEE PDF, check cross-references/citations/figure units, and verify claims against the manifest. The sibling Paper repository/Overleaf branch is a separate user-controlled destination; inspect its current state before any authorized manuscript update.

## 11. Concrete first assignment for the receiving Codex

Read `project_docs/PROJECT_STATE.md`, `MODEL_ARCHITECTURE.md`, `DECISIONS.md`, `OPEN_ISSUES.md` and this guide. Then inspect the single-section mechanics and the saved `NDI.slx` / `Pressure_actuationOnly.slx` interfaces without operating hardware. Produce:

1. A signed/unit-checked table of q_eq, Q, p_cmd, theta/phi, XYZ and all coordinate conversions, including the unresolved K definition.
2. One shared mechanical-terms function and hardware-free forward/inverse harness; preserve existing functions until verified replacements pass.
3. A proposed effective command map and clearly labeled allocation interface; identify which parameters/data are still missing instead of inventing physical pressure sensing.
4. A commissioning logging audit and one short researcher-recorded trial inspection before selecting chirp frequencies or fitting parameters.
5. A small, reviewable implementation plan for bounded output-error fitting and held-out validation. Controllers follow validated coupling and timing, not the other way around.

After significant work, update current state and specialized docs and append the pressure project's history. Record what was actually tested, which artifacts changed, and which physical assumptions remain unresolved.
