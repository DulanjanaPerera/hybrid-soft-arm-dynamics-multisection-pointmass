# Adaptive PAC integration

This is a new controller path. The existing nominal files and SLX models are
preserved. Add this folder and the project root to MATLAB path. Paste the
contents of `pacAdaptiveCommonPressureBlock.m` into a new MATLAB Function
block in your duplicated PAC model; use **Discrete, 0.05 s** execution.

## Inputs, in order

1. measured: double 4x1 `[l2;l3;dl2;dl3]`, metres and metres/second.
2. desired: double 6x1 `[qd;dqd;ddqd]`.
3. gravity_arm: double 3x1 from the base-rotation NDI class.
4. measurementValid: boolean orientationValid_control.
5. baseValid: boolean baseOrientationValid.
6. measurementTime_s: scalar double conditional successful-poll simulation time.
7. t: scalar double simulation time.
8. enablePressure: scalar boolean. Retain the independent downstream Kill.
9. lambda: double 2x1 Constant `[200;200]`, 1/s.
10. dampingGain: double 2x1 Constant `[20;20]`, N s/m (controller gain Ks).
11. enableAdaptation: scalar boolean Constant. False freezes estimates;
    true permits supervised learning. It does not reset estimates/history.

Keep gains constant within a run. Add dashboard controls only if necessary;
runtime changes to lambda are not covered by the constant-gain derivation.
Fresh duplicate timestamps cannot update estimates twice. Gaps over 0.075 s
skip history integration/learning. Restart the simulation to reset states.
Disabling pressure does not erase learned estimates; fault remains latched.

## Outputs

The first 17 outputs match the nominal wrapper in order. Connect output 1
directly to the existing pressure/Kill path, bypassing the old inverse block.
Log all outputs and the input gain/enable signals.

Additional outputs:

| Port | Name | Size / type |
| --- | --- | --- |
| 18 | stiffnessEstimate_NPerM | 2x2 double |
| 19 | dampingEstimate_NsPerM | 2x2 double |
| 20 | hysteresisForce_N | 2x1 double |
| 21 | adaptiveBias_N | 2x1 double |
| 22 | parameterEstimate | 10x1 double |
| 23 | historyState | 3x1 double |
| 24 | adaptationActive | scalar boolean |
| 25 | parameterAtBound | scalar boolean |
| 26 | adaptiveForce_N | 2x1 double |

Set port 13 status explicitly to uint8. Set ports 11,12,14,15,24,25 to boolean.
Other outputs are double. Parameter estimates/history correspond to the current
force evaluation; parameter updates take effect on the next sample.

## Adaptive model

Let T=[-1 -1;1 0;0 1], so muscle length changes are T*q. Use

    Khat = [k1+k3, k3; k3, k2+k3]
    Dhat = diag([d2,d3])
    Hhat = T' * diag(z) * [h1;h2;h3]
    bhat = [b2;b3]

Khat has positive independent springs and nonnegative shared coupling.
This restricted family preserves positive definiteness; it is not an arbitrary
symmetric stiffness matrix. Dhat remains positive. Initial K is exactly the
nominal 650*[2 1;1 2], and initial D is 40*I. History/bias weights start at zero.
Smooth length-bound forces remain fixed and outside adaptation.

    v = dqd - lambda.*(q-qd)
    a = ddqd - lambda.*(dq-dqd)
    s = dq-v
    tau = M*a + C*v + G + Khat*q + Dhat*v + Hhat + bhat - Ks.*s

IMPORTANT: Dhat*v replaces the nominal wrapper's D*dq. Even with adaptation
disabled, this changes controller behavior. Near a static hold its direct
position coefficient includes (Dhat+Ks)*lambda, not just Ks*lambda. At the
initial settings that is 12000 N/m. Do not interpret enableAdaptation=false
as reproducing the old nominal wrapper exactly. Compare adaptive on/off using
THIS same new controller to isolate the learning effect.

Parameter order is k1,k2,k3,d2,d3,h1,h2,h3,b2,b3. The 2x10 regressor uses
spring columns [q2;0], [0;q3], [q2+q3;q2+q3], damping columns [v2;0],
[0;v3], history columns T'*diag(z), and bias columns I.

Update normalized theta/scale using -Gamma*Yscaled'*s; apply per-parameter rate
caps and box projection. Settings/bounds are in pacAdaptiveConfig.m.
Updates occur only in status3, with valid/fresh inputs and no pressure/correction/
slew limiting, residual norm <=0.02 N, and norm(s)>0.005 m/s. There is no
leakage by default. Bias can integrate persistent mismatch; damping and bias
are correlated at static holds, so physical parameter identification is not
guaranteed. Static data alone cannot identify damping uniquely.

## History scope and stability limits

Three normalized Bouc-Wen n=1 feature states use measured muscle velocities:
z_dot=w-beta*abs(w).*z-gamma*w.*abs(z), w=T*dq/0.005 m,
beta=1, gamma=0.5. Bounded substeps limit numerical excursions. Force weights
adapt; alpha/beta/gamma of the internal state law DO NOT adapt. Those require
a separate nonlinear-state sensitivity/observer derivation. This is an adaptive
history-force representation, not online identification of all Bouc-Wen parameters.

Mass, C and gravity remain nominal. The paper's stiffness/damping adaptation
does not adapt inertia. Our bias/history extension is not a direct reproduction
of its theorem. Ideal, unprojected continuous-time regressor cancellation is
checked numerically; rate caps, deadband, history mismatch, sampling, filtering
and pressure constraints prevent claiming exact asymptotic/passivity guarantees
for hardware. Projections bound estimates, not the entire physical closed loop.

The existing trim forceBias stays captured once per enable cycle. It is distinct
from adaptiveBias. Hysteresis features and bias must NOT be re-zeroed at every
pose change; that would erase the history being represented.

## Offline checks and run analysis

    validatePacAdaptive
    validatePacAdaptiveTracking
    replayPacAdaptive('PAC_circle/static_phi_lambda_200_damp_20_sweep.mat')

The tracking check uses ideal continuous force and a chosen representable
mismatch. Replay holds the recorded plant response fixed and is counterfactual
numerical validation, not evidence of future physical tracking performance.
Replay may hit pressure rails as estimates grow against an unresponsive recorded
trajectory; learning freezes on those samples.

First physical comparison: existing 0.2 rad static ramp/hold, fixed gain pair;
test new wrapper with learning disabled, then enabled. Check engagement,
force transients, parameter bounds and limits before running the full sweep.
No new identification experiment is required before implementing this path.

    analyzePacAdaptive('your_recording.mat')
    analyzePacStaticSweep('your_recording.mat') % add PAC_circle to path

Keep parameterEstimate and adaptive diagnostic names exactly as listed.
The analyzer prints the parameter ranges, force contributions, active/bound/
limited/fault counts and plots estimate histories. Static analyzer remains
responsible for length/configuration accuracy. No script runs NI/NDI hardware.

parameterAtBound includes zero history weights at their lower bound initially;
this flag alone is not a fault. Inspect which coefficients remain at bounds.
The adaptive report/plots were verified with synthetic wrapper signal logs.

## 2026-10-13 - Activation-consistent fixed-mean PAC allocator
User authorized pacAllocatePressure change. Exact activation-map inverse now
uses requested+forceBias to derive effective pressure differences. Preserves1.6
bar mean; when effective muscle pressure iszero, its inactive commandedpressure
can decrease below.8 to satisfyfixedmean. Candidate accepted only withinexisting
.2-3bar rails ands/a +/-2.5trimcorrectionclamps. Infeasiblecandidate retains bounded
legacyfallback,nowlimited=true fornonzeroallocationresidual; no hiddenmeanincrease
or relaxedconstraints. Supervisor2bar/slew remains andmayproduce transientresidual.
Nominal activationthreshold.8 remains assumption; doesnot identifyphysicalpressure.
validatePacPressure includes feasibleinactiveexactforceand infeasibleflag checks;
validatePacAdaptive passes. No SLX edited byagent. Repeat1rad30s baseline,
lambda40/Ks20 andrecordallocationresidual/limits/updates beforehighergainorspeed.
