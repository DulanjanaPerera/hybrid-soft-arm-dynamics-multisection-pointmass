# Reusable technical notes

Current model overview is in MODEL_ARCHITECTURE; evidence in VALIDATION_RESULTS.
Historical derivations/discussions retain their original context in the log.

## Reduced coordinates, virtual work and common-pressure nullspace

For fixed backbone and symmetric geometry, `l1+l2+l3=0`.
Substitute `l1=-l2-l3` into actuator work:

`F1*dl1+F2*dl2+F3*dl3=(F2-F1)*dl2+(F3-F1)*dl3`.

Thus `tau=[F2-F1;F3-F1]`, and for equal nominal pressure-force gains,
`tau=f*[-1 1 0;-1 0 1]*p`, `f=pi*(0.0065)^2*1e5=13.2732 N/bar`.
The nullspace is span([1;1;1]). This algebra holds in the ideal bending model;
it does not prove common pressure leaves physical response unchanged.
Backbone load, PMA bulging, prestress and history may alter incremental bending.
Z displacement can arise from bending at fixed arc length and alone does not
identify axial shortening.

## Fixed-mean pressure coordinates and feasibility

Define mean c, symmetric differential s and antisymmetric differential a:

`p=[c-2*s/3; c+s/3-a; c+s/3+a]`.

Then `p2-p1=s-a`, `p3-p1=s+a`, so for incremental generalized force:
`ds=(dtau1+dtau2)/(2*f)`, `da=(dtau2-dtau1)/(2*f)`.
All retain mean c until clipping or asynchronous channel slew changes it.
The manual settings are absolute s/a, not accumulated changes; a trial ID is
an annotation, not a pressure setting or proof that a hold settled.

More generally let d2=tau1/f, d3=tau2/f. Every ideal solution is
`[p1;p1+d2;p1+d3]`. A fixed mean c gives `p1=c-(d2+d3)/3`.
For chosen active floor pmin and ceiling pmax, feasible p1 lies in
`[max(pmin,pmin-d2,pmin-d3), min(pmax,pmax-d2,pmax-d3)]`.
Check the interval; independently clipping channels changes delivered force.
A flexible common level offers more authority than a fixed mean.

## Deadzone versus mechanical hysteresis

The older allocator maps each positive effective pressure to command
`p_eff+0.8`; inactive channels command zero. Consequently arbitrarily small
activation produces a command jump near 0.8 bar. Mechanical Bouc-Wen force
cannot remove this allocation discontinuity. In an all-active regime with
identical offsets, the offsets cancel in pressure differences; do not add them
a second time around the empirical trim. Real channel offsets/gains may differ.
Initial/final ramps cross activation even when operating plateaus stay active.

Older candidate mechanical-force hysteresis used
`hdot=alpha_h*dq-beta_h*abs(dq).*h-gamma_h*dq.*abs(h)`, with
alpha_h=65.59, beta_h=49.04, gamma_h=-19.93, K approximately
[[2576.3,1271.0],[1271.0,1886.0]], D=diag([16.15,18.50]).
These came from the external older P2/P3-chirp simulation fit documented in
the October 5 log. Its [0,l2,l3]/direct F2,F3 convention differs when P1 is active.
They are candidate estimates, not current parameters. Adaptive coupled K must
remain positive definite; adapting K/D does not itself identify valve memory.

## Smooth bounds and apparent stiffness

The implementation augments each diagonal K entry by
`0.5*b*(2+tanh(mu*(q_i-upper))-tanh(mu*(q_i-lower)))`,
then multiplies Keff by q. Current common-pressure values are mu=2000,
b=1e6 and lower/upper=-0.035/+0.035 m. This is a smooth nonlinear restoring
penalty, not a hard state clip or measured physical limit. With symmetric
l1=l2=-l3/2, a -30 mm l3 plateau implies +15 mm l1/l2; no separate 15 mm cap.

A command-based static force-residual fit combines mechanical stiffness with
assumed pressure gain, deadzone, gravity, calibration and hysteresis. Changing
k alone cannot independently separate them. Pooled loading/unloading P1 means
and a phi-only PCHIP do not establish a history-independent constitutive law.
Standalone out-of-range NaN and block endpoint extension are different policies.

## Incremental anchoring and PD steady error

At trim engagement capture only static `tauAnchor=G+Keff*q`.
Allocate `deltaTau=tauDesired-tauAnchor` around empirical pressure trim.
Absolute tauDesired and ideal tauAvailable of the rate-limited command need
not match at trim; this is not automatically a force-reconstruction bug.
Anchoring the entire desired force would also cancel the initial PD correction.

PD correction vanishes in velocity at rest; a constant unmodeled bias generally
leaves position error unless gains/authority or a bias-compensation mechanism
remove it. Higher gain may encounter delay, slew, oscillation or saturation.
The latest residual was inside local pressure limits; the data do not identify
a stable new gain or integral design. Invalid-input returns currently bypass
normal slew; nominal smooth-output checks do not validate fault transitions.

## Measurement, timing and equilibrium interpretation

NDI orientation-derived lengths and XYZ are distinct measurements. Desired
task XYZ for circle plotting includes base/arc/tip geometry [0.055,0.17411,0.05] m.
Unwrap measured theta for trajectory comparison and use shortest signed angular
error; theta is ill-conditioned near straight. Signed U/V preserve bending
direction where nonnegative phi cannot.

pollStart is a practical acquisition-time proxy, not calibrated device frame
time; hostTime includes reply latency. Simulation-time replay and physical-time
delay analysis are different experiments. Sparse parameter log events do not
imply changed gains; hold their last values when aligning data. Mask invalid
orientation/pair/configuration samples and do not treat held displays as new data.
NDI startup preflight validates the saved mounting reference, not a bent pose
correction. Only one process can own COM12; ordinary hardware-model simulation
is not automatically hardware-free.

Final-window means are not proof of equilibrium. Compare slopes, ranges,
loading histories and repeated returns. Inverse analysis uses >=5 s physical
holds, last 3 s, speed <=0.20 mm/s and coordinate range <=0.50 mm for stationarity;
settling requires <=1 mm errors through the end for >=2 s. These are engineering
criteria, not identified sensor limits. A zero settling time can mean the target
was already satisfied. Local READMEs retain calibration/preflight and analysis
commands, thresholds and detailed output definitions.

## Earlier tapped-delay predictor and recursive NMPC

Audited old architecture: three pressures, backward pressure differences and
tanh(dP/0.05) loading indicators at delays 0..8, plus three previous lengths
for each scalar output: 84 features. Actual main_NARX uses separate l1/l2/l3
`fitnet([64 32],'trainbr')`, output history enabled and corresponding active-only
20 s chirp trials. Earlier [32 16] handoff architecture was an example.

Validation and offline predictor update history from recorded lengths, so they
demonstrate measurement-assisted predictions, not autonomous rollouts. The
equilibrium naming does not establish a physical equilibrium output.
For future NMPC, initialize real history, derive future pressure-rate/loading
features from candidate commands, insert predicted future lengths into output
taps, and keep real history unchanged during optimization. Future references
may be known trajectories or a held target; future measured lengths are unnecessary.
Use two independent outputs or enforce the length-sum constraint.
Audit alignment/sample period (default Ts=1 s is unsafe to assume), centered
Savitzky–Golay labels and normalization before block splitting for leakage.
Separate pressure/valve memory from passive mechanical hysteresis to avoid double
counting. A successful current one-step plot is insufficient for control rollout.
