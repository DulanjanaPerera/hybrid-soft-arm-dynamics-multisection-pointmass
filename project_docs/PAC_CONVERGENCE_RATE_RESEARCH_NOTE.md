# PAC gain relationships and convergence-rate research proposal

Date:2026-10-28. Status: literature comparison and conditional analytic results,
not an implemented controller, proven new robot theorem or novelty claim.
User proposes gain ranges for a desired tracking/adaptation convergence time.

## Verified source and mapping

Read/rendered pages1900-1901 of the supplied
Dynamic_Control_of_Soft_Robotic_Arm_An_Experimental_Study.pdf (Azizkhani et al.,
IEEE RA-L8(4),2023). Eq11 uses r=edot+Lambda e and feedback -KG r. TableI has
Lambda_AP=8,kg=8,kp=50^2,kd=100,GammaK=1e-5,GammaD=1e-2. Eq19 is
thetaHat_dot=-Gamma^{-1}Yp^T r: small tabulated Gamma means a large multiplier.
Their stiffness/damping update multipliers are therefore1e5/100, before unit,
coordinate and basis differences. The equality Lambda_AP=kg=8 is not a required
algebraic relationship. kp/kd belong to the separate PDFL benchmark in Eq22.
Ideal exact feedback linearization yields eddot+100edot+2500e=0, poles -50/-50.
Their setup executes control/encoder measurements at1kHz; optical camera16Hz.
Our PAC is intended at20Hz and uses inferred two-coordinate lengths rather than
their three encoder elongations. Actuation/sensing and regressor differences matter.

Current pacAdaptiveUpdate instead uses
thetaHat_dot=-diag(scale.^2.*adaptationGamma)*Y^T*s before clipping/gates.
Unloaded baseline multipliers are5e5 for stiffness and20 for damping (not simply
.5/.002); hysteresis/bias have separate regressors absent from their six-parameter
law. Numerical multipliers alone do not compare closed-loop performance.

## What can be assigned like poles

For e=q-qd,s=edot+Lambda e, exact e dynamics are edot=-Lambda e+s.
Only on s=0 is the decay exactly exp(-Lambda*t). Baseline Lambda18/16.5 gives
filter-only2% times ln(50)/lambda=.2173/.2371s; this is not a prediction that the
actual adaptive arm settles that fast. Persistent s drives ongoing error.

For an ideal matched plant, exact force delivery and continuous-time cancellation,
M*sdot+(C+D+Ks)*s=0. In a scalar frozen-coordinate approximation with C=0,
m>0 and actual d>=0, alpha_s=(d+ks)/m and the cascade poles are -lambda and
-alpha_s. Here d is actual damping, not the estimated feedforward coefficient.
These are local nominal design relations, not global poles of the real adaptive,
sampled, gated, pressure-limited system. Repeated poles produce t*exp(-lambda*t),
so multiplying an inverse pole by4 is not a rigorous complete settling-time bound.

## A conditional practical tracking-bound proposition (own derivation)

Assume scalar constant m,lambda,d,ks, exact s=edot+lambda*e, and an aggregate
residual w with verified |w(t)|<=Delta for all considered time:

  m*sdot+(d+ks)*s=w.

Let a=(d+ks)/m, S=Delta/(d+ks) and
F(t)=(exp(-a*t)-exp(-lambda*t))/(lambda-a), or t*exp(-lambda*t) when a=lambda.
Variation of constants and the triangle inequality give

  |s(t)| <= |s0|exp(-a*t)+S*(1-exp(-a*t)),
  |e(t)| <= |e0|exp(-lambda*t)+|s0|F(t)
            +S*((1-exp(-lambda*t))/lambda-F(t)).

Hence limsup|e| <= Delta/(lambda*(d+ks)). Choose a target epsilon strictly above
this floor, and bound initial e0/s0 on an operating set, before claiming an entry
time from the envelope. This relation explains how lambda,feedback damping and
residual mismatch jointly matter. It does not provide Delta or a Gamma-dependent
adaptation rate. With actuator lag, nonconstant inertia, filters, sampled data or
adaptation projection, its assumptions need extension/reproof. Nominal allocation
residual near0N does not measure/bound the physical w.

## Why current adaptation is not arbitrary pole placement

The ideal matched gradient proof uses V=.5*s^T*M*s+.5*thetaTilde^T*G^{-1}*thetaTilde,
G=diag(scale.^2.*gamma). Cancellation yields Vdot=-s^T*(D+Ks)*s, without a
negative parameter-error quadratic term. This alone does not imply Vdot<=-alpha V
or an assigned augmented-state decay rate. Tracking convergence and physical
parameter convergence are different claims. Excitation/identifiability are needed
for parameter-rate results; a fixed circle can leave correlated parameters.

The actual update additionally has a deadband, projection, component rate limits,
measurement/limiting gates and Euler updates. An adequate adaptation-availability
condition is needed for a uniform wall-clock learning deadline. If learning can
freeze indefinitely or required estimates lie outside the allowed set, no positive
Gamma alone guarantees the desired parameter deadline. Component rate clipping
is not automatically a Lyapunov-compatible projection. History-force states in our
code are Bouc-Wen states, not a composite-learning data-memory stack.

For bounded additive mismatch, linear -Ks*s produces a term
Vdot<=-k||s||^2+Delta||s||, not automatic cancellation for every nonzero s.
The source paper's Eq16->20 uncertainty domination should be rederived with
explicit assumptions rather than inherited as a zero-error guarantee. A practical
error tube is the appropriate target unless residual vanishing is justified.

## Candidate theorem target and proof work

Aim for practical time-to-tube on a compact operating set, with an identifiable
parameter basis, verified excitation/information lower bound, bounded actuator/
model/sensor error, feasible pressure/slew commands and bounded update outages.
A modified composite/prediction-error adaptation law may be required.

If a valid augmented Lyapunov function W including e satisfies
W>=c_e*||e||^2 and Wdot<=-alpha*W+beta, then

  W(t)<=W_inf+(W0-W_inf)*exp(-alpha*t), W_inf=beta/alpha.

When W0>c_e*epsilon^2>W_inf, a sufficient guaranteed entry-and-stay time is

  T_epsilon <= (1/alpha)*ln((W0-W_inf)/(c_e*epsilon^2-W_inf)).

This is a target proof structure, not a result established for our controller.
Specify all constants, the region of validity, information activation time and
feasibility checks. If excitation begins at t_e, account for t_e and W(t_e).
Bounded data/gate outages may require hybrid or active-time analysis. If a hard
exact convergence deadline is desired, distinguish finite/fixed/prescribed-time
designs; do not relabel ordinary exponential convergence as finite-time stability.
Pointwise e bounds imply suitable RMSE bounds, but passing last3-circle RMSE does
not establish a pointwise settling-time theorem.

## Prior art checked (primary sources)

- Slotine & Li, Composite adaptive control of robot manipulators, Automatica1989:
  https://www.sciencedirect.com/science/article/pii/0005109889900940
- Guo,Pan,Zheng,Yu, Composite learning control of robotic systems: A least squares
  modulated approach, Automatica111(2020)108612, DOI10.1016/j.automatica.2019.108612:
  https://www.sciencedirect.com/science/article/abs/pii/S000510981930473X
  Explicitly addresses balanced,easily tunable parameter convergence under IE.
- Dong,Hu,Akella, Composite Adaptive Control for Robot Manipulator Systems,
  EuroGNC2019, finite-excitation parameter convergence:
  https://eurognc.ceas.org/archive/EuroGNC2019/html/CEAS-GNC-2019-062.html
- Zhang et al., Neural networks-based composite learning control for robotic
  systems with predefined time error constraints, Neurocomputing608(2024)128414,
  DOI10.1016/j.neucom.2024.128414:
  https://www.sciencedirect.com/science/article/pii/S0925231224011858
  Prior work on preassigned error-region/time and IE-based learning.

The general idea of assigning adaptive convergence rates is established prior art.
Potential contribution is a verified sampled pneumatic-arm result combining
actuator/geometry uncertainty, pressure/slew feasibility, intermittent adaptation
and a measurable length-error time bound. Novelty and the full proof remain open.
Next analytic step: derive the actual augmented error system and check normalized,
gate-weighted regressor information/rank before proposing a gain-range theorem.
No controller/campaign changes, new hardware trial or gain recommendation here.
