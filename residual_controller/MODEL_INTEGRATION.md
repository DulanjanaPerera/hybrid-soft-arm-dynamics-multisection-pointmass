# Residual model integration and PAC parameter interpretation

2026-11-03. Algebraic implementation contract, not a new stability theorem.
q=[l2;l3], l1=-l2-l3. The frozen force GRU was trained with

    L0 = M*ddq+C*dq+D0*dq+K0*q+bound+G0 = B0*p+r0(q,dq,p,h).

K0=650*[2 1;1 2],D0=40I,G0 uses[0;0;-9.81]. The bound is the same
all-three-actuator force in forward/inverse/PAC. No identified-force labels
exist: r0 is trained through measured next-length/path prediction errors.

## Re-expression instead of silently changing the nominal definition

The controller's activation map is Fcfg(p), and its gravity term is Gg.
Define

    r_eff = r0 + B0*p-Fcfg(p) + Gg-G0.
    Lg = L0+Gg-G0 = Fcfg(p)+r_eff = B0*p+r0+Gg-G0.

This preserves the exact trained model when code uses a different nominal
activation/gravity decomposition. It does NOT teach the network a new gravity,
base orientation, pressure actuator coefficient or payload. Original operating
coverage still applies. Runtime requires the same force-per-bar coefficient;
the active model-on contract rejects a different coefficient.

All-active equal activation thresholds already cancel in force differences.
Below that branch the explicit activation lift avoids applying an additional
deadzone to a residual trained against B0*p. The correct pressure derivative
at fixed h is B0+dr0/dp, which is used by the bounded inverse. No true actuator
force sensor or physical model-error bound is introduced by this identity.

## PDFB and pressure calibration

Keep the existing virtual acceleration, including outer PID/integral terms:

    tauPDFB = M*a_virtual+(C+D0)*dq+K0*q+bound+Gg.

At a valid trim, define b=Fcfg(pTrim)+r_eff(pTrim,hAnchor)-staticAnchor.
The common pressure supervisor requests a blended model force tauBlend. Solve

    Fcfg(p)+r_eff(p,h_previous)-b = tauBlend

using bounded commands and the existing common-pressure policy. The anchor
and inverse use the same definition; otherwise a change in constant residual
could be credited incorrectly to trim. Bias is held after capture; it does
not get re-fit at every candidate. The state never advances during the solve.

## PAC keeps one nominal baseline, with adaptive corrections

Let theta0=[650;650;650;40;40;0;0;0;0;0], delta=theta-theta0,
v=qd_dot-lambda.*(q-qd), a=qd_ddot-lambda.*(dq-qd_dot), s=dq-v.
The existing force is EXACTLY

    tauPAC = M*a+C*v+D0*v+K0*q+bound+Gg-Ks.*s+Y*delta.

The regressor uses v for damping, q for stiffness, existing history states
for hysteresis and identity for bias. The source does not add a second K0/D0
term. It applies this force demand to the augmented inverse above. theta still
appears absolute in the old ports to preserve interfaces; delta is the
remaining learned-model correction, not an identified physical stiffness,
damping or hysteresis parameter. Fresh theta0 gives zero initial correction.
Existing gamma,scale,projection,rate caps and hysteresis-freeze behavior remain.

For an exact unconstrained augmented force delivery, define remaining physical/
calibration mismatch epsilon by Lg=Fcfg+r_eff-b+epsilon. The error equation is

    M*ds/dt + (C+D0+diag(Ks))*s = Y*delta+epsilon.

This explains why the existing gradient direction can continue learning the
remaining error. It does NOT establish epsilon is bounded, matched, constant-
parameter representable, or sufficiently small. Frozen RNN history/pressure
dependence, imperfect identification, sampled updates, pressure constraints,
rate caps, sensor dropouts and signed stiffness projections preclude importing
the old ideal passivity proof without a separate analysis. Adaptation gates use
augmented numerical allocation residual; that quantity is not epsilon or a
measured force error. Physical stability/performance needs separate evaluation.

The original model-off PAC path is delegated without semantic changes. PDFB
uses a new paired nominal/residual wrapper because the saved embedded inverse
retained a legacy bound/pressure path; both new conditions share the corrected
nominal definition and allocator. Record a fresh new-wrapper nominal baseline.
