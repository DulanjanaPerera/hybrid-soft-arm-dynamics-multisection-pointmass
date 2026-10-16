# PAC recordings and offline static analysis

Run from the project root:

```matlab
addpath('PAC_circle');
results = analyzePacStatic('PAC_circle/static_phi_test1.mat');
```

Paste the printed report into the discussion. The function saves a text report
and a compact analysis MAT beside the source file, and opens tracking/force
plots. It only reads recordings and never loads or executes a Simulink model.

Optional explicit constant-reference window and plot switch:

```matlab
results = analyzePacStatic('PAC_circle/static_phi_test1.mat',[40 60],false);
```

Automatic selection uses the longest nonzero constant theta/phi plateau,
and reports settled metrics over its latter half. For multiple holds, call
with each desired window. Gain estimates derive from logged -Ks*s and
s=dq-dqd+lambda*(q-qd), assuming desired rates have settled to zero; they
are not a replacement for explicitly logging gains. Missing optional signals
produce NaN; conflicting duplicate names fail rather than select arbitrarily.
Theta error is reported only where measured phi exceeds 0.05 rad.

## Latest static baseline (2026-10-11)

Reported lambda=200, Ks=20 trials at theta=0, pi/2, pi show
direction-dependent offsets without hold saturation. pi/2 has good bend
magnitude but a 0.16734 rad direction MAE. Keep both length and angle errors
in comparisons. See project_docs/VALIDATION_RESULTS.md for the compact table;
adaptation is pending and mismatch has not been uniquely identified.

## Combined static sweep

Paste pacStaticSweepReferenceBlock.m into the reference block. Ten trials:
theta=0,pi/2,pi each at phi=.1,.2,.3, then repeat theta=0/phi=.2.
Each trial: 5 s rise, 30 s hold, 5 s fall, 10 s zero-reference rest. Starts
at 25 s; last trial ends 525 s. Set Stop Time 535 s and check status=3 before
25 s. Zero phi means zero desired bend, not vented pressure. Controller remains
enabled at common pressure. Theta changes at zero desired phi; residual actual
bend may persist. Keep gains fixed for the whole run. Log testId/phaseId as
well as the existing reference, measurement and PAC diagnostic signals.

```matlab
addpath('PAC_circle');
r = analyzePacStaticSweep('PAC_circle/static_sweep.mat');
```

Prints a compact report and summary; saves TXT, MAT and CSV beside recording.
Uses logged constant-reference plateaus, not fixed times. Detailed per-hold
reports remain in results.holds. Exploratory feedback=deltaK*q+b fit describes
aggregate static mismatch, not uniquely physical stiffness. Fit requires at
least four usable holds with independent length-coordinate variation.

Workspace use: analyzePacStaticSweep() or analyzePacStaticSweep(out). Workspace calls do not save artifacts; MAT-path input retains report/analysis/CSV saving.
