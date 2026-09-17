# Continuum-arm dynamics progress log

Last updated: 2026-09-16. This file records the state observed in the local
MATLAB repository after the work described below. Read it together with
`CONTINUUM_ARM_DYNAMICS_HANDOFF.md` and inspect the current code before editing:
this log is a checkpoint, not a substitute for the working tree.

## Current repository state

- Branch: `compare-pointmass-standard`; observed HEAD: `07062c6` (`initial handoff`).
- The beta work and several validation artifacts are local changes. They have
  **not** been pushed or merged to `Ozi` or another branch.
- The working tree also contains the three reference PDFs and earlier analysis
  files. Preserve them; do not discard or overwrite unrelated user work.
- MATLAB is R2025a on Windows. MEX files use `.mexw64`. Generated build output
  lives under `C:\MATLAB_build`, outside OneDrive; MATLAB source and copied
  binaries in this directory are the working inputs.
- The point-mass MEX was rebuilt after the final source changes and was observed
  at `armS_dynamics_N3_entry_mex_mex.mexw64` on 2026-09-16 19:20 local time.
  The standard MEX was observed at `armS_standard_mex.mexw64` at 18:33.

## Modeling choices that must remain explicit

- Three sections, two independent length changes per section:
  `q = [l12;l13;l22;l23;l32;l33]`, `X = [q;dq]` (12 states). Local section
  input is `[0,q(2*n-1),q(2*n)]`.
- Translational kinetic energy only. Upstream rotations still affect point
  translation. Rotational kinetic energy is excluded.
- The standard model distributes each fixed section mass uniformly over
  normalized material position `xi` in `[0,1]`; it uses analytical integrals.
  The point-mass model places each mass at the selected backbone position
  `cog_xi(n)`, which need not be the integrated center of gravity.
- Preserve the implemented gravity convention, six-coordinate kinematics,
  13th-order Taylor HTM, stiffness law, damping, actuation, and input law.
  The distributed reference has no beta scaling.
- The saved **comparison** has `cog_xi = [0.5;0.5;0.5]`. The two older
  `runDynamicSimulation_armS_nume*` scripts currently initialize `loc = 0.1`.
  Do not confuse their current defaults with the saved comparison setup.
- A three-coefficient shared-beta candidate has now been fitted offline
  against the distributed model. It has **not** replaced the production
  all-ones default or been treated as experimentally calibrated.

## Standard distributed-mass model: completed

`armS_standard_core.m` uses the original `integratedPosition_nume` for `mu`
and the compact, algebraically equivalent S/F/E functions. Original S/F/E
functions remain for equivalence checks. The compact S/E comparison passed
70 cases at two geometries, with maximum relative errors below `5e-15`;
individual S and E source generation succeeded. Compact F had already passed
its corresponding comparison and isolated source-generation probe.

After switching all three calls, the assembled standard model passed its
three-pose symmetry, positive-definiteness, mass-derivative, Christoffel-skew,
and direct 16-node quadrature checks. The full standard MEX built in 178.6 s.
Its three MATLAB/MEX RHS checks had worst scaled error `5.731e-15`. A compiled
trajectory compared with the saved interpreted standard trajectory differed
by at most `5.733e-14 m` in coordinates and `9.378e-13 m/s` in velocities.
Rebuild if the standard source changes; a successful earlier binary does not
validate subsequent edits.

## Point-mass beta implementation: current

The CoG paper's translational kinetic-energy terms (equations (22) and (25))
introduce three coefficients. The user-provided block form, implemented per
section, is below. `A` is the upstream translational velocity Jacobian,
`B` is the upstream angular-Jacobian blocks applied to local mass position
`p`, and `P = p_q` is the local position Jacobian.

```text
M11 = m*(A'*A + A'*B + B'*A + beta_v1*B'*B)
M12 = m*(A'*P + beta_v2*B'*P)
M21 = M12'
M22 = m*beta_v3*(P'*P)
```

Only `M22` remains for section 1. `beta` is a `3 x 3` numeric matrix: rows
are sections 1 to 3, columns are `[beta_v1,beta_v2,beta_v3]`. It is constant
with respect to `q` during a simulation, so its derivative is zero, but its
multiplier remains in every product-rule derivative. See `POINTMASS_BETA.md`
for the runtime interface and example.

The inconsistency found in the previous code affected **both** assembly and
derivatives: `armS_core_N3_mex.m` assembled unit-coefficient blocks (and
symmetrized the assembled matrix), while `Mi_h.m` hard-coded beta values
independently. `Mi_h.m` now differentiates the same beta-dependent blocks with
complete product rules. Section mass multiplies both `M` and every derivative
slice. No post-hoc symmetrization is used to conceal derivative errors.
The Christoffel convention remains
`dM(i,j,h)=partial M(i,j)/partial q(h)` and
`C(j,k)=0.5*sum_h[(dM(j,k,h)+dM(j,h,k)-dM(k,h,j))*dq(h)]`.

Active compiled chain:

```text
runDynamicSimulation_armS_nume*.m
  -> armS_dynamics_N3_entry_mex_mex(t,X,L,r,cog_xi,mi,g,K,D,tau,mu,lKbounds,beta)
  -> armS_dynamics_N3_entry_mex.m
  -> armS_core_N3_mex.m
  -> Mi_h.m / christoffelSymbol.m
```

`build_armS_N3_mex.m` builds that 13-input interface, with `beta` as a runtime
input. Existing direct calls with 12 inputs to the **new compiled MEX** must
be updated. The interpreted `armS_dynamics_nume.m` path also accepts
`params.beta`, defaulting to `ones(N,3)` if absent. The interpreted entry and
core default to ones if beta is omitted. The older
`armS_dynamics_nume_N3.m`, `armS_dynamics_nume_v2.m`, and
`armS_dynamics_recursive.m` are alternate/legacy paths, not called by the
active scripts or MEX entry; they are not tunable-beta interfaces.

## Beta validation completed

`validate_pointmass_beta.m` tests three poses and two beta matrices:
`ones(3,3)` and
`[1.2 1.05 1.3; 1.1 0.9 1.2; 1.3 1.1 1.4]`. The latter is a test input,
**not a calibrated value**.

- Every finite-difference `dM(:,:,h)` check passed; worst relative error
  `6.053e-11`. The finite-difference displacement was `1e-7 m`.
- Mass-matrix symmetry was at roundoff; the largest reported
  `Mdot-2*C` skew residual was `2.327e-16`.
- Cholesky passed at all three tested poses for both matrices. Unit-beta
  `rcond(M)` was only `2.509e-7` at the straight pose, so conditioning needs
  attention even when Cholesky succeeds.
- Positive coefficients do **not** guarantee a positive-definite `M`:
  `[0.1,5,0.1]` repeated in all three rows gave Cholesky flag 4 and minimum
  eigenvalue `-16.90` at the reference bend.
- The rebuilt MEX and interpreted MATLAB paths agreed at both beta matrices
  to worst scaled RHS error `4.379e-13`. Changing beta in a call to the same
  binary changed the computed acceleration; no rebuild was needed.
- At all-one beta, the rebuilt MEX reproduced the saved corrected point-mass
  trajectory with maximum coordinate difference `2.063e-10 m` and velocity
  difference `5.072e-9 m/s`. Before rebuilding, the changed MATLAB RHS had
  also matched the previous unit-beta binary to `7.494e-13` scaled error.

Saved checks: `pointmass_beta_validation.mat`. No parameter optimization or
experimental identification has been done.

## Standard versus point mass at beta = 1

`comparePointmassStandardExperiment.m` used the saved corrected point-mass
parameters and 301 times from 0 to 5 s. Both compiled models used the same
initial state, solver (`ode15s`, `RelTol=1e-8`, `AbsTol=1e-10`,
`MaxStep=1e-3`), masses, geometry, damping, stiffness, and inputs. The
point-mass comparison location was `xi=0.5` in each section. The point-mass
free run reproduced its saved reference exactly in that experiment.

| Start and input | Maximum coordinate difference | Maximum tip-position difference |
|---|---:|---:|
| Reference bend, free | 2.520 mm | 261.3 mm |
| Asymmetric bend, free | 1.337 mm | 144.0 mm |
| Reference bend, differential | 2.527 mm | 261.3 mm |
| Asymmetric bend, differential | 1.346 mm | 144.4 mm |

The differential input added `[0.5;-0.5;0;0;0;0]` to the saved zero input.
The bent-case gaps changed little under that input. Coordinate difference
means the largest absolute difference across the six `q` values and sample
times; tip difference is Euclidean distance between kinematic tip positions.
These are **model disagreements**, not errors against physical measurements.
The straight free case was far closer. Results and plots are in
`comparison_pointmass_standard_experiment.mat`,
`comparison_pointmass_standard_errors.png`, and
`comparison_four_bent_cases.png`.

`diagnose_midpoint_beta_scaling.m` fitted one multiplier for each *entire*
section-mass contribution at the two bent initial poses. This deliberately
limited test is **not** a fit of the paper's three coefficients. Its best
scales `[29.4265,0.1339,1.1562]` only partly reduced mass-matrix gaps and
greatly worsened gravity agreement. Do not treat those numbers as candidate
translational betas or mass estimates.

## Further numerical checks already run

- `validate_strong_bend_standard.m` compared the standard core with direct
  quadrature at poses through 40 mm coordinate magnitude. At 20, 30, and
  40 mm symmetric bends, local Taylor-rotation orthogonality defects were
  `2.923e-6`, `4.491e-4`, and `1.698e-2`; relative mass-quadrature gaps were
  `2.781e-8`, `7.911e-6`, and `1.525e-3`. Derivative and Christoffel checks
  stayed numerically consistent. The 30 and 40 mm probes extend beyond the
  current ±20 mm length-limit settings. Results: `strong_bend_validation.mat`.
- `check_arm_energy.m` evaluated zero-input energy on the saved 301-sample
  trajectories with the implemented gravity convention and integrated
  nonlinear elastic potential. Total energy decreased at every recorded
  sample. Standard and point-mass energy drops were `1.16893 J` and
  `1.13934 J`; maximum residuals of energy drop versus integrated damping
  were `0.0904%` and `0.185%` of those drops. Output-time quadrature limits
  this check. Finite differences verified the potential gradients. Results:
  `arm_energy_check.mat`.

## Shared-beta tuning checkpoint

`POINTMASS_BETA_TUNING_WORKFLOW.md` contains the reproducible objective,
pose split, feasible cone, results, and timing method. The candidate shared
across sections is `[1.46573213,1.96512396,3.00000000]`; its third value
hits the current search bound. Normalized mass-matrix RMS fell from `0.32941`
to `0.080538` in training and `0.32971` to `0.13068` held out. Aggregate
held-out `C*dq` error fell from `7.9%` to `7.0%`. The four bent trajectories
improved, but their maximum tip differences remain `110–184 mm` with the
candidate. This is model-to-model matching, not experimental validation.
No section-specific nine-coefficient fit has been run.

A subsequent beta_v3 bound check kept all physical parameters and the pose
split fixed. Raising the cap from 3 to 5 gave
`[1.45574425,1.97259138,3.07558075]`; cap 7 gave the same result within
`9.816e-10`. Held-out normalized M RMS became `0.132242` (slightly worse
than `0.130683` at cap 3), while the four bent trajectories improved modestly
and minimum held-out mass-matrix `rcond` rose to `4.292e-4`. The exact wider
candidate passed derivative, skew, Cholesky, MATLAB/MEX, and sampled energy
checks. Both fits remain candidates; the runtime default is still all ones.
See `POINTMASS_BETA_TUNING_WORKFLOW.md` for the tradeoff and saved artifacts.

## Suggested next work

1. Reassess the physically relevant pose/velocity range and how to rank the
   cap-3 versus wider candidate, since held-out mass error and trajectories
   move in opposite directions. Neither has been adopted.
2. Keep the all-ones production default until a candidate is explicitly
   selected. Before adopting one, repeat derivative, skew, MATLAB/MEX, and
   energy checks at that exact beta over the intended operating range.
3. The CoG-specific section-energy fit is complete (see below). Its extra
   parameters are not adopted: held-out energy gains do not consistently
   carry over to mass, force, trajectory, and conditioning measures.
4. Both the fixed-beta location sweep and the subsequent location-specific
   energy fits are complete; retain their saved data as comparison baselines.
5. Consider denser output or solver-step-based damping integration if a tighter
   dynamic energy-balance result is needed. Strong-bend checks beyond the
   ±20 mm settings characterize numerical behavior, not physical validity.

Do not push, merge, delete backup stashes, or replace the standard reference
merely because this checkpoint exists. Inspect `git status`, active function
resolution, and the saved files in the current workspace before continuing.

## Fixed-beta mass-location sweep completed — 2026-09-17

User requested the same best betas held fixed while varying mass position.
This study uses the wider-bound candidate from `shared_beta_fit_bound5.mat`
at full precision: approximately `[1.45574425,1.97259138,3.07558075]`, shared
across sections. Production defaults remain all ones; no beta refit occurred.

- All three `cog_xi` values move together through `0.1:0.1:1.0`. Four cases:
  reference/asymmetric bend, each with free/differential input. All 40
  point-mass simulations and four fixed standard trajectories completed.
- Duration 5 s, 1501 recorded samples (300 Hz), ode15s tolerances 1e-8/1e-10,
  MaxStep 1e-3. Other physical parameters and production models unchanged.
- Saved signed coordinate histories and separate maxima for both coordinates
  per module, global section-tip vector/distance histories and maxima, and
  module kinetic/gravity/elastic/total energy differences. Also saved initial
  offset-corrected total energy and kinetic energy at the same standard state.
- Midpoint state reproduces the preceding comparison to 2.808e-17 m in q
  and 1.414e-16 m/s in dq at its original output times. Initial MATLAB/MEX
  RHS parity passed for all runs. All 31 mass checkpoints per run passed
  Cholesky; smallest sampled rcond 3.516e-6. No global SPD guarantee claimed.
- Midpoint minimizes worst-case tip peaks for all modules: 38.4, 109.9,
  176.9 mm. Module 3 local coordinate peak is slightly lower at xi=0.4
  (0.794 versus 0.848 mm), but its tip peak rises to 444.1 mm.
- Moving mass also changes gravity, which beta cannot correct at fixed q.
  Initial gravity mismatch at xi=0.1 is 22.8%/40.6% for reference/asymmetric
  bends; at midpoint 0.50%/3.30%; at xi=1 it is 29.1%/60.2%.

Read `cog_sweep_fixed_beta/ANALYSIS.md` for interpretation and 15 figures;
`ALL_MODULE_METRICS.md` contains all 120 case/location/module rows.
`cog_sweep.mat` retains full histories; `summary.json` retains summary metrics.
New reproducible helpers: `arm_module_observables.m`, `run_cog_sweep.m`,
`plot_cog_sweep.m`, `analyze_cog_sweep.m`. This is model-to-model comparison.
The later section kinetic-energy study below fits beta at each location
against this fixed-beta baseline on held-out data.

## CoG-specific section kinetic-energy fits — completed 2026-09-17

The user requested an energy-oriented fit following the fixed-beta sweep.
The image they supplied reproduced the fitting plan; the work was performed
against the actual MATLAB cores and saved split. New scripts:
`fit_section_energy_beta.m`, `validate_section_energy_beta.m`, and
`report_section_energy_beta.m`. Full output is in
`cog_beta_energy_fit/ANALYSIS.md`, with eight PNG/FIG figures,
`metrics.csv` (480 case/location/candidate/module rows), and MAT histories.

- For each shared mass position `xi=0.1:0.1:1`, fitted both three betas
  shared by sections and nominally nine section-specific betas to section
  kinetic energy at identical q,dq. All-ones and the previous fixed
  midpoint fit are the baselines. Used the original 12/8 training/held-out
  pose split and seeded, distinct velocity probes. Normalization comes from
  training section kinetic energies only, with a 5% cross-section floor.
- Section 1 has no upstream B, so beta_v1/beta_v2 cannot be identified or
  affect dynamics; they stay one. Only seven section-specific coefficients
  are active. Design ranks were [1,3,3] at all locations.
- Enforced the PSD beta cone and finite search bounds up to 100. At xi=0.1
  and 0.2, several coefficients hit 100. At xi=0.7 through 1.0, both new
  constrained fits return all ones. A beta correction constrained to add
  nonnegative kinetic energy cannot remove the excess point-mass energy
  there. CoG also changes gravity independently of beta.
- At xi=0.5, held-out normalized section kinetic-energy RMS: all ones
  0.08806, fixed midpoint fit 0.04563, new shared 0.03068, new sections
  0.03085. New section fit beats shared on training but slightly loses held
  out. Held-out normalized M RMS at the same xi: new shared 0.1233 versus
  sections 0.2772, a major tradeoff.
- At xi=0.5, maximum module-tip difference across four trajectories is
  176.9 mm for the fixed midpoint fit and 30.6 mm for section fit at 60 Hz.
  The new shared fit gives 39.9 mm. At xi=0.1–0.3, the section fit
  lowers held-out energy versus the new shared fit but worsens tip peaks.
  At xi=0.6, section-fit sampled fit-pose rcond falls to 2.036e-6 versus
  4.267e-4 for the new shared fit. Extra section coefficients are not adopted.
- Checked held-out `C*dq`, held-out normalized M, all six dM slices,
  symmetry and skew, fit-pose Cholesky, initial MATLAB/MEX RHS parity, and
  matched trajectories and module energy histories. Worst sampled dM
  relative error 3.345e-9; skew residual 2.504e-13. SPD is sampled only.
- Additional one-second MATLAB/MEX trajectory parity for section fits at
  xi=0.1, 0.5, 1.0 passed; largest scaled state gap 8.762e-10. These local
  timings are in `cog_beta_energy_fit/mex_parity.mat`.
- All-ones at xi=0.1–0.3 was too ill-conditioned for reliable trajectories:
  minimum fit-pose rcond at xi=0.1 was 1.071e-14 and first attempted run
  showed MATLAB/MEX mismatch. Its 12 trajectories are explicitly omitted
  (NaN in CSV), while its same-state fit/M/force metrics remain reported.

No production beta default, physical parameter, reference model, or MEX
binary changed during this fitting study. Do not treat model matching as
experimental validation. The new location-specific fits are exploratory,
especially at search bounds and at low sampled conditioning.

