# Shared translational beta fitting against the distributed model

Status: first shared-beta fit and staged beta_v3 bound check completed on
2026-09-16. This is a **model
matching exercise**, not validation against experimental arm measurements.
The distributed-mass implementation is fixed throughout. Section masses,
mass positions, geometry, gravity, stiffness, damping, actuation, six-coordinate
kinematics, and Taylor order are fixed to the saved
`comparison_after_C_correction.mat` setup. `cog_xi=[0.5;0.5;0.5]` and beta is
constant during each simulation. The three fitted beta values are shared by
all sections; no nine-parameter fit was performed.

## Reproduce in order

```matlab
validate_pointmass_beta        % gate: unequal non-unit dM, skew, SPD, MEX
fit = tune_shared_beta;        % saves shared_beta_fit.mat
report = validate_tuned_shared_beta; % saves shared_beta_validation.mat
timing = benchmark_shared_beta;       % saves shared_beta_timing.mat
```

The initial `validate_pointmass_beta` gate passed again before the fit:
maximum relative finite-difference `dM` error `6.053e-11`, skew residual
`2.327e-16`, and worst MATLAB/MEX scaled RHS difference `4.379e-13` across
unit and unequal non-unit beta. A positive-beta counterexample remained
indefinite, so positivity of entries was not used as an SPD test.

## Data split and objective

Training uses 12 fixed poses: straight, the saved reference bend, symmetric
`[+.01;-.01]` and `[-.015;+.015]` repeated over sections, and eight seeded
asymmetric poses in ±16 mm. Held-out validation uses eight distinct poses:
the symmetric ±20 mm bend, the previous asymmetric comparison bend, and six
seeded asymmetric poses in ±19 mm. The random seed and exact pose matrices
are saved in `shared_beta_fit.mat`. Velocities are absent from the M fit;
independent held-out velocities are used for `C*dq` evaluation.

For each pose, with standard mass matrix `Ms`, define

```text
D = diag(1/sqrt(diag(Ms)))
W = D/sqrt(norm(D*Ms*D,'fro'))
e(q,beta) = norm(W*(Mpoint(q,beta)-Ms)*W,'fro')
objective = mean_q e(q,beta)^2
```

The diagonal scaling makes each reference coordinate's inertia diagonal
equal before comparison. The Frobenius normalization gives every pose the
same reference-matrix norm. Thus a large raw matrix entry or one pose with
larger overall inertia does not dominate solely by scale. The reported RMS
is `sqrt(mean_q e^2)`, a dimensionless normalized matrix difference. Beta
enters `M` linearly, so the script precomputes its three matrix directions
at each training pose; held-out matrices are never in the optimization.

## Feasible set and guarantee

Bounds: `1 <= beta_v1 <= 3`, `0 <= beta_v2 <= 3`,
`1 <= beta_v3 <= 3`, with

```text
(beta_v2-1)^2 <= (beta_v1-1)*(beta_v3-1).
```

This is the positive-semidefinite condition on the Schur complement of the
coefficient matrix for the section translational energy,
`[1 1 1; 1 beta_v1 beta_v2; 1 beta_v2 beta_v3]`. It guarantees that each
section's *coefficient quadratic form* is nonnegative for arbitrary
upstream `A`, angular-position `B`, and local `P`. It does **not** guarantee
that the assembled six-coordinate mass matrix is positive definite or well
conditioned. Cholesky and `rcond(M)` are therefore checked at sampled train
and held-out poses. The chosen bounds keep this first search finite and
close to unit beta; they are modeling/search choices, not experimentally
identified physical limits.

## First fit and held-out checks

The shared fit converged to

```text
[beta_v1, beta_v2, beta_v3] = [1.46573213, 1.96512396, 3.00000000].
```

`beta_v3` reached its initial search bound and the cone boundary is active.
Do not interpret these values as final calibrated physical coefficients.

| Normalized mass-matrix RMS | Beta = 1 | Fitted shared beta |
|---|---:|---:|
| Training poses | 0.32941 | 0.080538 |
| Held-out poses | 0.32971 | 0.13068 |

Minimum sampled `rcond(M)` for beta one versus fitted was `2.509e-7` versus
`3.940e-4` in training and `4.302e-7` versus `4.020e-4` held out. Every
sampled mass matrix passed Cholesky. For independent held-out velocities,
RMS `||C*dq - Cstd*dq||` fell from `0.5072` to `0.4477` in generalized-force
units; aggregate relative error fell from `0.079` to `0.070`. Aggregate
relative error is `sqrt(sum ||difference||^2 / sum ||reference||^2)`, avoiding
unstable ratios at individual near-zero forces.

The four matched-input trajectories use the same saved parameters, initial
velocities, 0–5 s output times, and `ode15s` settings (`RelTol=1e-8`,
`AbsTol=1e-10`, `MaxStep=1e-3`). The differential input is
`[0.5;-0.5;0;0;0;0]` added to the saved zero input. Only beta changes.

| Start/input | Max coordinate difference, beta 1 → fitted | Max tip difference, beta 1 → fitted |
|---|---:|---:|
| Reference bend / free | 2.520 → 1.858 mm | 261.3 → 183.8 mm |
| Asymmetric bend / free | 1.337 → 0.998 mm | 144.0 → 110.2 mm |
| Reference bend / differential | 2.527 → 1.862 mm | 261.3 → 183.7 mm |
| Asymmetric bend / differential | 1.346 → 1.003 mm | 144.4 → 110.5 mm |

The remaining tip and coordinate differences are substantial. The shared
fit improves held-out mass, force, and trajectory measures, but it has not
made the point-mass model interchangeable with the distributed model.

## Staged beta_v3 bound check

The same objective, seed, training poses, held-out poses, and other beta
bounds were retained while raising only the upper bound on `beta_v3` from
3 to 5, then 7. `tune_shared_beta(5)` and `tune_shared_beta(7)` converged to
the same coefficients within `9.816e-10`:

```text
[beta_v1,beta_v2,beta_v3] = [1.45574425,1.97259138,3.07558075].
```

The beta_v3 cap is no longer active; the kinetic-energy cone boundary remains
active. The cap-5 fit is saved in `shared_beta_fit_bound5.mat`; the cap-7
replication is in `shared_beta_fit_bound7.mat`. Held-out comparisons:

| Measure | Initial cap 3 fit | Wider cap 5/7 fit |
|---|---:|---:|
| Training normalized M RMS | 0.080538 | 0.079922 |
| Held-out normalized M RMS | 0.130683 | 0.132242 |
| Largest held-out normalized M error | 0.265648 | 0.271305 |
| Held-out aggregate C*dq error | about 0.070 | about 0.070 |
| Minimum held-out rcond(M) | 4.020e-4 | 4.292e-4 |
| Reference bend/free max tip difference | 183.8 mm | 176.9 mm |
| Asymmetric bend/free max tip difference | 110.2 mm | 107.3 mm |

The widened solution improves trajectory and conditioning measures, but
slightly worsens the mass-matrix objective on held-out poses. This is a
tradeoff, not an unqualified win. Neither candidate has replaced the
all-ones production default. The widened trajectory/force results are saved
in `shared_beta_validation_bound5.mat`.

At the exact widened coefficients, `validate_pointmass_beta(betaMatrix)`
passed every sampled derivative slice (worst relative finite-difference
error `6.239e-11`), symmetry, `Mdot-2*C` skew symmetry, Cholesky at the three
test poses, and MATLAB/MEX RHS parity (`4.379e-13` worst scaled error). The
reference-bend zero-input fitted trajectory lost `1.17325 J` of mechanical
energy monotonically at the 301 output samples. Its maximum damping-balance
residual was `6.136e-4 J`, or `0.0523%` of the energy drop; this result is
limited by output-sample quadrature. See
`pointmass_beta_candidate_validation.mat` and
`arm_energy_check_candidate.mat`.

The next distinct experiment can vary the selected backbone mass location.
First hold beta fixed while moving `cog_xi` to measure location sensitivity;
then refit beta at each location on the same split. Keep those two effects
separate. CoG location was **not** changed in this bound check.

## Matched MATLAB/MEX timing check

`benchmark_shared_beta.m` solves the same first 1 s of the reference-bend
problem at the saved output times, after an RHS warm-up, with plotting and
code generation excluded. Three integrations per mode gave these medians:

| Execution path | Median integration time |
|---|---:|
| Standard MATLAB | 4.9827 s |
| Standard MEX | 0.1860 s |
| Fitted point mass MATLAB | 1.7369 s |
| Fitted point mass MEX | 0.0867 s |

Maximum state differences between matched MATLAB/MEX solutions were
`1.734e-12` for standard and `2.567e-16` for fitted point mass. This is a
small, one-machine timing check, not a general algorithmic benchmark.

## Decision about nine section-specific coefficients

At this first shared-fit checkpoint, no nine-coefficient fit was run. The shared fit does improve held-out scores,
while its held-out RMS (`0.131`) exceeds its training RMS (`0.081`) and
`beta_v3` hits a bound. Those observations call for checking pose coverage,
search bounds, and robustness before adding six degrees of freedom. A
nine-parameter model should be accepted only if a separately trained fit
improves *held-out* mass and `C*dq` scores and the four held-out trajectories
without unacceptable conditioning. Its improved training score alone would
not justify the complexity.

## Fixed-beta location sensitivity — completed 2026-09-17

The subsequent user-requested sweep holds the full-precision cap-5 candidate
fixed and sets all three `cog_xi` entries together to each of `0.1:0.1:1.0`.
All 40 matched cases completed (reference/asymmetric, free/differential),
with unchanged physical parameters and standard reference. No fitting was
performed. Production defaults are still all ones.

See `cog_sweep_fixed_beta/ANALYSIS.md`, its 15 PNG/FIG plots, full module
table, JSON summaries, and MAT histories. Reproduction scripts are
`run_cog_sweep`, `plot_cog_sweep`, and `analyze_cog_sweep`.

The midpoint minimizes worst-across-case tip peaks and kinetic-energy RMS
for every module on this grid. Module 3 coordinate peak alone is slightly
smaller at 0.4. Moving mass away from the midpoint increases gravity mismatch;
beta cannot alter gravity at the same q. The later per-location fitting used
the established training/held-out split and retained these fixed-beta
results as its baseline. Do not interpret simultaneous movement of all mass
positions as isolated per-section sensitivity or as experimental validation.

## Later section-energy fit at each CoG — completed 2026-09-17

See `cog_beta_energy_fit/ANALYSIS.md` for the full objective, all ten
coefficient matrices, held-out energy/M/force comparisons, matched-trajectory
plots, and numerical checks. The executable sequence is
`fit_section_energy_beta`, `validate_section_energy_beta`,
`check_section_energy_mex_parity`, then `report_section_energy_beta`. This
study uses the same 12/8 pose split,
distinct training/held-out velocity probes, and training-only energy scales.

Section 1 has only beta_v3 identifiable; its beta_v1/beta_v2 are inactive
and fixed at one. Thus the nominal 3-by-3 parameter matrix has seven active
values. The section fit improved training energy and some held-out energy
scores, but did not consistently improve held-out M or matched trajectories
versus a new shared fit at the same location. At xi=0.5, held-out normalized
kinetic-energy RMS was 0.03068 shared versus 0.03085 section-specific;
normalized held-out M RMS was 0.1233 versus 0.2772. At xi=0.1–0.3,
section-specific energy gains came with worse tip peaks. At xi=0.6,
conditioning worsened sharply for the section fit. No section-specific
coefficient set was adopted. The all-ones baseline was too ill-conditioned
to simulate reliably at xi=0.1–0.3, and those 12 trajectories are explicitly
missing from the report. Production defaults remain all ones.
