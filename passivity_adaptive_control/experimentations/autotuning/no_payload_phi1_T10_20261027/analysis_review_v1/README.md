# First actual-arm framework recordings: offline review

Reviewed 2026-10-27. Source files are the three user-acquired actual-arm MATs in
the parent folder. The assistant ran offline analysis only; no model/hardware was
executed or raw recording changed. `../reviewFirstTrials.m` reproduces per-trial
reports into a fresh output directory. SHA-256 and Dataset provenance are saved.

## Verified protocol and inputs

Smoke: reference circles 20-50 s (3 turns), full log 0-60 s. Baselines: circles
20-140 s (12 turns), full logs 0-150 s. All use phi=1 rad during circles and
period=10 s. All logged tuning/switch snapshots match, with no changes during
tracking. Initial estimates [650 650 650 40 40 0 0 0 0 0] and history [0 0 0]
match across runs, supporting consistent software initialization. Physical
preconditioning/calibration/trim details remain operator-reported/unspecified.
No payload is the user's stated condition, not independently measured by the logs.

- lambda=[18;16.5], Ks=[15;10].
- gamma k=.5, d=.002, h=5, bias=1.
- All adaptation and variable-common allocation enabled; hysteresis freeze=false.
- Pressure activation=.2 bar, minimum=.2, maximum=3, slew=6 bar/s, nominal mean=1.6.
- Bounds lower=[-100;-500;-100;1;1;0;0;0;-8;-8], upper=[6000;6000;6000;250;250;20;20;20;8;8].
- Rate caps=[850;850;850;20;20;5;5;5;.2;.2]; scales=[1000;1000;1000;100;100;1;1;1;1;1].
- Deadband=.005 m/s; nominal residual tolerance=.02 N; maximum accepted dt=.075 s.

The .2 value verifies the allocator's logged activation threshold. It does not
identify the actuator's physical deadzone or independently verify every downstream
deadzone block. There are no measured chamber pressures/physical force errors.

## Results

| Metric | Smoke | Baseline 1 | Baseline 2 |
|---|---:|---:|---:|
| Tracking valid coverage | 98.67% | 99.17% | 98.83% |
| Adaptation active | 97.33% | 98.33% | 97.67% |
| Pressure limiting | 0% | 0% | 0% |
| Latched faults, whole recording | 0 | 0 | 0 |
| Tuning changes during circles | 0 | 0 | 0 |

Late scoring uses 40-50 s for smoke and 110-140 s (last three circles) for baselines:

| Metric | Smoke | Baseline 1 | Baseline 2 |
|---|---:|---:|---:|
| theta MAE (rad) | .46684 | .15202 | .14836 |
| phi MAE (rad) | .15048 | .06248 | .06234 |
| q2 RMSE (mm) | 4.253 | 1.568 | 1.529 |
| q3 RMSE (mm) | 4.338 | 1.889 | 1.857 |

Baseline final-lap theta MAE is .13558/.13221 rad. On laps 10/11/12 it continues
falling (.1650/.1555/.1356 and .1649/.1480/.1322), so 12 circles do not establish
parameter or tracking convergence. These logs support a repeatable fixed-duration
learning comparison, not a claim of steady-state identification.

Nominal allocator replay matches logged flags exactly: zero allocation and slew
limits in both baselines. Pressure ranges are .698-2.638 bar and .715-2.637 bar;
peak commanded rate is about 1.11 bar/s, below the configured 6 bar/s. Nominal
residual is numerical roundoff (~1e-14 N), not measured force error. Gate-active
reconstruction has zero disagreements. Occasional invalid samples and recovery
dt skips explain inactive adaptation; there is no persistent eligibility bottleneck.

Early bounds contacts: k2 lower .375%/.250%, k3 lower .917%/.833% of circle samples;
k3 rate cap .381%/.171% of active transitions. No parameter touches a bound in
110-140 s. Raising caps is not supported as the immediate need by these runs.

Estimated hysteresis weights reach roughly [10.11 9.25 9.64] and [10.05 9.20 9.49],
below their upper bound 20; damping reaches about [47.27 48.17]/[47.31 48.05].
They remain compensation estimates. K is indefinite during some early learning
(minimum eigenvalue about -345 N/m in baseline 1), but final eigenvalues are
positive (minimum 101.62/111.47 N/m). This does not establish physical parameters
or a closed-loop stability proof.

No v10-like sustained roughness escalation is evident over these 120 s circles.
Late theta-rate SD=.334/.327 rad/s vs early 20-40 s .358/.346. Residual uneven
motion remains: late near-stop fraction about 5.31%/4.75%, phi peak 1.111/1.105 rad.
Middle-window rate SD is lower (~.28), so this is not uniformly improving
smoothness. Angle-rate diagnostics use consecutive valid samples, not independent
velocity measurements. Minor invalid samples are retained in coverage statistics.

## Decision

The raw-log-to-evaluator smoke check passed and the two baselines are closely
repeatable. All required signals are present. Their evaluator classification is
`unconfigured_targets`, with success=false: targets, gamma bounds and formal
operator manifests have not been agreed, so these are not optimizer confirmation
passes or an authorized automated campaign.

Use 12 circles with the last three scored as a candidate comparison protocol if
the aim is learning speed by 120 s. Agree tracking/smoothness limits and a bounded
gamma search box before optimizer-selected trials. If the aim is converged
performance instead, the continued learning needs a separate longer assessment.
Do not mix those objectives or change gains mid-run. The .2 activation setting
must remain fixed within this campaign; comparisons to .8 activation need matched
trials and cannot be inferred from these three recordings alone.

Artifacts: per-trial reports/inventory/lap CSV/gates/tuning snapshots/evaluation
MAT, comparison.mat, and baseline_comparison.png. Raw MATs remain unchanged.
