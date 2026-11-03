# First trained generalized-force residual model

Completed November3,2026, offline in MATLAB R2025a CPU. This is the canonical
first force-model release. No hardware/Simulink execution or controller edits.
Old position-correction models and raw recordings are unchanged. No commit/push.

## What was trained

    M(q)ddq+(C(q,dq)+D)dq+Kq+bound(q)+G(q)=alpha*A*p+r_tau

q=[l2;l3],l1=-l2-l3. The NEW GRU outputs two generalized residual forces in
NEWTONS, inserted into the mechanical force balance and integrated to obtain
q/dq. There is no learned position jump. This is the force-model placement
requested by the user, not relabeled weights from the earlier position GRU.

Architecture: GRU32 (reset-after gates, stable tanh equivalent of sigmoid),
linear force head plus affine input skip, smooth componentwise80N output bound.
Inputs: normalized q2,q3,v2,v3,p1,p2,p3,previous physical dt. Ten observed history
frames initialize hidden state; subsequent rollout inputs are integrated q/v
and known commands/time, not future measured motion. Inputs/normalizer contract
and weights are stored together in force_model.mat.

Nominal constants remain fixed: L=.17411m,r=.013m,mi=.1kg,vertical gravity,
K=650*[2 1;1 2],D=40I,three-actuator soft bounds at+-35mm withstrength1e6,
mu2000. Corrected baseline mode1,zero-activation pressure map with
alpha=pi*.0065^2*1e5 N/bar and A=[-1 1 0;-1 0 1]. No additional mechanical
parameter fit, payload estimator or identified chamber-pressure model.

## Data use and model selection

Only October19 and October27 raw recordings supplied training/development
data. The actual ten old-source SHA256 hashes were checked. Day identities
and whole trials are retained. Independently within each day/channel stratum,
the middle listed trial is reserved for internal selection; other trials train.
No random adjacent-row train/selection split. Existing invalid/gap guards,
causal v filter and .60s warmup are retained from the audited data preparation.
Neither legacy acceleration nor inverse-pressure residual labels are targets.

Training normalizer:59,630 mature TRAIN rows from the older days only.
Gradient-fitting windows:2,883 (cap360/run),24 integrated steps each.
Internal selection:160 windows,40 steps each,20windows across eight usable runs.
Overlap between windows does not create independent experimental replicates.

October27 StepWave contributes only3 gradient windows and no selection windows;
October27 combinePressure contributes zero complete windows. This follows the
unchanged timing/history rules, not source exclusion after examining fit errors.
All ten sources are audited, but not all recordings supply substantial usable
training coverage. See ../dataset_20261103_v1/coverage.csv.

Two completed fits were compared ONLY using older internal selection scores:

| Fit | Learning rate | Selected epoch | Old selection macro length RMSE |
|---|---:|---:|---:|
| train_20261103_v2 | .002 | 5 | 3.0332mm |
| train_20261103_v3 | .0003 | 38 | 2.1103mm |

The first fit worsened after epoch5; its best checkpoint remains archived.
The lower-rate fit ran40epochs (2,440updates, batch48), gradient clipping1,
selected epoch38 by the older trial score. The final model is an exact-byte
copy of that checkpoint release; no retraining after November2 evaluation.
Old matched baseline macro scores on the same40step selection windows are
nominal6.8260mm and constantvelocity1.0377mm. The force model improves the
nominal model but does not yet beat constant velocity on that aggregate.

## Objective and numerical implementation

Length-output error through24 integrated steps, normalized by.005m, plus
.001 times squared residual-force rate normalized by100N/s, plus.0001 times
squared residual magnitude normalized by80N. Differences are between adjacent
times in each valid sequence, never unrelated shuffled batch examples.
This implements the retained smoothness proposal inside training. No output
low-pass or plot smoothing is applied. Penalty strengths were fixed before
November2 evaluation. No measured force or differentiated acceleration is used
as force ground truth.

Analytic mass/gravity polynomials are an exact rational rearrangement of the
EXISTING nominal PCC Taylor model. They are not learned or calibrated. Batched
coefficients match armS_single_core to1.56e-15 relative error in the test suite.
Coriolis terms use the same Christoffel convention. Only fixed vertical gravity
is implemented in this training module; arbitrary base attitude is not qualified.

The nonlinear model is integrated with second-order L-stable SDIRK, four
substeps/observed interval, three fixed-point stage updates and five Newton
iterations for bounds. It replaces the old locally frozen affine transition.
Numerical solves are approximate; no claim of exact continuous-time integration.

Neural gradients use MATLAB autodiff. To avoid an impractically large implicit
solver graph, local q/v/force solver Jacobians use batched central differences
and reconnect to the neural graph through an exact-forward affine tangent.
Thus forward trajectories remain exactly the implemented nonlinear solver;
this is not stopped-gradient training. Derivatives are first-order numerical
approximations, not higher-order automatic differentiation. Perturbations:
q1e-6m,v1e-5m/s,force1e-4N.

Verification: residual force/pressure-map equivalence checks units/sign/RHS;
direct force-through-solver AD versus finite differences relerror6.19e-10;
network-through-rollout gradient probe relerror1.37e-10; changing future targets
does not change predictions. SDIRK versus independent ode15s test positions at
four substeps differ by.0336/.0297/.0049mm for three chosen transients. Doubling
substeps on one older selection window/run gives path differences.002-.076mm,
maximum point difference.174mm. These are numerical checks, not real-model
accuracy or force-identification evidence.

## November2 held-out evaluation

The user reserved November2 for evaluation. It was loaded by the evaluator
only AFTER selecting/freezing the release and writing its checksum marker.
It supplies no fitting, normalizer statistics, early stopping or candidate
selection. However, these recordings were inspected/evaluated previously in
the project, so they are not an untouched blind confirmation set.

All five recordings: StepWave,SineWave4,combinePressure,circleSweep,helix.
50,078 original strict eligible anchors; no relaxed masks or invalid filling.
Matched comparisons use the corrected nonlinear nominal with zero residual,
constant velocity and the new force GRU. Previous position-model tables use
different propagation and must not be silently merged with this table.

Mean per-run vector length RMSE over l2/l3,mm:

| Predictor | One native sensor step | .5s path | 2s path | 5s path |
|---|---:|---:|---:|---:|
| Corrected nominal, zero residual | 4.7352 | 6.9508 | 7.0438 | 6.6372 |
| Constant velocity | .06235 | .18575 | .89271 | 3.4062 |
| New force-residual GRU | 2.1510 | 1.9373 | 1.8375 | 2.1267 |

Each rollout horizon has60 matched initial windows/model,12/run, known
recorded pressure/query sequences, observed initial history only, no measured
state resets. Start pools can differ across horizons. Run-level path error is
the square root of mean per-window path MSE; final macro averages five runs.
540 model/window propagations complete with zero inherited coordinate-domain
guard failures. One-step examples also have zero numerical-guard failures.
This guard is not the physical actuator bound or a stability/safety theorem.

The new model materially improves the nominal baseline. At5s it improves the
constant-velocity macro score by37.6%, but constant velocity remains much better
one-step/.5s/2s. In particular, the force model's2.15mm one-step error is too
large to call the initial model ready for inverse/PAC/PDFB use. Do not present
only its favorable5s comparison or claim it matches the old position GRU.
Initialization, state/velocity consistency, capacity and optimization remain
possible development issues; their causes have not been isolated here. Later
model revisions must use older development data for tuning, not these scores.

One-step learned residual peaks range9.7-18.0N across runs; none approaches the
80N head bound. These are effective compensation forces, not measured or
uniquely identified physical forces. Commands are not chamber pressure; lengths
are inferred from orientation; physical NI-write times and exact actuator lag
are unavailable. Setup/calibration equivalence between days is not established.

## Use and limitations

The frozen model is stored in force_model.mat. Required source is in
source_snapshot. Evaluation/plots/numerical inputs are in
../evaluation_20261103_v1; its report and FIG/vectorPDF/300dpiPNG retain both
favorable and unfavorable results. Old datasets/position weights unchanged.

Future inverse pressure calculation must account for the residual's dependence
on current pressure/history and choose a common-pressure policy for three
actuators/two coordinates. Subtracting a command-dependent residual without
solving that dependence is not a complete inverse. No inverse-pressure or
controller deployment is included. Bounded outputs do not prove passivity or
closed-loop stability. More development/validation is needed before arm use.
