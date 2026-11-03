# Generalized-force residual model

Current November3 checkpoint: the three-seed1x64/1x128/2x64 comparison favors
1x128 on the guarded short-horizon score. Candidate and complete report:
capacity_comparison/campaign_20261103_v1. Keep fine_tuning/campaign_20261103_v1/
release/force_model.mat as the64-unit reference because128regresses long circle
free predictions. Weights are frozen; controller integration remains pending.
Use streaming/sgCell with cfg.widths, not the legacy hiddenUnits-based cell.
See project_docs/RESIDUAL_CONTROLLER_HANDOFF.md before integration.

Earlier completed stages and their separate state contracts follow.

November3 streaming follow-up sequence is complete under streaming/REPORT.md:
six architecture/horizon fits, burn-in comparison and three seeds per finalist.
New observed-state TBPTT contract and separate free prediction; selected GRU64
with .5s horizon. Frozen November2 one-step .563+/- .068mm and .5/1s paths
.361/.427mm improve nominal; CV remains better through1s. Offline bounded
pressure inversion works on model-generated targets, with timing outliers and
no physical-force/controller qualification. New research checkpoint in
streaming/release_20261103_v1. This parent release and its code/results remain
unchanged; do not mix the two model state/parameter schemas.

Completed first fit/release/evaluation November3. Canonical release:
release_20261103_v1/REPORT.md and force_model.mat. Selected train_v3 epoch38
using older trials only. November2 held-out evaluation is in
evaluation_20261103_v1. It improves nominal prediction but remains weak at
short horizons: one-step2.151mm,2s1.837mm,5s2.127mm macro length RMSE.
Constant velocity is better short-term; force GRU better at5s. No inverse or
controller qualification. Full reports retain the limitations and negative data.

Implementation begun November3,2026 following explicit user authorization to
train on October19+October27 and use November2 only for evaluation. Earlier
position-correction models remain unchanged. No Simulink/hardware execution.

Model: M(q)ddq+(C(q,dq)+D)dq+Kq+bound(q)+G(q)=alpha*A*p+r_tau.
New GRU32 plus affine skip outputs two forces in N, bounded smoothly at80N.
Observed q/v/command/dt history initializes hidden state; subsequent rollout
uses integrated q/v, never future measured response. Current pressure is an
input, so future inverse-pressure use will require solving that input dependence;
simply subtracting an unknown command-dependent residual is not a complete
inverse implementation. No inverse/controller integration is supplied here.

Fixed corrected all-three-actuator bound nominal (mode1), L=.17411m,r=.013m,
mi=.1kg,K=650*[2 1;1 2],D=40I,vertical gravity,zero activation pressure map.
forceResidualPhysics is an analytically reorganized EXACT original PCC Taylor
mass/Coriolis/gravity model, not a fit/surrogate. Its coefficients come from
rational integration; numerical agreement with armS_single_core is1.56e-15.
Second-order L-stable SDIRK integrates the dynamics with nonlinear implicit
stages. This replaces the locally frozen affine transition used by old screens.

Training loss: integrated length error over24steps, residual-force-rate penalty
on adjacent physical-time samples, small force-magnitude penalty. No noisy
acceleration/force labels. Force-cell and GRU gradients use MATLAB autodiff;
solver gradients use batched central-difference Jacobians reconnected to the
graph with exact forward values (first-order training derivatives). This avoids
a prohibitively large implicit-solver autodiff compilation graph. It is not
stopped-gradient propagation. Finite-difference network/solver checks pass.

Whole-trial internal selection is drawn independently within each older day.
No November2 data, normalization or metrics enter fitting/early stopping.
Existing gaps/maturity/causal-filter policy retained. Day2 pressure recordings
have few/no complete windows; coverage.csv makes this explicit. New normalizer
uses older TRAIN rows only; old cached nominal labels are never training targets.

Primary commands from repository root (fresh output folders required):

    addpath(pwd,fullfile(pwd,'residual_dynamics_learning','force_model'));
    cfg=forceResidualConfig();
    prepareForceResidualData(pwd,dataFolder,cfg);
    probeForceResidualTraining(dataFolder);
    trainForceResidualModel(pwd,dataFolder,trainingFolder);

Canonical preparation: dataset_20261103_v1. Initial backward-Euler checks_v1
failed transient accuracy, marked DO_NOT_USE. SDIRK checks_v2 pass coefficients,
force-placement/units/sign, nonlinear-ODE convergence and force gradients.
The abandoned full-autodiff accelerator compilation probe was stopped; it
performed no fitting. Raw files, old predictor weights and hardware untouched.

During training, training_history.csv and checkpoints are written every epoch.
force_model.mat and FROZEN_BEFORE_NEW_DATA.md are created only after selecting
the checkpoint on older internal trials. External evaluation must happen after
that freeze; no test-informed refit. November2 has been inspected historically,
so call it held-out evaluation, not an untouched blind confirmation experiment.
