## 2026-11-03 - HIL paper reviewed: proposed supervised residual-model loop

Reviewed Marquez et al2023 local HIL/TD3 paper. It learns a locomotion POLICY,
not residual dynamics; its reported gains do not establish model accuracy.
Proposal: reuse hardware-data/replay/batch-update structure for supervised GRU
model learning, contiguous recurrent sequences and held-out prediction tests,
with PAC/PDFB retained. Start between trials; neural weights frozen within each
comparison. No actor-critic or ROS rewrite required for this proposed version.
Detailed source/limitations/sequence: project_docs/HIL_RESIDUAL_MODEL_LEARNING_NOTE.md.
Primary next step remains consistent PAC/PDFB model integration/comparison;
no implementation, hardware run, new fit or automatic deployment authorized.

## 2026-11-03 - Three-seed GRU capacity comparison completed

User authorized quick offline1x64/1x128/2x64 comparison. Three seeds113/211/307, matched original1200-update stage plus1200 old-only continuation updates; unchanged physics/loss/masks/normalizer.
Six new stage1 and eight new continuation fits; exact existing64-unit fits reused by hash. Old TRAIN gradients only; old/Nov3 holdouts select checkpoints and architecture; Nov2 evaluation excluded from selection.

| Architecture | Parameters | Mean score (mm) | Sample SD | Seed113 | Seed211 | Seed307 |
|---|---:|---:|---:|---:|---:|---:|
| 1x64 | 14162 | 0.58531 | 0.02495 | 0.58716 | 0.55949 | 0.60928 |
| 1x128 | 52882 | 0.48684 | 0.02226 | 0.46154 | 0.49553 | 0.50344 |
| 2x64 | 38930 | 0.53814 | 0.05209 | 0.50814 | 0.50799 | 0.59828 |

Predeclared winner: **gru128**. Frozen representative candidate: candidate_force_model.mat.
SHA256 `8f2e0bb55ea9399b62abd5abe0b8f6dd15ab0470a400ad3e593e1edcd8f9f545`.
The existing fine-tuned GRU64 research release remains unchanged. This candidate is an offline result, not deployed or controller-qualified.

Representative seed113 fixed in advance, not chosen by best seed. Current fine-tuned GRU64 release retained; separate frozen candidate archived.
GRU128 improves 4/10 fixed long examples relative to GRU64, while others worsen; these correlated examples are not independent statistical trials. Circle free RMSE worsens on Nov2 0.623->0.930mm and Nov3 1.005->1.366mm. It improves several sine/helix examples. The short-horizon winner is therefore a candidate for measured-state model augmentation, with mixed long-horizon accuracy. Retain the current GRU64 reference; do not infer improved physical tracking from short prediction scores.
All nine cases have common-anchor evaluation; three representatives have whole-record one-step, fixed long motion predictions, and208-case bounded numerical inverse checks.
Report: residual_dynamics_learning/force_model/capacity_comparison/campaign_20261103_v1/REPORT.md.
Figures FIG/PDF/PNG, predictions/checkpoints/source/hashes and verification retained. No physical force/closed-loop qualification; fresh confirmation and PAC/PID FB integration pending. No hardware/SLX/migration/Git changes.

## 2026-11-03 - Force GRU64 fine-tuning completed

User authorized fine-tuning. Two matched1200-extra-update arms from exact frozen
GRU64/B seed113: old October TRAIN only;16old/16Nov3 stream replay. LR1e-4,
fresh Adam, unchanged architecture/.5s horizon/burn-in/loss/filters/normalizer.
Nov3 whole-trial split:28train/9selection,2264new train chunks,80new holdout
anchors; original128old selection anchors/splits/59,630-row normalizer unchanged.
No interpolation. Both candidates improve old and new internal scores.
Predeclared .5old+.5Nov3 score (each .7oneStep+.3halfSecond), old regression
guard10%, chose old_only BEFORE Nov2 evaluation: composite.71126->.58716 versus
mixed.59696. Old internal.54411->.41699; Nov3.87840->.75734. Mixed improves
new score more (.68348), but old score less (.51044). Keep both candidates;
this is one optimization seed, not evidence of universal old-only superiority.

Selected research model: residual_dynamics_learning/force_model/fine_tuning/
campaign_20261103_v1/release/force_model.mat, SHA256 1c38edc36ba6ff4a41da8f974a88ca1c5ea775655426da91124a35b255a9c870.
Full report: same campaign/REPORT.md. Selected gradients use older TRAIN data;
Nov3 holdouts were used for checkpoint selection, not gradients in that arm.
Nov2 external one-step macro vector RMSE .53588->.36143mm; mixed.44580mm.
Nov3 heldout one-step .93307->.75485mm; mixed.69078mm. Values use fixed original
seed113, not earlier three-seed mean. CV remains better one-step.
Full circle q-vector RMSE: Nov2 .9440->.6229mm; Nov3 1.2636->1.0047mm.
Selected model improves allfive Nov2 motion examples, but some Nov3 long
windows regress slightly; mixed replay has larger long-window regressions,
despite better new one-step scores. Allnegative results retained. Nov3 helix
comparison uses HELDOUT trial2 (4.0294mm original), not prior TRAIN trial1
(2.8787mm); do not mix those comparisons.30 complete long model/window paths.

75one-step and30long-window table rows independently recomputed<5e-15mm;
source/model/data/splits/normalizer and original releases preserved. Selected
inverse self-consistency208/208, max equation error3.57e-5N, state-commit and
bounds pass; no physical force/control qualification. Seven FIG/PDF/PNG bundles,
full CSV/MAT/checkpoints/optimizer states/source hashes archived. Current
controller/PAC adaptation integration plan still pending; no SLX/hardware,
migration, commit or push performed. No fitting or checkpoint selection on Nov2.

## 2026-11-03 - Primary goal clarified: model augmentation for PAC and PID FB

User explicitly intends to use the improved dynamics MODEL with existing PAC
and PID feedback controllers and test whether tracking improves because the
model is more accurate. This is the primary integration target. Detailed plan:
project_docs/RESIDUAL_CONTROLLER_COMPARISON_PLAN.md. Four primary conditions:
PAC nominal / PAC augmented; PID FB nominal / PID FB augmented. Compare model
off/on within each controller with matched gains, references, filters, pressure/
common-pressure allocation, load, calibration, initialization and duration.
Freeze neural weights; recurrent state evolves. Match PAC adaptation settings
and reset policy. Report tracking/convergence/smoothness/limiting and independent
model accuracy, not just one-step prediction error. A separate retuned-controller
study can follow but cannot isolate model benefit by itself.

Before integrating, reconcile the residual's fixed nominal K0/D0/B*p definition
with each actual controller. PAC already estimates K/D/hysteresis/bias; derive
consistent compensation/adaptation before stacking the neural force term.
Pressure-dependent residual requires bounded inverse calculation with fixed
pre-current hidden state during candidate search and one selected-state commit.
Verify actual PID/integral path (repository also uses PDFB/PDI/FB names), actual
pressure/common-pressure operating coverage and runtime settings. Original
passivity is not automatically inherited. This is a recorded design/test plan,
not a completed controller derivation, SLX change or hardware authorization.
No continuation training, controller edits, hardware, commit or push this turn.

## 2026-11-03 - November3 frozen-model plots and continued-training recommendation

User requested same circle/sine/helix plots for November3 and asked which newer
day to use for continued training. Complete seven FIG/PDF/300dpiPNG bundles,
MAT/CSV/source/hashes: residual_dynamics_learning/figures/
motion_predictions_20261103_v1/frozen_model_predictions_v1/REPORT.md.
Same frozen FORCE GRU64/B seed113, no fitting, changed normalizer or imputation.
Separate measured-state one-step and free recorded-input predictions, nominal
shown separately. New helper reproduces ALL saved November2 circle one-step
anchors/outputs<1e-12m; all15 initial free steps match measured replay<1e-10m.
All15 model/window free paths complete, four helix labels stay NaN/unscored.

Free q-vector RMSE mm, GRU/CV/nominal: circle1.264/1.948/11.625;
sineP1 1.718/9.767/6.742; P2 1.379/15.735/9.700;
P3 1.688/19.562/7.670; helix2.879/13.169/11.808. Correction to an intermediate
commentary range: sineP2 is1.38mm, not1.69mm. Selected source intervals:
circle phi_reference1,159-169s; sine[6.6-19.2,102-114.6,197.4-210]s;
helix6.6-64.95s (59.047physical seconds). These are chronological examples,
not independent repeated trials. One-step CV still better; reconstructed PCC
XYZ is not sensor XYZ. No controller tracking or unique physical deadzone claim.

Logged Nov3 circles10s versus Nov2 20s; sine~12.575s versus~31.4s; helix10s/turn
and peak requested phi2.5 versus20s/turn,phi1.5. Speed, amplitude and date all
change, so larger circle/helix errors cannot be attributed to speed alone.
Actual59,630 old TRAIN feature rows audited: some new circle/helix q/v extend
beyond old observed extrema (~2-3.5% per affected feature); q2 minimum old
-14.98mm versus newcircle-22.18mm/helix-20.73mm, positivev2 old.01322m/s versus
new~.0209m/s. Large training z-scores alone do not prove outside-range input.

Recommended NEXT offline experiment, not executed: fine-tune frozen GRU64 on
old October19/27 TRAIN plus November3 TRAIN, retaining older data replay and
whole-trial/condition internal holdouts. Keep November2 evaluation-only for
frozen before/after comparison, acknowledging it is previously inspected, not
blind confirmation. Lower learning rate, same architecture/loss/history and
original normalizer initially; declare fresh optimizer state (release weights
are not a full Adam resume). Compare extra-old-only updates against old+Nov3
at matched budgets, to separate new data benefit from just optimizing longer.
Do not combine imputation/filter/capacity changes in the first continuation
experiment. If Nov2 enters training it loses held-out status; use another
cohort/fresh day for independent evaluation. A future fresh repeat/day is needed
for final confirmation. Neither current day establishes2s-circle performance.
Sources/model/raw/old Nov2 plotting archive unchanged; no hardware, SLX,
training, migration, Git commit or push performed.

## 2026-11-03 - Requested November2 motion plots; November3 NDI gap audit

User clarified main plotting cohort is November2. Frozen FORCE GRU64/B seed113,
not the older position-correction GRU. Seven FIG/PDF/300dpiPNG sets plus full
MAT/CSV/source/hashes: residual_dynamics_learning/figures/
motion_predictions_20261102_v1/REPORT.md. Circle phi_reference1,259-279s;
one sine period/channel; helix6.6-124.95s (119.85physical seconds), loading and
unloading after observed warmup. Separate exact saved measured-state one-step
and free recorded-pressure predictions, nominal shown separately. Fifteen free
model/window paths complete; first steps reproduce frozen results<1e-10m.
Length-vector RMSE mm, GRU/CV/nominal: circle.944/2.253/11.481;
sineP1 2.336/9.481/6.070; P2 1.778/17.480/9.656;
P3 1.796/21.731/7.559; helix.763/4.740/8.972. Chronological examples exceed
.5s training and prior1s evaluation; this long-window GRU advantage coexists
with CV winning short-term. No universal or controller tracking claim. PCC XYZ
is reconstructed from inferred q, not independent sensor XYZ. Six missing
labels remain NaN/unscored; free model carries on known pressure/time without
measurement injection. One-step retains original reset/recovery warmup.
Source/model/raw/export hashes unchanged; previous release preserved.

All five November3 files exported/audited separately, no fitting or new model
score evaluation.57 invalid NDI rows/55 episodes:53 one-frame, two two-frame.
Separate explicit preview fills53 one-frame gaps at<=.12s valid-to-valid
physical span, keeping raw values/validity and imputed masks, q1=-q2-q3; both
longer episodes remain missing. Artificial-mask proxy reconstruction vector
RMSE.040-.098mm acrossfive records (max individual one-frame error.4855mm).
Actual dropout truth remains unknown; overlapping masked blocks are correlated.
Artifacts/source/CSV/FIG/PDF/PNG in figures/motion_predictions_20261103_v1/
short_gap_audit/REPORT.md. Preview never entered training or frozen evaluation.
Proposed future training experiment: bounded offline inputs, explicit imputed
flag, measured-only targets, long-gap resets, unchanged measured-only/causal
validation. Two-sided interpolation uses the future endpoint: masking only
filled targets does not remove leakage through recurrent inputs/state.
Alternative: gap-aware physics carry with loss only at real observations.
No interpolation policy deployed, retraining, SLX/hardware, migration or commit.

## 2026-11-03 - Streaming force-GRU sequence completed

Explicit sequence executed: six architecture/horizon fits, one burn-in comparison,
four additional seed fits. New observed-state TBPTT primary path and independent
short free-prediction state; no future target inputs. Old October19/27 whole-trial
splits and normalizer unchanged; all finalists frozen before November2 replay.
Selected B,64-unit single layer,0.5s horizon,.5s burn-in;
fixed representative seed113, chosen config by three-seed mean older score.
Two-layer32+32 did not win under the common1200-update budget; not a general
depth result. November2 macro one-step 0.5629 +/- 0.0684mm; .5/1s path
0.3613 +/- 0.0148/0.4272 +/- 0.0183mm. Matched CV one-step
0.0621mm. CV remains better at all tested horizons through1s. All numerical
guard failures 0; per-phase/coordinate/
force-rate/seed results retained. New mode differs from parent release, so do
not attribute all changes to width. No post-test tuning or new November3 pooling.
Offline bounded, pressure-dependent model inversion: 128/128 known
targets feasible, max equation residual 4.743e-05N, median/p95 latency
5.19/13.93ms excluding acquisition/control scheduling.
Trial hidden states remain fixed; one selected-state commit passes. This is
self-consistency, not physical force or closed-loop tracking qualification.
Code/state/gradient/horizon/split/hash checks and independent saved-metric
recompute pass. Prior release unchanged. Reports/weights/CSV/MAT/FIG/PDF/PNG:
residual_dynamics_learning/force_model/streaming/REPORT.md. No SLX/hardware,
migration, commit or push.

## 2026-11-03 - First trained generalized-force residual, two old days only

Explicit user request executed. Canonical model: force_model/release_20261103_v1
under residual_dynamics_learning; detailed report there, held-out results under
force_model/evaluation_20261103_v1. NEW GRU32 outputs two generalized forces in
N inserted into Mddq+(C+D)dq+Kq+bound+G=alpha*A*p+r. No position jump or relabeling
of old position-network weights. Fixed corrected-mode1 nominal, no calibration
detour. Ten initial observed frames, then integrated q/v and recorded commands.

October19+27 source hashes checked, whole-trial internal holdouts retained on
both days.59,630 mature TRAIN rows normalize inputs;2,883 training windows
(24steps),160 selection windows(40steps). Day2 StepWave3/0 train/selection,
combined0/0 from unchanged gap guards; other8runs360/20 each. No new-data fit,
normalizer or early stopping. Two completed rates.002/.0003 compared using old
macro score; selected.0003 epoch38 (2,318 updates) from40epochs,2.1103mm versus
older nominal6.8260/CV1.0377. Earlier high-rate fit diverged after selectedepoch5;
best weights retained. Incomplete train_v1 marked DO_NOT_USE.

Loss is integrated q output error plus modest force-rate and magnitude terms;
no measured-force or derived-acceleration labels. Exact analytic PCC mass/G
reorganization matches existing coefficients1.56e-15. L-stable SDIRK nonlinear
integration, numerical local solver Jacobians with exact forward values plus
neural autodiff. Direct/neural gradient finite-difference checks pass; actual
4-substep ODE fixture errors<.034mm; older-trial doubled-substep path differences
.002-.076mm, maxpoint.174mm. Numerical checks do not establish physical accuracy.

Exact model SHA256 b5cfd7d0cd13b2b9067a35d810b8ae166e032e96421b6a30c55d703039af268f
frozen before loading November2 sequences in evaluator. Five recordings,
50,078 one-step anchors,540 matched .5/2/5s model/window propagations, no numerical
guard failures. Mean per-run vector q2/q3 RMSE mm:
one-step nominal4.7352/CV.06235/forceGRU2.1510;
2s paths nominal7.0438/CV.8927/forceGRU1.8375;
5s nominal6.6372/CV3.4062/forceGRU2.1267.
New force model improves nominal; CV still better short-term, forceGRU37.6%
better at5s. Do not overstate or omit poor one-step performance. No refit after
new-data scores. Prior-inspected November2 is held-out evaluation, not blind.

One-step residual peaks9.7-18N, no approach to80N head bound. These are effective
model compensation, not uniquely identified/observed physical force. Pressure
proxies, inferred lengths, missing NI-write physical times and setup differences
remain limitations. Current-p dependence means inverse pressure computation
requires solving that dependence/common-pressure allocation. No inverse/PAC/PDFB
integration, force/passivity/stability qualification, hardware or SLX run.
Weights/config/old-source/split hashes, checkpoints, numerical CSV/MAT/FIG/PDF/PNG
and source snapshots retained. No commit/push. OI-11 tracks further development.

## 2026-11-03 - Requested sine plots with CV/GRU separated from others

Canonical bundle figures/sine_predictors_20261103_v1. Existing SineWave4
frozen forward predictors only; no fitting or architecture selection. User
requested grouping CV+GRU separately from nominal/ridge/kinematic-ridge. Four
FIG/vectorPDF/300dpiPNG sets: each group as one-period free predictions and
complete-record one-step predictions. Rows P1/P2/P3, columns l2,l3,phi.
Sine is a pressure excitation, not a desired phi trajectory.

Actual cycle-label spacing defines ~31.425s simulation period. First mature
continuous one-period windows: P1 trial1 40.8-72.25s, P2 trial3 214.6-246.05s,
P3 trial6 527.35-558.8s, physical31.88/31.86/31.85s. Same strict gap guards;
initial measured10frames only, predicted subsequent q/v and recorded inputs,
no true-state resets. Different channel start phases/repeats limit comparisons.
All15 model/channel windows complete numerical guard, not physical stability.

Path q2/q3 vector RMSE mm: CV9.604/16.748/20.685; GRU.951/3.899/2.401;
nominal7.018/10.551/9.433; nominal+ridge3.376/2.724/4.714;
kinematic+ridge2.954/3.092/4.851 (P1/P2/P3). GRU is stronger than CV, but has
spikes and both ridges outperform it for P2 here. Grouping does not establish
universal ranking. Predicted nominal/GRU chatter is not measured arm behavior;
no new physical-cause diagnosis made. Same legacy local-affine context and
approximate pressure/physical-time alignment, not full nonlinear ODE testing.

q/phi roundtrip, first-step agreement with existing frozen predictions
(max2.52e-10m across15pairs), independent CSV path score recomputation and
raw/export/normalizer/weight/source hash checks pass. Complete source/traces/
scores/protocol retained. No training/force-model/SLX/hardware/commit/push.

## 2026-11-03 - Requested phi1 circle plots for five existing predictors

User asks to SEE circle behavior after clarifying predictor meanings. Used only
the five listed models, no added architecture/fit or controller test. Logged
phi_reference=1 selects trial3 of November2 circleSweep,259-359s. Measured phi
.05237-.23015rad across five circles; requested1 does not mean actual1.
First full20s circle (20.23physical seconds), common observed initial10frames,
then predicted q/filteredv only and recorded pressure/time sequence. Path vector
length RMSE mm: nominal11.489, constantvelocity2.253, nominal+ridge3.825,
kinematic+ridge4.750, nominal+GRU1.096. All finish numerical guard. GRU closest
on this onecircle in length/PCCXY; constantvelocity phiMAE smaller but lacks
rotation. Long20s ridge drift limits conclusions from previous5s averages.
No universal architecture choice or control tracking success inferred.

PCC XY reconstructed consistently from measured inferred and predicted lengths,
not independent NDI XYZ. Angles/XY show desired reference separately. One-step
100s plots use measured history eachsample and retain GRU excursion. Fullcircle
plots have no state resets. Report/source/CSV/MAT/hashes/FIG/vectorPDF/PNG under
figures/circle_phi1_predictors_20261103_v2; final one-step variant has _labeled
suffix. No training/newpredictors/SLX/hardware/commit/push changes.

## 2026-11-02 - Frozen model evaluation on five fresh experiments

Canonical report: residual_dynamics_learning/validations/fresh_20261102_v2/REPORT.md.
All five supplied exports exactly reproduce observations, commands and model
predictions from their current raw MAT files. Original normalizer, model banks,
ablation weights and inference sources hash-verified unchanged during evaluation;
50,078 strict guarded anchors. All November2 data designated fresh_test, no fit.

Mean per-run independent-length vector RMSE (mm): one-step nominal 4.9183,
constant velocity .06235, ridge .05976, MLP .41796, GRU .34447, kinematic ridge
.05961. Ridge's small one-step improvement is not universal (sine worse).
Step loading/unloading and circle/helix phases show modest improvements;
static/relaxation intervals can favor constant velocity. Per-phase/trial scores
are retained rather than inferring success from a pooled result.

Autonomous transition-model rollouts use only observed initial motion history,
predicted subsequent q/filtered v, known recorded pressure/time sequences and
unchanged gap/warmup policy. Twelve deterministic windows/run/horizon; 60 per
horizon/model; .5/2/5 physical seconds, no measured-state resets. At 2/5 s path
RMSE: constant velocity .8930/3.4076, ridge .7268/1.6204, GRU 1.1265/1.5112,
kinematic ridge .6795/1.8064 mm. No numerical-guard failures among 1,440
model/window propagations. This is not safety, full-circle or controller proof.
No new winner selected from fresh results; historical selection remains ridge.

Old mode0 nominal context explicitly retained for learned weights; corrected
mode1 nominal independently scored and nearly identical on measured states.
Locally frozen-coefficient transitions are not full nonlinear ODE validation.
NI-row pressure ZOH and approximate physical action timing remain limitations.
Raw model versions 1.46 step/combined versus 1.47 sine/circle/helix retained.
Unrecorded calibration/load/setup equivalence is not established.

Conditional inverse full ridge .000559 bar versus pressure-only AR .000555
mean three-pressure vector RMSE again shows command-history shortcut; observed
preview and prescribed common pressure do not constitute a deployed inverse.
No learned pressure applied; MLP/GRU pressure violations remain material.

Conclusion: useful multi-second history prediction, mixed evidence for the
mechanical prior, and no trained generalized-force residual yet. Next implement
versioned corrected-nominal plus learned force augmentation with integrated-q
output-error fitting. Preserve these predictors/baselines. Fresh cohort is now
evaluated, not untouched confirmation for later revisions. Retained numerical
artifacts, hashes, FIG/vectorPDF/300dpiPNG and source. V1 incomplete evaluator
attempt marked DO_NOT_USE. No training/SLX/hardware/commit/push changes.

## 2026-11-02 - SineWave4 acquisition review with raw status

Same-BX diagnostics present/aligned,10/12721 one-frame losses allchannel3:
rawtipstatus2=MISSING, base1, volume/partial/brokenflags0, system0, freshpaired
frames andoverrun0atloss. Mask4098tip=2+4096. Replylength87->55matches32byte
omitted tiptransform. Not a quaternion/angle rejection of available pose;
underlying physical cause unproven, CRCunchecked. Precedingphi.648-.884rad,
4loading/6unloading losses. Separate7overruns/5hostgaps>.10s and12straight
theta-unobservablevalid samples retained. Use unchanged gap/window guards,
not invalid filling. Corrected forward-bound behavior observed in modelphi
peaks2.5268radall3; otherinlinecodeinstallationunverified. No learned predictor
scored, fit, runtime/SLX edit or hardware/model run. Retained report/events/
figures/source under audits/validity_sine4_20261102_v1.

## 2026-11-02 - NDI diagnostics implemented for the next acquisition

User authorized NDI output extension. See NDI_DIAGNOSTICS_RESULTS.md and
ndi_diagnostics/README.md. Passive same-BX32x2 output (BaseRotation23,
Position_only20), raw tracker/port words, quaternion/frame/host rejection
conditions and packet metadata. Existing acceptance/holds/pressures unchanged;
checksum not verified.32 synthetic replay cases plus named export/decoder/
legacy real-log checks pass. No hardware/Simulink/reader setup or SLX edits.
Operator logging/new record required to diagnose actual losses; old MATs lack
these fields. No training/relabeling. Proposal-only status superseded for
source, not for unverified physical integration.

## 2026-11-02 - SineWave3 recurring tip losses; richer NDI output recommended

Raw20261102/SineWave3.mat reviewed offline after operator generator-move report.
14losses/12721samples (13ch3,1ch1), same tip-invalid/base-valid/NaNXYZ with fresh
paired frames andoverrun0. All14 lastone.05s sensor interval;78nativefalse rows.
Preceding valid ch3phi spans.186-2.202rad in several groups; exact missing pose
unknown, no interpolation.11/13 have other valid samples>1s away within3mm of
preceding validXYZ,21-116examples. This does not prove a unique blind location.
Raw source status bits absent; cannot identify volume/interference/cable or host
rejection uniquely.15/20/14losses across3runs are not controlled causal evidence.

User asks whether to add NDI diagnostics. Recommend same-BX raw sensor status
and complete portword, volume/partial/broken flags, quaternion norm and host
rejection mask, with existing sensor/frame/time/error signals. Preserve existing
ports/validity/hold/pressure behavior, avoid extra device queries/per-step disk
I/O. Current parser discards packet/system status; do not claim CRC validation.
Diagnostic runtime implementation remains proposed, not authorized/performed by
this question. No reader/driver/SLX/hardware changes. Full report, event/pose
CSV, FIG/PDF/PNG and source hashes in audits/validity_sine3_20261102_v1. Final
canonical pose plot is a2D projection; failed3D export attempt preserved.
No model fit/learned prediction evaluation; logged nominal stilllegacy.

## 2026-11-02 - Standing reminder for other PAC and PDI/FB controllers

User explicitly asks to retain the three-actuator bound correction for future
PAC and PDI/FB controller work. Recorded permanent decision D-12 and pending
migration issue OI-10; updated current state/architecture/installation guidance.
When the same reduced coordinates are used, align forward/inverse/controller
force paths and audit embedded copies. Keep35mm, include l1; distinguish actual
model versions and preserve old weights/results. No claim of all-model migration.
This turn changes documentation only, not controller code, SLX or hardware.

## 2026-11-02 - Equal actuator bound correction implemented in source

User authorized correction; limit already+-35mm, no increase. Shared
armS_single_boundForce reconstructsl1=-l2-l3 and maps identical actuator force
laws through J'F. Current RHS defaultsboundMode1; PAC fixed/adaptive force and
common-pressure inverse source share this helper. Added paste-ready collection
rlInverseModelBlock.m, preserving3inputs/5outputs/invalid hold/allocation. No
SLX edited: saved collection inline inverse must be replaced manually on a fresh
rebuild to match updated externalforward. No hardware/Simulink run or model fit.

Numerical results under residual_dynamics_learning/validations/actuator_bounds_20261102_v1:
potential-gradient relative error5.11e-10; cyclic mapping6.53e-17; legacy RHS
3.28e-16. Static forward/inverse/PAC agreement including implied-l1 bound passes,
both-mode nominal inverse roundtrips and greybox causality pass. Existing fixed/
adaptive PAC numerical and inverse dropout/fault suites pass. Device-free same
3bar ramp/hold corrected equilibria2.526837rad inallthree channels; legacy
3.005137/2.526837/2.526837. Retained trajectories, verification/source/hashes and
FIG/PDF/300dpiPNG. Counterfactual is not actual-arm validation or model learning.

Historical semantics preserved: RHS params.boundMode0, residualNominal/
residualFrozenStep optional5thmode default0, greyboxSimulate optional4thmode
default0, calibration/nonlinear benchmark explicitly0. Legacy two-cohort
preparation records0; screening guard rejects explicitmode1. Against old source,
legacy frozen step error0 and nominal relative difference1.44e-16. Existing
weights/datasets/reports untouched; new corrected labels/context/training need
a versioned model. Concrete structural fix, not more nominal-parameter tuning.

Additional standalone Coder checks failed after C emission during toolchain
initialization/build, including shorter-path retry; saved status/transcript.
Do not claim a passed C library or Simulink build. No compiler installation
settings altered. See RL_data_collection/ACTUATOR_BOUND_CORRECTION.md for manual
installation and limits; no commit/push.

## 2026-11-02 - Simulated phi asymmetry traced to omitted l1 bound force

User reports channel1 computedphi~3rad versus~2.5rad onchannels2/3. Recorded
des_phi confirms3.005162/2.526838/2.526837rad, identical in both20261102 sine
files. Actual orientation-derived peaks are1.952/1.935/2.494rad (SineWave),
1.951/1.932/2.494 (SineWave2), so distinguish model prediction from measurement.
Report/artifacts: residual_dynamics_learning/audits/bound_asymmetry_20261102_v1.

Current saved forward chart sets bounds+-35mm, strength1e6,mu2000. Existing
armS_single_dynamics adds bound forces only toq=[l2;l3], not impliedl1=-l2-l3.
At channel1 model peak,l1=39.067mm,l2=l3=-19.534mm; existing boundforce~zero.
Channel2/3 extending coordinate32.849mm already incurs~6.020N soft restoring
force. No channel-specific phi cap exists. This is an identified asymmetric
model term if intended to represent equal actuator limits; do not interpret
it as measured directional physical stiffness. Proposed three-actuator force
mapping J'F(l)=[F(l2)-F(l1);F(l3)-F(l1)] is not implemented or simulated here.
Correction would need versioned baseline/revalidation; old benchmark results
remain records of the existing asymmetric model. No learned predictor scored,
nominal simulation run, controller/helper/SLX changed or hardware executed.
Retained peak/plot-data CSVs, FIG/PDF/PNG and source provenance.

Operator says channel3 brings tip near generator; repeats underway with generator
moved back. Record as setup change/hypothesis, not proven dropout remedy. Repeat
straight calibration/preflight because current reference/gravity handling assumes
fixed NDI frame; preserve prior runs. New recording not yet supplied/reviewed.

## 2026-11-02 - Actual SineWave/SineWave2 configuration-validity drops

User supplied RL_data_collection/20261102/SineWave2.mat and requested both sine
recordings. Read both original out logs without model execution. New report:
residual_dynamics_learning/audits/validity_20261102_v1/REPORT.md. Each12721 sensor
rows. Invalid configuration intervals15/20 (0.117915%/0.157220%); every one is
one0.05s simulation-time sensor interval. All match tipValid0/baseValid1,
orientationValid0, nonfinite relativeXYZ, fresh paired frames, overrun0 and
tip raw errorNaN. Native false rows80/107 reflect faster/variable solver logs,
not separate sensor failures. Native checks pass with1ns floating-timestamp
tolerance, no physical lag compensation or invalid interpolation.

Losses cluster in channel3 sine excitation:14/15 and19/20, trials5/6; remaining
one each is channel1 at14.05/14.15s. SineWave2 example429.70->429.75s. Therefore
configuration drops originate in unavailable tip acquisition, not a separate
XYZ-IK failure with valid input. Current configuration helper rejects nonfinite
XYZ, no fit-accuracy threshold; chart holds prior output while flagfalse. Actual
NDI status bits are not logged, so do not claim a unique physical cause such as
out-of-volume. No repeated/mismatched frames or flagged overruns at these drops.

Saved event/summary/native CSVs, exporter tables, original-source hashes/current
source snapshots and two FIG/PDF/300dpiPNG diagnostics; see FIGURE_INDEX.md.
Do not relax flag; mask losses and avoid crossing windows. Full quality/model
qualification remains separate. No new data used for fitting or predictive
evaluation, no source/hardware/SLX/controller/commit/push changes. Development
audit failures retained in verification.txt; final recorded checks pass.

## 2026-11-02 - Validity-flag source trace (specific event not yet identified)

User clarified the signal is orientationValid; no recording/event was specified,
so no event-specific physical cause is established. Read current
AuroraReader_s2_BaseRotation, ndiSupportFnc, exportRlTrainingData and preparation.
Raw sensor validity rejects device-invalid/nonfinite poses, nonfinite/repeated
frame IDs, and mismatched base/tip frame IDs. orientationValid requires both
sensor valid bits, available calibration and valid quaternion orientation
(each quaternion finite, four elements, norm within0.02 of1). Lost orientation
frames hold the last theta/phi/length estimates while orientationValid=false;
smooth-looking held outputs are not new valid measurements. overrun is logged
separately and does not directly set orientationValid=false in the reader.
thetaObservable=false near straight is distinct from orientation invalidity.

Exporter dynamics_observation_valid additionally requires fresh paired frames,
finite phi and covered finite three-pressure commands. transition validity
also requires a valid predecessor and acceptable host interval. Prepared S.valid
adds base validity, finite q, command-age/bound and positive trial checks;
physical/simulation gap checks and0.60s warmup govern segments/maturity/windows,
not all the same Boolean. These acquisition/preprocessing flags are not residual
model error or training/validation partition assignments. No flag changed,
threshold relaxed, hardware run or cause attributed without the actual log.

## 2026-11-02 - Proposed fresh-data test of the current frozen predictors

User asks whether current residual models can be tested with fresh repetitions.
Recommendation: yes, first as offline forward prediction on newly collected
data using the established operator-run collection procedure. Start with one
fresh StepWave and one fresh SineWave recording, retaining complete original
trial/channel sequences, existing excitation settings/limits and full out.
Keep payload, mounting, reset/calibration and collection settings documented.
Record both sensor physical timestamps and simulation/command times with all
validity flags; intended sampling does not establish actual physical intervals.
Do not change excitation amplitude/speed just to improve model scores. Additional
repeats and circle/helix/combined coverage can follow the initial quality check.
This is a proposal, not a hardware action or an already accepted run schedule.

Freeze current ridge/MLP/GRU weights, existing training normalizers/filters and
simple/kinematic baselines before new-data evaluation. Score one-step plus
longer autonomous rollouts at matched physical horizons, with gap guards and
failure counts; preserve plots and per-run data. This tests external prediction
and loading/unloading response, not controller improvement or novel desired-
motion inverse command generation. Keep the fresh data out of fitting until
that frozen test/report is finished; subsequent use for fitting changes its role.

Inspected residualPredict/residualBenchmarkRollouts/prepareResidualData: frozen
inference exists, but preparation is pinned to the original ten-file manifest
and rollout selection is day1-validation-only. A separate fresh-cohort evaluation
adapter is still needed; do not append new files to the training manifest or
refit normalization. Current forward residual corrects position transitions,
not force; no direct injection into the PAC/PDFB ODE or learned-pressure hardware
test is proposed. No hardware/controller/SLX changes, fitting or commits here.

## 2026-11-02 - Scope clarification and paper-artifact retention

User challenged the value of further nominal-model calibration because model
mismatch was already established. Clarified its limited purpose: ask whether
bounded parameter corrections/simple actuator memory explain enough mismatch
to simplify learning. Result does not establish a unique physical cause, and
calibration was overemphasized relative to the residual-learning objective.
Do not require an accurate nominal model before learning its missing dynamics.
Freeze an explicit baseline; proceed toward a history-dependent generalized-
force residual fitted/evaluated through integrated length predictions. This
is the next intended implementation, not an already fitted force-residual model.
Prior step3 outputs were next-position residuals. Test longer autonomous held-out
rollouts and loading/unloading with nominal, simple and data-driven baselines;
then consider shared model augmentation for PDFB/PAC. Timing/quality findings
remain applicable but are not a reason for indefinite calibration work.

User also requests retaining potentially useful paper plots. Added
project_docs/FIGURE_INDEX.md and a versioned model_diagnostics_v1 bundle under
residual_dynamics_learning/figures/20261102_model_diagnostics_v1: four saved-table
summary plots with editable FIG/vector PDF/300dpi PNG; original example rollout
PNG preserved; numerical data/source reports/code/captions/provenance retained.
Figures are diagnostic paper candidates, not evidence of residual/control success.
No legacy raster is represented as vector data. Earlier PAC plot files remain
in their original folders; this is not a complete retrospective PAC figure audit.
Local artifact retention is separate from committing. No training/hardware/SLX,
manuscript, commit or push changes for this documentation/export work.

## 2026-11-02 - Static/dynamic mismatch, constrained mechanics and effective memory

Executed offline using only October19 prepared data. Reproducible report and
artifacts: residual_dynamics_learning/calibration/greybox_20261102_v2/REPORT.md.
15 StepWave training holds fit K/bias;10 active training windows fit D and then
shared first-order command memory at delays0/.1/.2s. Mass/geometry/force gain
fixed. 30 whole-trial-disjoint validation windows, six per waveform, full existing
nonlinear dynamics integrated without future-state resets. Native command clock
mapped approximately to sensor physical time; no measured chamber pressure or
NI-write wall time. October27 not loaded in this calibration.

Mean ~2s vector pathRMSEmm: original5.13569,mechanical3.36946,
memorydelay0=3.32361,delay.1=3.31408,delay.2=3.29525; no domain failures.
Matched persistence1.33761 and constantvelocity1.05490 are substantially better.
CV also wins at~.5/1s. Do not interpret34% improvement over original as successful
missing-dynamics identification or control qualification. All are prediction
errors, not tracking length-RMSE results; windows are not independent repeats.

Fitted k=[1349.866;630.474;1008.796]N/m, b=[.810708;1.821780]N,
D=[250;250]Ns/m BOTH at upperbound. Training-only shared-D40/100/250/500/1000
profile monotonically lowers cost5.3833/4.9337/3.9428/2.7434/1.5689; no refit or
expanded-bound validation performed. Excess damping may suppress inaccurate
predicted motion; these are effective coefficients, not physical estimates or
PAC recommendations. Selected memory delay.2s is grid edge, T.01066s nearlower
bound/subsample; delays0/.1 giveT.11242/.04330s with nearly flat training cost.
No uniquely identified actuator delay/timeconstant; memory adds only~2% beyond
mechanics. Bounds were not widened based on validation.

Four static validation holds: original mean vector equilibrium error7.604mm ->
5.583mm. Nineteen diagnostic training holds9.343->5.638mm. Fit used15loading
holds, not independently matched unloading holds. Significant equilibrium error
cannot be solved by pure input lag alone. Candidate causes include effective
force map/common-pressure dependence, history, calibration/geometry and reduced
coordinates; no cause uniquely established. Next inspect those residuals and
timing, then design targeted operator data if coverage insufficient before
flexible residual fitting. No arm experiment queued, controller/SLX changed,
hardware run, commit or push. Test passes: original-RHS equivalence, future-q
invariance, delayed-input causality, fit/validation trial isolation. Exact source
snapshots/results saved; incomplete first calibration attempt marked INCOMPLETE.

## 2026-11-01 - Step3 offline screening, trained predictors and ablations

Canonical report residual_dynamics_learning/benchmarks/screen_20261101_v2.
Separate forward/inverse ridge, MLP32/16tanh and GRU16 trained on11,785Oct19
anchors(stride3), validation13,977; seed73, CPU, <=20epochs, patience4. GRU
history10frames resets per window. Not a broad architecture/seed/memory search.
Transforms/main validation selection frozen before external evaluation24,079
Oct27anchors; no test fitting. V1 repeated with doubleQR after single-precision
normal-equation warning, same protocol/grid. Exact QR reproduction error0,
lambda1e-6 both tasks. Weights/training histories/source retained.

Forward output is next POSITION residual around a frozen local nominal step,
not acceleration/force. Macro validation vector RMSEmm: nominal4.98454,
persistence.10352,CV.06128,ridge.05802,MLP.41469,GRU.30156. External nominal
5.42676,CV.05193,ridge.08939,MLP2.47972,GRU2.80554. Kinematic-only ridge
val.05807/ext.04887; nominal prior benefit not established by these results.
Tiny external step/combined subsets and timing shift limit generalization claims.

Autonomous validation100rollouts/model,10steps~.5s, observed initial history
only, known sampled inputs/times, no true-state resets. Macro final RMSEmm:
frozen nominal8.20,CV.37,ridge.33,MLP3.10,GRU1.59; no.05mguard failures.
Actual nonlinear armS_single_dynamics separately integrated viaode15s with
integrated velocity, same100windows: path6.6868/final6.9045mm, no guard failures.
This does not establish full-circle or control performance.

Conditional inverse pressure-command macro RMSEbar validation/external:
nominal.50588/.46741,persistence.00985/.00906,ridge.000614/.000626,
MLP.03613/.02779,GRU.03294/.02764. Validation-motivated shortcut ablations,
separately frozen: pressure-history+prescribed-common AR.000621/.000424,
motion/preview withoutpressurehistory ridge.32888/.40082. Tiny pressure errors
therefore do not establish missing mechanics or novel-motion inverse capability.
Ridge external max bounderror.00166bar, none>.01; MLP/GRU max.07308/.08071bar,
19.86%/28%examples violate>.01. Outputs not clipped to hide errors.

validateResidualScreening passes train-only transforms, target-free inference,
exact QR reproduction and constraints. Nominal/feature tests pass. Oct27 has
now been evaluated, including numerical repeat/ablations: fresh independent
confirmation is needed for later claims. No trained model qualified for arm
testing. No SLX/controller/hardware/commit/push changes. Next refine prior/
actuator-memory and motion-conditioned identification/excitation, rather than
simply increasing neural size. Full protocol, limitations and per-run tables
in benchmark report.

## 2026-11-01 - Step2 preparation and integrity validation

Prepared all ten raw recordings in a new versioned output. Valid output:
residual_dynamics_learning/datasets/two_day_v2. Preserved day IDs, hashes,
original native command streams and sensor/control clocks. Bounded previous
command alignment, physical-gap ceiling.10s, simulation gap.075s, command
age.075s, trial/channel/split/invalid segmentation; no invalid filling.
Causal derived velocitytau.20s and accelerationtau.30s, physical-time warmup.60s.
Window history10frames and inverse preview5frames, same-segment mature guards.
October19 whole-trial development:35,350train/13,977validation anchors;
October27 held-out24,079. Total73,406 matched forward/inverse candidate anchors.
Overlapping windows are not statistically independent trials.

Data pipeline first attempt transcribed bound-penalty signs incorrectly and
excluded all inverse anchors, failing completion. It was quarantined under
two_day_v1/DO_NOT_USE.md, no raw changes. Corrected helper verified directly
against existing armS_single_dynamics inside/near/outside bounds and feasible
inverse/forward roundtrip. V2 labels are unclipped nominal differential-pressure
residual candidates and acceleration candidates from DERIVED acceleration;
prefer measured-q rollout/output-error checks before claims of model quality.
Nominal infeasibility retained as flag, not filtered away to hide model mismatch.

All logged effective-pressure force maps verified. Historical inverse
pressure-force convention: October19 step/sine/combined require.8bar activation
to match logged tauAvailable; seven remaining runs match0. Errors<1e-6N after
documenting per-run convention, zero-activation mismatch10.61858N for those
three. This is model-code history, not physical-pressure identification.
Prepared labels are recomputed against one declared zero-activation baseline,
not heterogeneous logged inverse predictions. Source snapshot/dependencies saved.

Timing quality: day1 median physical dt~.0501s. October27 sine~.0709s,
step~.0885s, combined~.0877s, circle/helix~.0502s. Long gaps>.10s:
step2001,combined1194; complete guarded windows only46/26 respectively. These
remaining subsets are not representative of those whole experiments. External
sine2952,circle11339,helix9716. Policy unchanged; do not claim matched physical
speed from identical simulated waveform names. Actual NI-write wall times not
logged, so precise physical input delay remains unverified.

validateResidualData checks every window's same segment/split and whole-trial
isolation, recomputes training-only means, tests changing current pressure cannot
alter inverse history and future measured q cannot alter forward inputs.
All pass; normalizers never fit to external day. No training, architecture
selection, test predictive evaluation, SLX edits, hardware or commits performed.
Next step3 is offline nominal/memory/windowed-MLP/GRU benchmarking, respecting
locked test split and candidate-label/timing limits. See dataset README.

## 2026-11-01 - Step1 interface specification and two-day audit

Read all ten out.logsout recordings using auditResidualInterfaces, without
running/opening hardware models. Output audits/interfaces_20261101 under
residual_dynamics_learning contains runs.csv, signals.csv, physical command
candidates, interface_audit.mat, source hashes/provenance and saved SLX scripts.
All ten hashes distinct; five runs per acquisition day confirmed from stored
TimingInfo:2026-10-19 root and2026-10-27 subfolder. Model names RL_dataCollection
andRL_dataCollection_v2 respectively. No cross-day merge or assumed identical
calibration/load. Each has StepWave,SineWave,combinePressure,circleSweep,helix.

Physical actions: final three NI pressure command proxies under Pneumatic Valves;
not requested/nominal predicted/effective pressures. State response q from
orientation-inferred length channels2/3, dq derived. Nominal dynamics found in
armS_single_entry/dynamics/core and extracted collection-model scripts. Saved
inverse differs on invalid samples(v1zero versusv2hold/invalid flag). Historical
logged identity does not prove current SLX snapshot matches historical bytes.

Documented forward residual acceleration2-vector and conditional inverse
three-pressure residual. Nominal A=[-1,1,0;-1,0,1] has common-pressure nullspace.
Exact rational arithmetic verifies specified B: A*B=I2 and common mean of B*d=0.
Thus initial inverse design recommends externally specified common pressure
plus two differential residuals, reconstructing all three commands. Do not
label allocation ambiguity or hardware saturation as uniquely mechanical error.

Current benchmark decision: nominal/zero, simple input-memory model, delayed-MLP,
smallGRU; LSTM only if justified. No trained winner or performance claim.
Targets/preview availability and day-aware held-out splits in MODEL_SPECIFICATION.
Training-ready windows/labels are step2, not produced here. No extra hardware
data collection, control changes, commits or pushes made in this step.
