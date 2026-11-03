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

## 2026-11-03 - First force-residual model training and held-out evaluation

Canonical residual_dynamics_learning/force_model/release_20261103_v1/REPORT.md.
New force GRU32 trained October19+27 only, old-only whole-trial selection, new
normalizer older TRAIN only. Selected epoch38/lr.0003; frozen checksum before
November2 evaluation. Force RHS sign/units equivalence, existing nominal
coefficient agreement1.56e-15, direct and network-through-integrator gradients,
target-free propagation and old-only split/hash checks pass. Four-substep SDIRK
ODE fixture errors<.034mm; doubled-substep older path differences<.077mm.
No noisy acceleration/force target used. Initial failed numeric/train attempts
are marked DO_NOT_USE; completed candidates/checkpoints preserved.

50,078 new-data anchors,540 matched rollout paths, all numerical guards pass.
Macro q2/q3 vector RMSE mm, nominal/CV/forceGRU:
one-step4.7352/.06235/2.1510;2s7.0438/.8927/1.8375;
5s6.6372/3.4062/2.1267. Improvements versus nominal and CV at5s coexist with
weak short-horizon force-model accuracy; no controller/inverse qualification.
Held-out data previously inspected, no post-evaluation refit. Source/weights/
CSV/MAT/FIG/vectorPDF/PNG retained. No Simulink/serial/DAQ/hardware execution.

## 2026-11-02 - Frozen residual models on fresh exported recordings

Five raw MAT re-exports exactly match supplied observations/commands/model
predictions. Frozen model, normalizer and source hashes unchanged during test;
50,078 gap/whole-trial/mature-history anchors, forward target-free and inverse
current-pressure exclusion checks pass. No training or hardware/Simulink runs.
One-step mean per-run q2/q3 vector RMSE: nominal 4.9183 mm, constant velocity
.06235, ridge .05976, kinematic ridge .05961. At 2/5 s path RMSE, ridge
.7268/1.6204 versus constant velocity .8930/3.4076 mm. GRU 1.1265/1.5112;
kinematic ridge .6795/1.8064. All 1,440 offline model/window propagations avoid
the numerical-domain guard; not a stability/physical-safety claim. Per-phase/
trial/run metrics, failure counts and deterministic windows retained.
Physics-prior value remains mixed; inverse pressure-history shortcut persists.
Versioned report/source/CSV/MAT/FIG/PDF/PNG: residual_dynamics_learning/
validations/fresh_20261102_v2. Incomplete v1 marked DO_NOT_USE. November2 now
evaluated; no force-residual training or controller qualification claimed.

## 2026-11-02 - Passive same-BX NDI diagnostics

32 synthetic reader-body cases passed with actual decoder/cache: legacy
non-wall-clock output equality, status availability, quaternion/frame/calibration
reasons, straight/wrap cases, one poll, incomplete exception/recovery and output
metadata. Named32x2 export/decoder pass; actual old SineWave3 observations
unchanged. No serial/DAQ open, reader setup/preflight or Simulink model run.
Physical timing/port integration and CRC integrity unverified. Source/fixtures/
transcript: ndi_diagnostics/validation_20261102_v1; wiring: ndi_diagnostics/README.md.

## 2026-11-02 - Three-actuator nonlinear bound-force correction

Soft bound remains+-35mm, strength1e6,mu2000; corrected omittedl1 contribution
with shared J'F mapping. Forward, inverse and PAC source consistency checked
offline, no SLX edits/hardware/Simulink execution. Force potential-gradient
relative error5.11e-10, cyclic virtual-work symmetry6.53e-17, historical-mode
RHS3.28e-16. Existing fixed/adaptive PAC and inverse fault checks pass.
Same3bar synthetic ramp/hold equilibria: legacy3.0051/2.5268/2.5268rad ->
corrected2.5268rad allthree; not experimental tracking/model-accuracy evidence.
Historical screening transition matches archived source exactly; mode0 retains
old dataset/model semantics, current defaultmode1. Coder toolchain step failed
after C emission; successful compiled integration not established. Collection
embedded inverse requires manual paste of rlInverseModelBlock.m before rebuild.
Full results/source/plots: residual_dynamics_learning/validations/actuator_bounds_20261102_v1.
Instructions: RL_data_collection/ACTUATOR_BOUND_CORRECTION.md. No commit/push.

## 2026-10-28 - Final smooth candidate10; campaign budget exhausted

Imported smooth_candidate10_00g_phi1_T10_12cycles_r01.mat as smooth_candidate10_r01.
Exact original-filter replay max error0; full tuning/software initialization match,
no live edits/protocol mismatch; physical reset per standing user attestation.
valid_infeasible: late110-140 l2/l3 RMSE[2.880204,3.166379]mm,vector4.280366mm,
velocity.003572203m/s,roughness.016077616m/s^2. Both2mm limits fail; other caps
pass. Coverage99.5%,pair99%; tracking coverage99.083%,adaptation98.167%,tracking
faults/pressurelimits0. Lastlap length[2.6794,2.9826]mm still above targets.
Smoothest measured,~16.7% less rough than mean baseline, but unacceptable length
tradeoff under current criteria. No claim of isolated gamma causality/convergence.

Adapter returns stop/evaluation_budget_exhausted,0remaining,success=false. No new
candidate generated. Best qualified screened reference remains original baseline
gamma.5/.002/5,roughness.01937096/.01922323 with length[1.5684,1.8887]/[1.5290,
1.8569]mm. Feasible03/04/07/08 are rougher;08 has best length accuracy. Wider
exploratory09/10 did not improve the constrained smoothness objective.
Historical data were used in training,not held-out confirmation; future repeats
or search require budget review. Saved final_budget_review,completed_search_review
CSV/MAT/tradeoffplot and reproducible reviewSmoothBudget.m. No raw/model/controller
changes,hardware execution,live gains,commit or push by assistant.

## 2026-10-28 - Smooth candidate09 rejected; final search point queued

User supplied smooth_candidate09_00g_phi1_T10_12cycles_r01.mat. Verified original
measured-filter restoration by exact recorded-input replay,max error0. Import
smooth_candidate09_r01 matches full tuning/software initial state, no live edits
or protocol/manifest mismatch; physical reset per standing user attestation.
valid_infeasible: last3circles lengths[2.527911,2.752004]mm fail2mm, vector3.736825mm;
velocity.003515224m/s, fixed-raw roughness.022076340m/s^2 pass caps. Coverage99%,
pair98%; tracking coverage99%,adaptation98%,faults/pressurelimits0. Lastlap
length[2.3255,2.5333]mm still fails; not smoother than baseline~.0193 and not
an improvement in accuracy. No causal claim about any individual gamma.

Queued smooth_candidate10 from feasibility-weighted EI: k.25306817650802854,
d.00069737816290234929,h2.6573766973344188. Normalized prediction mean.8121412,
std.1511505 corresponds to roughness~.02030/.00378m/s^2; mean not better than
baseline, acquisition includes uncertainty. Nine observations/one slot left;
after next import stop/review budget and constraint-margin modeling. No automatic
extra runs or formal confirmation from historical screenings. Artifacts and
status_after_smooth09/status_after_queue_smooth10 saved. No controller/model
edits, live gain write, hardware execution, commit or push by assistant.

## 2026-10-28 - Broader smoothness-primary campaign prepared

User requested increasing the proposed search ranges and testing performance.
prepareSmoothRangeCampaign creates new immutable campaign_smooth_v2 with
gammaK[.1,1],gammaD[.0005,.008],gammaH[1,10],original measured filter.20/.30s,
2mm coordinate ceilings and fixed-raw-length roughness objective (normalizer.025).
Velocity.01m/s,roughness.025m/s^2,zero limiting and95% coverage retained. Eight
prior old-filter runs verified by exact measured-velocity replay before reuse;
changed-filter trial excluded. All are historical training,not held-out validation.
Source snapshots archived; old campaign/records untouched. Standing user reset
attestation preserved, no independent physical calibration verification claimed.

New config optimizer.continueSearchAfterFeasible=true explicitly continues
objective improvement after a feasible screening; default false preserves prior
confirmation-first behavior. Unit/integration tests cover continued search,
queueing and pending-observation status. Manual/evaluator/optimizer/adapter
regressions pass. This scheduling flag does not change GP fitting equations.

Rescoring makes baseline repeats feasible with best mean roughness (~.0193);
candidate03/04/07/08 also feasible,05/06 infeasible. New GP suggestion uses
feasibility-weighted EI with a feasible incumbent: smooth_candidate09 gammaK
.90453500157774203,gammaD .0022428161211506688,gammaH2.8192356276278843.
Predicted normalized mean.7256975/std.1169877 (~.0181424/.0029247m/s^2), heuristic
feasibility.584496 is not a calibrated probability or physical guarantee.
Budget retained at10total observations:8reused+2new search slots, confirmation
requires review beyond that. Original filter must be restored before the next
user-run trial. Saved initial_status and RUN_NEXT.md. No controller/model edit,
live gain write, hardware execution, commit or push by assistant.

## 2026-10-28 - Actual-arm replacement filter worsens roughness; restore original

User submitted filter_smoke_00g_phi1_T10_12cycles_r01.mat. It is a12circle run
with candidate08 gamma(.75,.008,7.5),not suggested baseline. Compared against
candidate08 using identical fixed raw-length analysis over110-140s. Exact new
filter velocity replay error0; old-source replay max difference.008204m/s,
confirming the tested implementation. No live tuning changes, tracking faults
or pressure limits; new tracking coverage99.083%,scoring99.5%,pair99%.

Old/new qRMSEmm [1.172517,1.302891]/[1.189358,1.416503], vector1.752804/1.849609;
velocity RMS.00300670/.00357242m/s; roughness.02396971/.03408892m/s^2 (+42.2%,
new exceeds.025). Secondary theta-rate SD.39579/.52262rad/s,near-stop8.08/10.10%,
phiPeak1.1673/1.2111. Supports reverting to original measured.20/.30s; added lag
is plausible but one run does not establish causation. Raw inputs filtered offline
had appeared smoother; that did not predict improved closed-loop physical motion.

Prepared same-port pacMeasuredDerivativesOriginal.m from saved original source,
updated smoothness-first draft to select it, and marked experimental source as
not selected. Desired filter remains.10/.15s. Saved actual_arm_comparison reports;
changed-filter data excluded from old optimizer campaign, ledger remains8/10.
No new active campaign/candidate, SLX edit, hardware execution, commit or push.

## 2026-10-28 - Filter replay and smoothness-first preparation

User requests smoother velocity/acceleration filtering, wider search and higher
priority for smoothness, then explicitly approves2mm RMSE ceilings for each
length. Read saved PAC_autotuning_manual.slx XML only: measured filter backward
differences with velocity tau.20s/acceleration.30s; desired.10/.15s. PAC force
uses measured q/dq, not measured ddq. No SLX load/run/edit. Saved measured source
replays baseline1 and candidate08 logged dq exactly with pollStart_s timestamps;
controller-time substitution differs, so clock provenance is material.

Prepared same-port standalone pacSmoothMeasuredDerivatives: two-pole velocity
tau.15s/stage with exact ZOH propagation and acceleration tau.35s, frame/gap
handling, raw length passthrough. Earlier.10/.20s prototype worsened replay
non-fundamental content and was not selected; its summary retained. Selected
candidate replay reduces velocity non-fundamental RMS about7.6/11.0% and
acceleration19.5/23.5% on baseline1/candidate08 respectively, but adds roughly
.07-.08s velocity-estimate harmonic lag. These are same-input filter effects,
not measured reductions in physical jerk or proof filters caused the roughness.
Artifacts under filter_review_v1: saved source, comparison.csv, replayMAT/PNG,
reproduction function and detailed smoke instructions. No raw data changed.

Added configurable objective (old configs remain vector length RMSE) and fixed
raw-length analysis velocity for the new draft. Fixed analysis taus measured.20,
reference.10 stay constant regardless of live filter settings. Regression verifies
that changing only logged velocity/reference-velocity channels cannot game the
fixed-analysis score. Draft minimizes velocity-error-rate RMS, with2mm coordinate
constraints, existing velocity/limiting settings and provisional roughness.025.
Proposed gamma bounds k[.1,1],d[.0005,.008],h[1,10]; damping upper not increased.
No active new campaign, old immutable config/ledger altered or gain applied.

validatePacSmoothMeasuredDerivatives passes ramp/invalid/repeated/gap/synthetic
high-frequency checks. Evaluator/manual/optimizer/adapter regressions pass after
objective/scoring changes. Next physical step is user-run filter-only smoke at
baseline gamma.5/.002/5 then12circle assessment if acceptable. Physical filtering
performance and2s-speed behavior remain unvalidated. No hardware execution,
controller/SLX edit, commit or push by assistant.

## 2026-10-28 - Candidate08 best observed; acquisition/budget review

Imported candidate08_00g_phi1_T10_12cycles_r01.mat as candidate08_r01: exact
tuning/software initialization, no live edits or protocol/manifest mismatch,
physical reset per standing user attestation. Valid_infeasible on both1mm length
limits. Late110-140 qRMSE[1.172517,1.302891]mm/vector1.752804mm, velocity
.002993780m/s, roughness.023912515m/s^2, coverage99.5%, pair99%. Velocity/
pooled smoothness/zero-limiting pass. Tracking coverage99%, adaptation97.333%,
faults/limits0. Early k2/k3 lower occupancy .958/4.083%; k3ratecap .942% active
transitions. No upper bounds. Vector15.4% below07 and27.9% below mean baseline,
but still not an accepted/confirmed result. Lap12 roughness.026547 exceeds .025
individually while the agreed pooled3lap score passes; lastlap length[1.1433,
1.2684]mm still above1mm. No retrospective acceptance-rule change.

Eight observations/two budget slots remain, no candidate qualifies for2fresh
confirmations yet. Raw adapter suggestion seek_feasible_region saved without
queueing candidate09: k.26075453527454223,d.0010338149762252346,h3.9294515860178212,
predicted objective2.69475mm/std.63161mm, heuristic feasibility.49796. Prediction
does not indicate improvement over08. Identified current acquisition limitation:
all observed labels infeasible, so it favors uncertain feasible-region exploration
without learning individual constraint margins. Recommend reviewing margin-aware
selection and budget before further hardware trials, using a new version if the
algorithm changes. No algorithm or immutable configuration changed in this turn.
Saved REVIEW_AFTER_CANDIDATE08.md/status_after_candidate08 and trial artifacts.
No new proposal event, live gain write, model/hardware execution, commit/push.

## 2026-10-28 - Candidate07 best observed; first fitted-surrogate proposal

Imported candidate07_00g_phi1_T10_12cycles_r01.mat as candidate07_r01, exact
tuning/software initialization, no live edits or protocol/manifest mismatch;
physical reset per user attestation. Valid_infeasible: late110-140s qRMSE
[1.236682,1.662630]mm/vector2.072130mm, velocity.002995711m/s, roughness
.022585094m/s^2. Length targets fail; other configured limits pass. Coverage98.5%,
pair97%; tracking coverage98.625%, active97.042%, tracking faults/limits0 and
rate caps0, no upper bounds/stiffness lower contacts. Final-lap qRMSE[1.2090,
1.5526]mm. Best observed vector score6.3% below03; not repeated confirmation.

Seven objective observations/6 unique points complete initial design. Adapter
queued candidate08 reason seek_feasible_region: gammaK .75,gammaD
.0079999999999999967,gammaH7.4999999999999991. This is upper-bound exploration,
not a predicted improvement: GP mean2.56631mm/std.60053mm and heuristic weighted
feasibility.49866. No feasible incumbent exists; acquisition favors uncertain
regions, and binary labels do not measure individual constraint margins. High
gammaD near this box edge previously coincided with05 roughness violation; multiple
gain changes prohibit isolated causal conclusions. Caveat recorded in run guide.
Three observations remain including two confirmations. Saved status_after_candidate07
and RUN_CANDIDATE08.md. No live gain write, hardware/model execution, controller
edit, commit or push. No performance/stability guarantee for this exploratory point.

## 2026-10-28 - Candidate06 reviewed, candidate07 queued

Imported candidate06_00g_phi1_T10_12cycles_r01.mat as candidate06_r01 with
exact tuning/software reset, no live edits or protocol/manifest mismatch;
physical reset/calibration per standing user attestation. Valid_infeasible:
late110-140 qRMSE[2.633751,2.903274]mm/vector3.919903mm, velocity.003443642m/s,
roughness.019350617m/s^2. Both1mm limits fail; velocity/smoothness/zero-limiting
pass. Coverage98.833%, pair97.667%; tracking coverage98.958%, active97.917%,
tracking faults/pressure limits0, no rate caps or upper bound contacts. Small
early k3 lower occupancy.417%. Lap12 length[2.4442,2.7215]mm continues learning
but does not meet fixed-duration objective. Smoother than05, less accurate than03;
no isolated causal conclusion about one gamma. Candidate03 remains best observed.

Queued final initial-design candidate07 gammaK .35438912779204967,
gammaD .0030892129738313034,gammaH6.8259366991000565. Six observations/four
remaining including confirmations. Saved status_after_candidate06 and
RUN_CANDIDATE07.md. No live gain write, hardware/model execution, controller
edit, commit or push. Next acquisition may fit the GP after all6 design locations.

## 2026-10-28 - Candidate05 perceived response-speed follow-up

User felt05 reached the target faster. First-lap l2 RMSE05 3.2578 vs03 3.9665mm
supports a coordinate-specific initial improvement; l3 worsens6.2193 vs5.6682mm,
combined7.0209 vs6.9182mm. A descriptive10s-frequency sin/cos/DC least-squares
fit of desired/measured lengths over valid samples gives early20-50s lag03
[.67360,.84411]s vs05 [.66912,.93106]s; late110-140s03 [.18533,.22565]s vs05
[.29203,.32149]s. Early fits include changing adaptation and are approximations,
not actuator identification or proof of transient crossing time. Data do not
establish overall faster response for05 despite initial l2 improvement. Current
last3circle objective measures sustained tracking, not first arrival/crossing.
No controller/campaign mutation, hardware run or new proposal in this follow-up.

## 2026-10-28 - Candidate05 rejected on roughness and length; candidate06 queued

User supplied candidate05_00g_phi1_T10_12cycles_r01.mat. Exact proposed/full
tuning, software initial states and10s/phi1 reference match; no live edits or
protocol/manifest mismatch. Physical reset/calibration per standing attestation.
Imported candidate05_r01 valid_infeasible. Late110-140 qRMSE[1.836021,2.068622]mm,
vector2.765894mm, velocity.003724713m/s, roughness.033585356m/s^2. Both length
limits and .025roughness limit fail; velocity/zero pressure-limiting limits pass.
Scoring coverage99.333%, pair98.667%; tracking coverage99.167%, active98.333%,
tracking faults/limits0, no stiffness/damping bounds and no parameter rate caps.
Roughness lap7 .027405,lap10 .033727,lap12 .037364; learning-length error falls
while roughness worsens. Near-stop10.30%,phiPeak1.196 are secondary diagnostics.
Candidate03 remains best measured; candidate05 vector error25.1% higher. Three
gamma changes and different-day run prohibit isolated causal claims about gammaD.

Queued candidate06 initial_log_space_design: gammaK .42192531104075343,
gammaD .0020323562614174683, gammaH2.5725609289380671. All other settings fixed;
five observations/five remaining including confirmations; success=false. Saved
status_after_candidate05, per-trial report/artifacts and RUN_CANDIDATE06.md.
No live gain application, hardware/model execution, controller edits, commit/push.

## 2026-10-27 - Candidate04 reviewed, candidate05 queued

User supplied candidate04_00g_phi1_T10_12cycles_r01.mat; imported candidate04_r01
with exact candidate tuning, matching software initialization, no live edits or
protocol/manifest mismatch (physical reset per standing user attestation).
Last3circles110-140s qRMSE[1.692167,1.940342]mm/vector2.574559mm, velocityRMS
.003027893m/s, roughness.020728518m/s^2. Coverage98.833%, pair97.667%, limiting0.
Both1mm length limits fail; velocity/smoothness/limiting pass. Tracking coverage
99.125%, adaptation98.25%, no tracking faults/limits or upper bound contacts.
Rare early k2/k3 lower contacts .708/ .958%, k3 ratecap .594% active transitions.
Late length error16.4% above candidate03: best observed remains candidate03;
simultaneous gamma changes and one run do not isolate individual gain causality.

Next proposal is initial_log_space_design, candidate05 gammaK .28570361215890433,
gammaD .0077863906974842481, gammaH3.5952162270591348. Bounds unchanged and all
other settings fixed. Four observations/six remaining; success=false. Archived
status_after_candidate04 and RUN_CANDIDATE05.md. No live gain application,
hardware/model execution, controller edits, commit or push by assistant.

## 2026-10-27 - Actual candidate03 imported; candidate04 queued

User supplied candidate03_00g_phi1_T10_12cycles_r01.mat (SHA256
93f1c6298709cde0f694639628e975b01b75a3a98d42d2797daa30e27ac9ff71).
Campaign intake verifies exact full tuning, software initial estimates/history,
reference10s/phi1 and no live edits; physical reset/calibration follows standing
user attestation. Imported screening trial candidate03_r01, valid_infeasible.
Last3circles110-140s qRMSE[1.425487,1.690009]mm/vector2.210915mm, velocityRMS
.00314726m/s, roughness.02076540m/s^2, coverage98.833%, pair97.667%. Both1mm
coordinate targets fail; velocity/smoothness/limiting pass. Tracking coverage
98.875%, active97.75%, pressure limiting0, latchedfault0, no upper bound contact.
Early k2/k3 lower contacts .917/1.875%; k3 ratecap .725% active transitions.
Lap12 length[1.3507,1.6040]mm and roughness.022242: continued learning with some
roughness growth; no assertion of convergence or guaranteed stable future motion.
Vector score about9% below baseline mean but one candidate run does not establish
repeatability, adaptation causality or physical parameter identification.

Proposal-only next action is initial_log_space_design. Queued candidate04 with
gammaK .58034718095284454, gammaD .0037139224909066579, gammaH4.4053683649622428.
All other settings fixed. Ledger now3 observations/7 remaining, success=false,
pending candidate04. Saved status_after_candidate03 and RUN_CANDIDATE04.md.
No live tuning write, model/controller edits, hardware execution, commit or push.

## 2026-10-27 - Agreed manual campaign prepared, first proposal queued

User approved k gamma[.25,.75], d[.001,.008], h[2.5,7.5], 12-circle runs scoring
last3, zero limiting, two fresh confirmations and full reset. Trial schedule is
20s startup+120s circles+5s ramp-down+5s pressurized hold, manual pressure kill.
agreedEvaluationConfig now records limitedFraction=0 in addition to approved
length1mm/velocity.01m/s/roughness.025m/s^2 limits. No angular acceptance limits.

prepareLengthCampaign.m created immutable campaign_length_v1, source/model
snapshot and baseline candidate, then imported both existing baseline recordings
with user-attested reset metadata. Software initial states and logged settings
already verified; model snapshot not claimed embedded historical provenance.
Both baselines classify valid_infeasible, not successful, because q RMSE>1mm.
Budget explicitly chosen conservatively as10 total observations including these2
and confirmations, leaving8 new observations; extension requires review.

Proposal-only adapter queued candidate03 from seeded initial log-space design:
gammaK .645028694505006, gammaD .0010090123623767366, gammaH 5.8362570107825071.
Bias gains and every other input remain fixed. Status configurationReady=true,
awaiting_observation, success=false. Exact tuning/provenance saved in candidate
artifact and RUN_NEXT.md. No controller/SLX edit, live gain write, physical run,
commit or push. Single-run batching is not implemented or advised without verified
controller/history/supervisor reset plus consistent physical preconditioning.

## 2026-10-27 - Length-based threshold policy

User requested implementing length-primary acceptance with optional angular
limits. New pac-tune-2 config requires q2Rmse_mm/q3Rmse_mm, vector velocity RMS,
velocityErrorRateRmse_mPerS2 and limitedFraction (0..1); vector length RMSE remains
the objective and its limit is optional, as are theta/phi/theta-rate limits.
pacTuneTargetPolicy centralizes required/active limits in evaluator and manual
campaign readiness. Optional NaN angular measurements/observability do not veto
valid length acceptance unless the angular limit is enabled. Length velocity
valid-pair coverage remains required. Version-1 snapshots keep original angular
requirements; campaign files are not rewritten. Numerical limits are still unset.

Added tests for accepted length-only targets with unavailable angles, active
angular coverage rejection, either-coordinate failure, missing required length
limit, roughness rejection and version-1 policy compatibility. Evaluator/manual/
optimizer/adapter offline regression suites pass. Roughness is RMS of consecutive
velocity-error slopes on valid ungapped samples, not independent acceleration;
filter/sample-period changes require rebaselining. Report tables identify required,
optional disabled/unset and exceeded limits.

Existing actual-arm recordings reevaluated under analysis_length_policy_v2 without
changing raw files: tracking/length results unchanged; late vector velocity RMS
.0030014/.0029658 m/s and velocity-error-rate RMS .019338/.019177 m/s^2, pair
coverage97.33/98.33%. These diagnostics are not chosen acceptance bounds. No
controller/SLX/hardware changes, real optimization run, commit or push performed.

## 2026-10-27 - First actual-arm framework smoke/baseline review

User supplied smoke_00g_phi1_T10_3cycles_r01.mat and two
baseline_00g_phi1_T10_12cycles repeats under experimentations/autotuning/
no_payload_phi1_T10_20261027. Offline reader/evaluator succeeds with all required
signals, verified 3/12 reference turns, 10 s period, phi1, tracking starts20s,
ends50/140s and recorded smooth-return interval to60/150s. Logged configuration
identical and constant: lambda18/16.5, Ks15/10, gamma k.5/d.002/h5/bias1;
mechanical/common enabled, hysteresis freeze false. Pressure activation/min .2,
max3, slew6. This verifies nominal activation configuration, not physical deadzone.
Software initial parameter/history states match; physical reset/trim metadata
remain unverified beyond user statements. No payload per user.

Smoke coverage98.67%, active97.33%, no faults/limits. Baselines coverage99.17/
98.83%, active98.33/97.67%, no faults or logged/replayed pressure limits. Late
110-140 thetaMAE .15202/.14836, phiMAE .06248/.06234, qRMSEmm [1.568,1.889]/
[1.529,1.857]. Final-lap theta .13558/.13221; continued learning on laps10-12
means repeatability at a fixed duration, not convergence. Gate mismatch0;
no final-three-lap parameter bound contacts. Early k2/k3 bounds and k3 rate cap
contacts rare. Pressure roughly .698-2.638/.715-2.637bar, peak rate~1.11bar/s.
Nominal residual~1e-14N is not physical force error. H weights~9-10 by140s,
d~47-48; no late bound bottleneck. Late sampled theta-rate SD .334/.327rad/s,
near-stop5.31/4.75%, phiPeak1.111/1.105: residual roughness/overshoot remains,
without v10-like sustained escalation. K signed/temporarily indefinite does not
identify cause or prove stability; final minimum eigen101.62/111.47N/m.

Reproducible reviewFirstTrials.m and analysis_review_v1 reports save settings,
provenance/hash, windows, gates, lap metrics, evaluation MAT and comparison plot.
All three classify unconfigured_targets, success=false; operational smoke check
passes but no formal optimizer confirmation or real campaign is certified.
Next: define fixed-duration objective, thresholds/bounds and operator metadata.
No raw/model/controller edits, hardware execution, commit or push by assistant.

## 2026-10-27 - Manual workflow and optimizer checkpoint

User requested committing completed milestones B/C. Checkpoint includes the manual
campaign/manifest/ledger, bounded proposal-only GP adapter, four passing offline
test suites, labeled historical/synthetic demonstrations and handoff documentation.
No new implementation changes since the recorded validation runs. Unrelated paper
drafts, backups, duplicate presentation ZIP, RL exports and temporary build/PDF
files remain local. No hardware/model execution or remote push performed.

## 2026-10-27 - Proposal-only optimizer adapter (milestone C)

User requested the next development step. Added bounded fixed-kernel GP search
for grouped gammaK/gammaD/gammaH in normalized log space, with seeded baseline/
initial design, finite-pool expected improvement, explicit feasibility frequency,
repeat/noise handling and budget/distinct-point limits. The core invokes no
objective callback. bayesopt and a Statistics license were detected locally, but
the proposal-only implementation uses base MATLAB numerics/RandStream. New manual
campaigns archive optimizer settings. Existing campaigns without that snapshot
must be versioned before optimizer use; unconfigured historical campaigns block.

pacTuneSuggestNext reconstructs evidence from ledger entries and recommends
propose/observe/confirm/stop/blocked. pacTuneQueueSuggestion records the reviewed
proposal and diagnostics with an event-count/config-hash stale-snapshot check
under the existing writer lock. Confirmation is never optimizer training data.
Valid comparable screening scores train the objective; known comparable latched
faults supply constraint evidence without partial tracking scores. Missing data,
early end without a verified fault and protocol mismatches remain unknown.
Explicit checks prevent live edits hidden by a fault classification from becoming
comparable fault labels. All imports still count against campaign budget.

validatePacTuneOptimizer passes bounds/log mapping, exact baseline, pending-point
avoidance, global RNG isolation, GP/repeat statistics, missing/fault separation,
MAT replay and budget cases. validatePacTuneOptimizerAdapter passes baseline
queueing/provenance, pending observation, stale ledger rejection, held-out
confirmation, fault/protocol distinction and deterministic ledger replay.
validatePacTuneManualCampaign and validatePacTuneEvaluator regressions pass.

reports/optimizer_synthetic_v1 contains 30 synthetic observations with seeded
noise, explicit repeats, one known simulated fault and one missing-data case,
plus raw CSV, replicate mean/std/SE, MAT proposals and a scientific summary plot.
This is software validation only; gamma bounds/numbers are illustrative and not
arm recommendations. Best individual noisy score is not confirmed performance.
Fixed-kernel latent uncertainty and weighted feasibility scores are not calibrated
physical guarantees. No real campaign configured, model/hardware executed,
controller/SLX/raw-data changes, commit or push in this development turn.

## 2026-10-27 - Manual tuning campaign workflow (milestone B)

User requested committing milestone A and executing the next development step.
Committed evaluator/reports/docs as e6fb02c on codex/pac-autotuning-framework.
Implemented milestone B under autotuning: manual campaign template, recorded
baseline extraction, bounded grouped-gamma proposal/review, operator manifest,
offline observation import, per-event append-only ledger, candidate/trial tables,
restart/status and screening/confirmation selection. Campaign MAT SHA-256 is
verified on each read; an exclusive writer lock and atomic event rename protect
normal writes. This is not a signed/tamper-proof database or automatic recovery.

validatePacTuneManualCampaign passes using explicitly synthetic MAT copies in a
temporary directory: candidate bounds/IDs, matching full tuning, two confirmation
repeats after screening, longer confirmation windows, metadata mismatch rejection,
withdrawal after failed confirmation, live-gain rejection, duplicate IDs/raw-file
hashes, exhausted budgets, restart, lock rejection, configuration-change detection
and campaign report export. A scalar ledger-table issue found on first ingestion
was corrected with explicit AsArray table conversion. validatePacTuneEvaluator
also passes unchanged v8/v9/v10/parser/coverage/fault regressions.

Saved reports/manual_demo_v1 exercises create/propose/import/status on historical
v9. Purpose is offline_demo, actual reset/payload details are unverified and
performance targets unset, so configurationReady=false and success=false. Its
illustrative bounds are software demonstration values, not gain recommendations.
No new experiment, model launch, controller/SLX edit, optimizer or hardware adapter.
Original MAT files were not modified. New milestone B source/docs/artifacts remain
uncommitted after the requested initial commit. Next milestone C is optimizer
integration after defining its API and testing against synthetic objectives.

## 2026-10-26 - PAC tuning framework milestone A (offline)

User authorized a new branch and framework development. Created
`codex/pac-autotuning-framework` from checkpoint 63f679a. Implemented MATLAB
intake/evaluation/report/replay functions under `passivity_adaptive_control/autotuning`.
No controller, tuning vector, SLX or raw recording was changed. No hardware model
was launched. No commit/push for this development step has been requested.

`validatePacTuneEvaluator` passes: duplicate-name selector enforcement, time-last
matrix decoding, singleton and gap alignment, constant Inf tuning comparison,
missing-data rejection, unset-target rejection and explicit roughness constraint
failure. Real-data checks match all v8 saved theta/q2 lap metrics within 1e-8;
v9 actual period 2 s, 65/1200 active and 1071/1200 limited, exact replay 1071 slew
target mismatches with zero allocation limits/flag mismatches and zero gate-active
disagreements. v10 period 4 s, fault 104.85 s, matched 36-44/52-60 s windows show
theta MAE .44945 -> .34942 while theta-rate SD .87627 -> 1.16627 rad/s.

Sparse constant switch logging initially caused false missing values; corrected
with an explicit constant-input lifetime policy for the three logged adaptation
switches, only when every native value is identical. Other inputs do not extrapolate.
Reports record which inputs were extended. Tuning changes are detected on native
timestamps. Bound/rate flags account for stored pre-update output timing.

Reviewed reports: `autotuning/reports/historical_reviewed_v1` (inventory, per-lap
metrics, tuning events, gate counts, diagnostics, evaluation MAT and summary PNG).
v8 is protocol_deviation (live gamma edits); v9 is invalid_data at the default
95% coverage rule (about 94.83% valid tracking), also with live edits; v10 is
aborted_or_faulted. These are evaluator classifications, not new conclusions
about experimental usefulness. All have success=false; performance targets
remain unset. Plots omit faulted/low-coverage tracking cycles and mark fault time.
No optimizer, hardware adapter, manifest/ledger, reset verification or repeated
confirmation is implemented yet. See the framework plan/README for next work.

## 2026-10-26 - v10 actual4s: cumulative compensation/roughness and separatefault
Offline1rad_10s_gainTesting_v10.mat:0-109.6s,4scirclereference(pi/2rad/s),
phi1,starts20s,recordingendsbefore30intendedcircles. Feedbacklambda18/16.5,
Ks15/10fixed. Gamma initiallyk.5,d.02,h5,bias1;d->.002at78.7s,k->.1at93.5s.
No gamma/feedback change atreported55sonset. Recordedcaps850/20/5/.2,
slew6bar/s,hupper20,lowerk2-500;otheruserprovidedoptionsconfirmed.
At55:validmodel/status3,noslewlimit/dropout. Estimatedd grows40/40at20s,
68.806/75.519at40,85.835/90.123at55,103.93/109.12at78.7,~105/110late.
Scale100makesgamma_d.02effective200,10x .002effective20. Loweringgammad
slowsupdatesbutdoesnotresetd. No damping/historyratecapcontacts;springratecaps
are rare overall. No configuredlimit-changeat55;h~3Nvs20bound.
k2firsthits-500at62s,afterreportedonset;thusboundcontactnotinitialtrigger.
K indefiniteintermittentlyfrom29.7s,andcontinuously55.2throughrecordingend,
minEigen-312.628N/m. Correlateswithchangingcompensation,butindefiniteKalone
doesNOTproveclosedloopinstability;noexactcausalidentificationfromoneclosedloop
log. Knownunmodeledactuator/derivativefilter dynamicsandchangedcommandslew
allowanceareplausiblecontributors. No actuatorpressuremeasurementsavailable.
Growingwithin-cycleoscillationcoexistswithdecreasingmeanerror:36-44thetaMAE
.44945,thetaRateSD.876rad/s,nearstop6.58%,phiMax.9998;52-60thetaMAE.34942,
SD1.166,nearstop16.77%,phiMax1.276;60-68thetaMAE.32410,SD1.413,
nearstop24.34%;84-92thetaMAE.27484,SD1.490,phiMax1.3787;94-102thetaMAE
.24428,phiMax1.4012. Ratesareapparentsampledangleincrements,notindependent
velocity. Do notusemeanerroraloneasstabilitycriterion. Forcecomponents
52-60:dampingRMS6.460/6.535N,history3.104/2.852,explicitfeedback.981/.646,
inertia.147/.153. Adaptive dampingforceisDh*vmodelcompensation,notKs damping.
Offlineallocatorreplay:20-44noallocationorslewlimiting;44-52oneslewlimit;
52-60sevenslewlimited;60-78.7allocation1,slew25;78.7-104.85allocation5,
slew59(withoverlap). Notv9-stylepersistentearlyblocking;limitinggrowslater.
Finalshutdownisaseparateidentifiedevent:tipinvalid104.75/.80/.85,
status4at104.75/.80,status5from104.85,norecoveryaftertipvalid104.90.
Enablepressureandkilllogged1,basevalid;age.05/.10/.15s;grace.10s exceeded.
95of96faultsamplesstillmodelvalid:forcefunctioncanbeevaluatedwhilepressure
supervisorfaulted;metricsmustrequirestatus3forcontrolledtracking. Mainaudit
lap22(104-108)includesfault/unloadingandisnotpuretrackingperformance;use
onset_audit.txtstatus3-onlywindows. Nohardwareorcontroller/modelchanges.
Suggesteddiagnosticfreshcomparisongamma_d.002fromstartup,other4ssettings
fixed,ifuserchooses;notvalidatedcure. CoupledSPDpolicyandactuatordynamics
needdesignreviewbeforeautomaticoptimization;donotrelaxboundsorlimitinggates.
Artifacts v10_analysis:audit,lapCSV,onset_audit,onset_timeline.png/fig,
parameter/gain/switchplots. Scripts tmp/auditPacGainTestingV10.m,
tmp/auditPacV10Onset.m. No commit/push.

## 2026-10-26 - User-requested repository checkpoint
User authorized staging/committing important project files,without push. PAC
SignedK2,HysteresisFreeze,Tunable,Adaptive and Pressure offline validations
passed. SignedK2 now also checks helper acceptance of user lower bound-500;
source block remains at-100. RL pressure/reference waveform validations also
passed offline. No model/hardware execution or fresh standalone codegen claimed.
Checkpoint includes saved user PAC model snapshot, raw circle recordings,
reproducible audit scripts/reports and presentation artifacts. Saved model is
not claimed compiled/tested by these source-only checks. Separate open-loop
collection tools/recordings are checkpointed independently. Unrelated paper
drafts,model .original backups,generated C build trees and duplicate presentation
ZIP remain local. Automatic optimization remains a proposal,not implemented.

## 2026-10-26 - v9 actual2s circles: slew limit blocks tracking and adaptation
Offline1rad_10s_gainTesting_v9.mat despitefilename:0-89.55s,30x2scircles
20-80s,phi1. Medianloggeddesiredangularratepi rad/s confirmsperiod2s.
Fixedlambda18/16.5,Ks15/10. Gamma initiallyk.5,d.02,h5.8,bias1;
k->.05at63.4s,h->.1at69.7s. Caps650/10/2/.2,lower1/-500/1,...,
hupper20,pressure2bar/sslewunchanged. No controllercodechangesperformed.
Exactgate reconstruction1200circlesamples:65active5.4167%;1071limited89.25%,
62invalid,2dt,nomanualdisable/deadband/independentresidualgate. Whole64invalid
54events/alltipinvalid,nofaults. Adaptactive falls12.75%first20stracking to
1.75%middle/last20s;someentirelapsnoupdates. Bias/history compensation remains
poorlylearned;reducinggammamidrundidnotremovelimit.
OfflinepacAllocateVariableCommonreplayon1138validstatus3circlesamples:
allocationinfeasible0,target-commandmismatch1071;ORmatchesloggedlimitedexactly.
Allocatorresidualpeak8.05e-15N,butpostslewcommandresidualpeak9.385N(nominal,
notphysicalmeasurement). Commandsmax[2.30469 2.55624 2.41175]bar,below3rails;
min[.988246 .8 .8]. Allcommandsreach2bar/s;fractionconsecutivecircleincrements
atabsoluteslewcap27.36/69.89/54.46%. Thus pressurecommandrate,notstaticforce
range,dominates thisrun. No inference chamberpressurefollowstheslewcommands.
Circlevaluebounds:k3=1on987samples82.25%,h1=0on1196(99.67%),h3=0on1178
(98.17%),h2=0on28;no uppercontacts/no k1,k2,d,biasboundcontacts. Veryfew
eligibleupdatesandphase-dependentprojectionlimitinterpretationoftheseweights.
Among65activecircleupdates,k2ratecap1,k3cap7,d2cap1,b2cap1;historycaps0.
These arenot thedominant causeofstalledlearning. Endtrackingparams
[1506.20 64.347 1 52.340 42.346 0 .91825 0 -.00143 .03171];do notsubstitute
afterreturnfinalvalues. EstimatedK remainsSPD:circleminEigen65.346N/m.
Final-lapthetaMAE1.067384rad,phi.275316,length8.069/9.760mm. Lastthreecircle
theta~1.08,phi~.28;no convincing sustained learning afterinitiallaps.
Recommendation:donotincreasegamma/capsordisableallocationgatesasfix. To
reachphi1/2s,characterizeactuatorbandwidthandcommandslewcap,thenreviseonly
withphysicalevidence;otherwise reduceamplitude/speed. Progressivefixedsetting
speedtestsorvalidatedactuatordynamicsneeded. Existinggatescorrectlyblock
updateswhenrequestedforceisnotappliednominally. No model/hardwarerun.
Artifacts v9_analysis:audit,lapCSV,allocation_and_bounds.txt,
two_second_slew_diagnostics.png/fig,gain/parameter/switchplots. Reproduction
tmp/auditPacGainTestingV9.m andtmp/auditPacV9Allocation.m. No commit/push.

## 2026-10-26 - v8 saturation audit and 2s-circle headroom
Offline1rad_10s_gainTesting_v8.mat:30x10scircles,phi1. Fixedlambda18/16.5,
Ks15/10. Gamma initiallyk.7,d.002,h5.8,bias1;k->.5at221.65s,d->.02at271.45.
Caps650x3/10x2/2x3/.2x2;lower1/-500/1/1/1/0x3/-8x2;hupper20. Allinputlogs
verified. Exactgate reconstructioncircle5769active96.15%;off62invalid+62dt+
107deadband,no pressure/slewlimitingorindependentresidualgate. Whole65invalid
singleevents/alltipinvalid,nofaults. Fastinitiallearningtheta<.2bylap8;
lap10theta.150591,phi.062523,length1.582/1.762mm;lap20.097513/.062452/
.943/1.340;lap30.099420/.061626/.911/1.359. Lateerrorlargelyplateaus.
Valuebounds onstoredpreupdatecirclesamples: k2=-500on24/6000(.4%);
k3=1on93/6000(1.55%);k1none. hzeroon16/5/17samplesearly;no huppercontact,
max12.178/11.440/11.652vs20. No damping/bias valueboundcontacts,no upper
boundcontacts foranyparameter. Final3cycles290-320:noneof10atbounds.
Returnholdisdistinct:k1lower91samples,k3lower135samples;do notmixwithcircles.
Ratecapactualcontacts:onlyk3,25activeupdatesbefore221.65s(~.42%wholeactive);
allotherparametercaps0;no ratecapcontactsafter221.65orinthefinal3cycles.
Thus remaininglateerrornotvalue/ratesaturation. Deadband107samples1.78%circles
notdominant;nominalresidualnotphysicalforcevalidation. Pressurecommandsmin.8,
max[2.67315 2.62012 2.48632]bar;floorcontactsnotnecessarilyinfeasibleforces.
Latevsv7:thetaRateSD.37460v.67930rad/s,nearstop3.952v27.428%,phiMax1.15769v
1.32153,limiting0v20.167%;positionaccuracyworsebutmotionmuchlessrough.
EstimatedKindefinite435circlesamples(431valid),7.25%,minEigen-347.062N/m;
notSPD/passivitycertifieddespitenofaults.
2sgoal:atfixedphi1,referenceangularrate5xandacceleration25xvs10s;nominal
13mmlengthamplitudeimpliespeak40.84mm/sand128.30mm/s^2. .05sPACgives40samples
percirclevs200. Simple5xcompressionofv8pressurecommands(notclosedloop2s
prediction)requirespeak~7.62/9.57/6.80bar/s. Inlatewindow~34.6/42.6/39.2%
commandincrementswouldexceedcurrent2bar/s. Actuator/commandbandwidthandlag
thereforemustbeaddressedalongsidelearning;raisinggamma/capsaloneno2sguarantee.
Recommendfixedv8candidateandprogressivespeedcomparisonwithsmoothreference
transitions/heldconditions,characterizingpressure/velocity/limitingandSPDfirst;
no autonomoushardwaretests,limitrelaxationorcontrollerchangeperformed.
Artifacts v8_analysis:audit,lapCSV,bound_and_rate_audit.txt,v7_v8comparison,
parameter/gain/switchplots;tmp/auditPacGainTestingV8.m,tmp/auditPacV8Limits.m,
tmp/comparePacV7V8.m.

## 2026-10-26 - v7 dynamic gamma sweep: fast errorlearning,large-d roughness
Offline1rad_10s_gainTesting_v7.mat,30x10scircles,phi1;lambda18/16.5,Ks15/10
fixed. Springcaps650(v6/v5:450),historyupper20,otherloggedbounds/historysame.
Gamma timeline k/d/h (biasgamma1always):initial.4/.0002/1.5;
51.35:.5/.0002/1.5;59.55:.5/.0002/1.8;81:.5/.0002/5.8;
114.6:.5/.002/5.8;126.1:.7/.002/5.8;170.15:.7/.02/5.8;
233.15:.7/.02/20.8;248.85:.7/.02/5.8;254.95:.7/2/5.8;
272.6:.7/.2/5.8;281.3:.7/.05/5.8;301.1:.7/.002/5.8.
Allchangesmidrunwithretainedstates,notindependentgaintrials. h5.8period
coincideswithfastererrorreduction:lap7theta.490726,lap10.303962,
lap14.154844,lap20.102254;gamma_kand_dalsochange,socausaloptimumunproven.
Final-laptheta.063644,phi.051364,length.787/1.013mm;lowerpositionerrorbut
latewithin-cycleroughnessworsens. Whole62invalidsingleevents/alltipinvalid,
nofaults. Exactgate reconstruction,circle93.5%active;189limited,61invalid,
53dt,87deadband. Deadbandnowblocks87samples,notallinactiveperiods.
At254.95d~49.952/52.396;gamma_d.02->2(100x) and d estimates reach
106.869/117.879by272.6,112.055/124.241by319.95. Latergammareductiondoesnot
reset learned d. Prior180-220thetaRateSD.33735,slow4.124%,phiMax1.1033,
limiting0;254.95-272.6SD.47959,slow13.256%,phiMax1.1861,limiting5.099%;
301.1-320SD.66565,slow25.474%,phiMax1.3215,limiting20.370%.
Final3290-320v7vsv6:apparentthetaRateSD.6793v.32025rad/s,slow27.428v4.103%,
phiMax1.3215v1.0855,l3velocityPeak18.102v12.270mm/s,limiting20.167v0%.
Thesearevalidsampleangleincrements,notindependentvelocitymeasurements.
AdaptiveDh*v ismodelcompensation,notidenticaltoexplicitKs dampingfeedback.
Gamma_dhugejumpandpersistenthighdstronglyassociatedwithroughness/limiting;
notjustgatingorimprovedconvergence. h20.8spans~1.6cycles,noisolatedbenefit
established;hmax12.559/10.995/12.086below20. k2neverhits-500;Kindefinite227
circlesamples(225valid),minEigen-172.699;notSPDcertified.
Suggestedfreshfixed-gaincomparison:lambda18/16.5,Ks15/10,gamma_k.5,
gamma_h5.8,gamma_d.002,bias1,retainingother v7settings. Candidateonly,not
validatedoptimum;nofurtherlivegammajumporcontinuingfromhigh-dstateforcomparison.
Artifacts v7_analysis:audit,lapCSV,gamma_stage_audit,gamma_and_damping_history,
v6_v7roughness/plots;tmp/auditPacGainTestingV7.m,tmp/auditPacV7Stages.m,
tmp/comparePacV6V7.m. No controller/model changes orhardwareexecution.

## 2026-10-26 - v6 springgamma.8: no clear benefit,more signedK excursions
Offline1rad_10s_gainTesting_v6.mat:30x10scircles,phi1. Springgamma.8verified
(v5.4);historygamma1.5unchanged. Fixedlambda18/16.5,Ks15/10fromstartup,
unlikev5earlygainchanges. Historyupperbounds20(vs10);nonbindingbecausemaxh
9.353/8.352/8.400. Otherloggedtuningmatchesv5. No sourcegainchangeperformed.
Gate reconstructionexact:circleactive98.1667%,55invalid+55dt,no limiting,
deadbandorindependentresidualgate. Whole57invalidsingleevents,alltipinvalid,
nofaults. Springratecaps450reached15.80/15.45/10.10%activeupdates(v5about
.21/.36/.53%). Thisisstrongerwithin-cycleupdateclipping,notproofbetterlearning.
v6 thetaMAE lap10.445369,20.302774,30.192586vs v5.443142/.286515/.188067;
phi.118820/.073161/.060439vs.131882/.072830/.064654. Finallength1.867/2.148
v1.884/2.201mm. Earlyrunsdifferfeedbackhistory/preconditioning;notclean
springgammaisolation. Late290-320s thetaRateSD.32025v.31694rad/s,nearstop
fraction4.103v3.242%,max2.232v1.544rad/s. l3velocityRMS2.302v2.368mm/s,
peak12.270v11.786mm/s;phimax1.0855v1.1069. Mixedchanges,notclearoverallgain.
EstimatedK indefinite774circlesamples(761modelvalid),12.9%vs v5 44samples
(.733%). CircleminEigen-288.615N/m;indefinitenessrecursuntil~239.65s,
notonlystartup. Wholemin-419.531;finalreturnk1=k3=1,k2=-420.528.
NoobservedfaultsdoesnotrestoreSPDguarantee. k2reaches-500early.
Recommendedcomparison:fixedfinalfeedbackfromstartupwithspringgamma.4,
historygamma1.5,andotherwiseexactv6settings;avoidraisingcapsorgammafurther.
Proposalonly. CoupledSPDconstraintremainsunimplemented andwouldbeaseparate
controllerdesignchange. Artifacts v6_analysis:audit,lapCSV,v5_v6roughness/plots;
tmp/auditPacGainTestingV6.m,tmp/comparePacV5V6.m. No model/hardwarerun,
controllercodechange,commit/push.

## 2026-10-26 - v5 smoother feedback settings,ongoing learning
Offline1rad_10s_gainTesting_v5.mat,30x10scircles,phi1. Actualgammaunchangedfrom
v4(history1.5),allothertuningunchangedexceptlivefeedback:lambda28/25.5initial,
28/16.5at40s,18/16.5at59.5s;Ks25/1initial,25/10at70.25s,15/10at88.2s.
Thus smootherlatebehaviorcannotbeassignedtogammachange;v2v3remainsgamma
comparison. After88.2feedbackfixedthroughend. Lap8thetaMAE.488633,phi.153605
->lap20.286515/.072830->lap30.188067/.064654;continuinglearningconsistent
withadaptation,notisolatedcausalproofwithoutdisabledrun. Final3meanlaptheta
v5.195919v v4.115522,phi.064290v.052724,length1.927/2.263v1.211/1.453mm.
Late290-320s smoothnessimproves: apparentthetaRateSD.316941v.504951rad/s,
nearstopfraction3.242v14.576%,max1.54431v3.28536rad/s;ratesaresampledvalid
angleincrements,notindependentvelocitymeasurements. l3velocitypeak11.786v
15.732mm/s,velocityerrorRMS2.368v2.564mm/s. Phi max1.10689v1.19559,
>1.05fraction7.589v11.765%. No circlelimiting(0v137samples),nominalresidual
gateindependent0;adaptactive97.933%(5876/6000),off62invalid+62dt,no deadband.
Whole64invalidsingleevents/alltipinvalid,nofaults. Gate reconstructionexact.
Historyweightsmax8.927/8.435/8.595,below10butstillgrowing;capsnotreached.
Promisingsmoothness/accuracytradeoff;livechanges/preconditioninglimitattribution.
Nextcomparisonshouldstartfinal lambda18/16.5,Ks15/10 fromfreshfixedbaseline
beforemoregammatuning,withsameinitialization/reference/payload. Proposalonly.
SignedKindefinite44validcirclesamples23.55-25.7s,minEigen-219.928N/m;
alsofinalreturn/hold. NoSPDguaranteerestored. No code/model/hardwareexecution.
Artifacts adaptCommonP/1rad_10s_gainTesting_v5_analysis:audit,lapmetrics,
v4_v5_roughness.txt,v4_v5_last3cycles.png/fig,parameter/gain/switchplots.
Reproduction tmp/auditPacGainTestingV5.m andtmp/comparePacV4V5.m.

## 2026-10-26 - v4 secondKs.1->1: small tracking gain,no clear smoothness cure
Offline1rad_10s_gainTesting_v4.mat:30x10scircles,phi1. Loggedlambda28/25.5,
Ks25/1 constant;allothertuningmatchesv3(historygamma1.5,caps2,lowerk2-500).
Gate reconstruction exact;circleadaptactive96%(5760/6000),off137limiting,
54invalid,42dt,7deadband. Whole58invalidsingleevents,alltipinvalid,nofaults.
Final-lapthetaMAE.112240,phi.058125,length1.228/1.481mm vs v3.118949/.052894/
1.200/1.543. Toavoidcherrypickinglastlap,last3mean-per-laperrorsv4theta.115522,
phi.052724,length1.211/1.453mm vs v3.121886/.052813/1.241/1.541:
theta~5.2%better,phiunchanged,lengthsmallimprovements,separaterunvariationunknown.
Last3velocityRMSEl3 2.564v2.699mm/s(~5%better);peak15.732v15.917mm/s.
ApparentthetaRateSD.50495v.51419rad/s(~1.8%better),slowfraction14.576v14.555%
(unchanged),maximum3.285v2.591rad/s(worseisolatedpeak). Phi max1.19559v1.22038
(lowerpeak),but>1.05fraction11.765v9.459%(longerovershoot);notuniformimprovement.
Last3pressurelimiting3.333v3%;wholecircle2.283v1.967%. Historyweightsmax
9.114/8.446/9.264below10;damping/historyratecapsnotreached. Gamma/capincrease
notindicated;secondKschangehasmodesteffect,notproofdampinghypothesisfalse.
EstimatedKindefinite71circlesamples(70valid),23.5-25.8and34-35.15s;
circleminEigen-243.848N/m;finalreturnmin-307.089. NoSPDguaranteerestored.
Suggestednextstep:retainfixedbaselineandinvestigatemeasured-velocityfiltering,
pressure-command/actuatorresponseandphase-lockedburstsbeforesweepingmoregains.
Artifacts adaptCommonP/1rad_10s_gainTesting_v4_analysis:audit,lapCSV,
v3_v4_roughness.txt,v3_v4_last3cycles.png/fig,gate/gain/parameterplots.
Reproduction tmp/auditPacGainTestingV4.m andtmp/comparePacV3V4.m.
No controller/model changes or hardware/model execution,commit/push.

## 2026-10-26 - v3 historygamma1.5 vs v2gamma1: faster learning,roughness persists
Offline1rad_10s_gainTesting_v3.mat,30x10scircles,phi1. Constantlambda28/25.5,
Ks25/.1. All logged tuning matchesv2 except historyweightgamma1->1.5; other
gammas/ratecaps/bounds/initials/historylaw unchanged. Payload/preconditioning
notestablished;separaterunsnotrepeatabilitystatistics. Gatedupdate reconstruction
exact. Circleupdates96.0167%(5761/6000);off118limiting+64invalid+54dt+3deadband.
Whole67invalidsingleevents,alltipinvalid,no faults. History/dampingratecapsnever
reached;peak historyrates.185/.198/.202N/s vs2caps. Springcaps .436/.520/1.292%
activeupdates. Historyweightsmax9.149/8.490/9.254,below10butapproachingbounds.
v3 thetaMAE lap10.283318,20.171381,30.118949 vs v2.325568/.229879/.168024;
phi.065628/.055170/.052894 vs.078437/.061870/.054390. Final lengthRMSEv3
1.200/1.543 vs v2 1.643/1.881mm. Finaltheta~29%better;phi~3%better. Firstlap
thetaMAE<.2 is17vslap24,~70searlier;single-runcomparison,notguaranteedspeedup.
Late290-320s roughness: l3filteredvelocityRMSE2.699v2.764mm/s;peak15.917v
16.067mm/s(almostunchanged). ApparentthetaRateSD.51419v.51998rad/s;
fractionapparentthetaRate<20%desired14.555v14.726%;maximum2.591v3.337rad/s.
Rate useswrappedangleincrements/.05s onconsecutivevalidsamples,notindependent
mechanicalangularvelocity. Phi max1.22038v1.16837;>1.05rad9.459%v6.25%valid
samples. Latepressurelimiting3%v2.333%;l3residualfrombestfit10sfundamental
1.003v.864mm (diagnosticfitonly,nofittedtracking). Thus averagepositionerror
improvesandlargestthetaburstreduces,butroughnessdoesnotclearlyresolveandphi
overshootincreases. Raisinggammafurthernotindicatedbythistradeoff;separate
secondKs/actuatordynamicsinvestigationwithfixedlambda/gammaisnextproposal.
SignedKindefinite69validcirclesamples23.65-25.75and33.95-35.2s,
circleminEigen-211.167N/m;returnholdwholemin-362.647. NoSPDguarantee restored.
Artifacts adaptCommonP/1rad_10s_gainTesting_v3_analysis:audit,lapCSV,
v2_v3_roughness.txt,v2_v3_last3cycles.png/fig,gain/parameter/switchfigures.
Reproduction tmp/auditPacGainTestingV3.m and tmp/comparePacV2V3.m.
No controller/model edits or model/hardwareexecution; no commit/push.

## 2026-10-26 - v2 catch-up burst near151s,not an error-zero stop rule
User observes phase-dependent roughlengthtracking/abrupttheta slowdown,mainlyl3.
Read-only inspection v2 at150.8-151.35s: modelValidtrue/status3 throughout;
one limitedsample150.95s withP2+0.1/P3-0.1bar per.05s and residual
[-.0965;.2769]N. No pressurezero/invalidholdatthisevent. Pressurecontinues
varying,allcommands~1.1-2.4bar,awayfromactivationthreshold. Filteredl3velocity
peaks.0168m/sat151 vsdesired.0039;then.0002by151.35while desired.0022.
Phi .9672at150.8 ->1.1342at151.05 ->1.0027at151.35(desired1).
Theta .2448at150.9 ->.4593at151 ->.4848at151.05,then~.486-.504through151.35:
burst then plateau,not evidencewholearmstops. Bothinferredlengthcoordinates
change;theyarenotindependentlengthsensors. Controller hasnoif(error==0)stop:
v=qd_dot-lambda*e becomesqd_dotwhen e=0 and force terms remain.
Catch-upovershoot/braking withmodel/actuatorlag is plausible,notprovenmechanism.
Known derivativefilter lag can delay velocityfeedback; rawposeartifactversus
physicalstick-slip remainsunresolved. CurrentexplicitKs=[25;.1] makes -Ks*s
secondcomponentweak:at151sfeedback~[-1.438;-.00168]N. This isnotfullclosedloop
damping;Mlambda andothermodeltermsalsoaffectresponse. Highlambda can increase
transientdemand;do not assumeincreasinggamma removeswithin-cycleroughness.
Nextgamma comparison shouldretainfeedback/reference/filtersettings andassess
velocityspikes,phiovershoot,thetavelocity,andlimitingaswellasMAE. A separate
modestsecondKscomparisonmaytestdampinghypothesis;notimplementedorvalidated.
No code/model changes or hardwareexecution.

## 2026-10-26 - Fixed-gain gainTesting_v2: gradual error reduction
Offline1rad_10s_gainTesting_v2.mat:0-330s,30x10scircles,phi1,payloadunverified.
Allsettingsconstant:lambda28/25.5,Ks25/.1;gamma.4x3/.0002x2/1x5;
caps450x3/10x2/2x3/.2x2;scales1000x3/100x2/1x5;lower1/-500/1/1/1/0x3/-8x2.
Adapt/common true,hfreeze false;optionsdeadband.005,tol.02,maxdt.075.
Gate reconstruction exact.6000circlesamples:5771active96.1833%;105limited,
69invalid,55dt;no independentdeadband/residual/disable/unexplained. Whole
70invalid69events,69tipinvalid,no faults;one modelinvalidsample not accounted
for by tipinvalid alone. Do not attribute everyinvalidsample totipsensor.
ThetaMAE lap2.435323,lap10.325568,lap20.229879,lap30.168024rad;
phi.138773/.078437/.061870/.054390;lengthRMSElap2:3.989/4.191mm,
lap30:1.643/1.881mm. Theta/phi each~61%lower lap30vslap2. Five-lapmean
theta1-5:.413708,6-10:.351219,11-15:.293805,16-20:.247350,
21-25:.208504,26-30:.177861. phi:.132535/.086546/.070767/.062998/.059407/
.055848. Continuedthetaimprovement;phidecay slows;no settledfinalfloorproved.
Springratecap fractions activeupdates.687/.285/2.044%;damping/history0%;
biascaps.251/.184%. Peak historyrates.122/.132/.141N/s vs2N/scaps;peak
dampingrates.095/.098 vs10. Raising caps cannot accelerate these uncapped
updates;deadbandreductionnotindicated. Gradient gains are potential controlled
tuning variables,not evidence of guaranteed improvement from larger gamma.
Historyweightmeanslap10[3.216 3.004 3.167],lap20[5.628 5.240 5.572],
lap30[7.275 6.778 7.307];stilladapting,notat10upperbound. Slowingthetaerror
coexistswithongoingcompensation;noidentifiedphysicalparameterconvergenceclaim.
k2cycles;min-500duringlap1,final-afterreturn-217.498;k3hits1. EstimatedK
indefinite81circlesamples,79modelvalid,intervals23.4-25.9and33.9-35.35s,
minEigen-257.412N/m. Observedboundedtracking/nofaultsdoesnotcertifySPD.
No controllerchanges. Artifacts adaptCommonP/1rad_10s_gainTesting_v2_analysis:
audit,lapCSV,gatefigure,parameter/gain/switchplots;tmp/auditPacGainTestingV2.m.
Possible next experiment: modest historygamma1->1.5 withallotherinputsfixed,
compareearlylapsandlimiting/bounds;proposalonly,notimplementedorhardwaretested.

## 2026-10-26 - rad_10s_gainTesting gate and gain audit
Offline0-330s,30x10scircles,phi1;payloadnotestablished. All tuning logged,
history options namedhystoryOptions. Gamma[.4 .4 .4 .0002 .0002 1 1 1 1 1],
scales[1000x3 100x2 1x5],caps[450x3 10x2 2x3 .2x2],lower
[1 -500 1 1 1 0 0 0 -8 -8];these differ from earlier supplied history settings.
Mechanical/commonenabled,hfreeze false,deadband.005,residualtol.02,maxdt.075.
Actual lambda starts21/20.5;changes110.35:21/18.5,122.20:28/18.5,
166.10:28/12.5,201.65:28/8.5,219.30:28/38.5,280.55:28/28.5,
321.90:28/25.5. Ks starts25/17;changes101.85:25/15,144.65:25/1,
176.60:25/10,243.75:25/100,250:25/.1,272.20:25/.01,298:25/.0001,
308.85:25/.1. Gains change midcycle; adaptive states evolve throughout.
Gate reconstruction matches all6601 adaptationActive values. Disjoint first
failure counts for6000circle samples:5534active(92.233%),358limiting,
67invalid/model,41updateinterval;zero independentdeadband,residual-only,
disable or unexplained. Whole run866inactive includes399supervisornotcontrol,
68invalid,358limited,41dt;no faults. 68invalidsamples59events,alltipinvalid.
Kssecond100 interval243.75-250:123/125sampleslimiting,active.8%;thetaMAE.314676,
phi.152063 vs prior219.3-243.75:.209407/.079024 and next250-272.2:
.213122/.059409. Higherfeedback can block learning via limiting;not too-small
gamma. Springratecap fractions amongactiveupdates.401/.140/1.622%;othercaps0.
Final-laptheta.187799,phi.060273,lengthRMSE1.967/1.904mm. Late improved
tracking does not establish optimal gains because adaptation/order/speedphase
confound partial intervals. Need fixed-gain repeated comparisons.
Signed K indefinite on29validcircle samples23.95-25.35s,minEig-87.667N/m;
whole-runminEig-170.793 duringreturn/hold. modelValid doesnotcheckKSPD.
No observedfaults isnotSPD/passivityproof. k2min-318.435 neverhits-500;k3hits1.
k2 cycles,final-afterreturn-171.791;do nottreatthissinglevalueasconstantcirclefit.
Artifacts adaptCommonP/rad_10s_gainTesting_analysis: audit,lapCSV,gatefigure,
parameter/gain/switchfigures;reproduction tmp/auditPacGainTesting.m.
No controller/model edits,hardware/model execution,commit/push.

## 2026-10-26 - User-authorized signed k2 experiment
User confirmed k2 reaches its bound even after manually lowering it to1 and
requested -100 N/m. Updated pacTuningBlock lower entry2 to-100 and relaxed
pacTuningConfig sign assert only for k2,to lower(2)>=-100. Other constraints,
initial estimates,adaptation law/rates,wrapper ports unchanged. Legacy positive
defaults retained. A bound hit demonstrates gradient direction,not that physical
k2 requires negative stiffness. Combined K SPD requires k1+k3>0 and
k1*k2+k1*k3+k2*k3>0;relaxed box can violate this. No coupledSPD projection or
invalid-K fault added. No new passivity/stability claim for signed configuration.
validatePacSignedK2 passed negative update/clamp,outwardhold,inwardrecovery and
examples showing signed k2 can yield either SPD or indefinite K. Existing
validatePacTunable/Adaptive/HysteresisFreeze passed; their positive-default SPD
check does not certify the signed tuning box. No model/hardware execution,
SLX edits,commit or push. TUNABLE_BLOCK documents integration and limitations.

## 2026-10-23 - 100 g / 10 s circle presentation figures
Follow-up: plotPacAdaptationPresentation.m exports allten parameterEstimate
histories,lambda/Ks,and adaptation/common/update Boolean histories in sameformats.
Dottedgainchange markers94.00/108.60s. Usesretainedparameters,notmatrixoutput
dropouts. Tuninggamma/scales/ratecapsabsentfromsavedlogs,notreconstructed.
Read-only 1rad_10s_tr2_100g.mat,0-330s,thirty10s circles,phi1rad,userreported100g.
Added plotPacPresentation.m for PNG300dpi,vectorPDF/SVG,editableFIG exports;
five figures per window: XYZtime,pressure,length,theta/phi,3DXYZ. Full recording
and290-320s finalthreecircles saved in adaptCommonP/presentation_1rad_10s_tr2_100g.
XYZ follows existing sensor coordinate convention: measured loggedposition,
desiredCCgeometry[.055 .17411 .05]m includes base/tipsensoroffsets; no fitted
alignment,phase shift or smoothing. Physical module endpoint is not substituted.
Invalidsamplesmasked;wrappedtheta brokenatwraps;lengths labeledchanges,pressure
labeledcommands. Actual Ks23/10->23/30at94s;lambda18/16.5->21/18.5at108.60s.
Final3circleMAEtheta.373505,phi.100557rad;XYZRMSE21.925/21.729/18.501mm;
lengthchangeRMSE3.378/3.295/3.698mm. These metrics are descriptive,not controlled
payload/adaptation comparisons. No controller/model changes or hardware execution.

## 2026-10-16 - Trial3 hysteresis-only freeze physical recording reviewed
Offline adaptCommonP/1rad_10s_tr3.mat:0-330s,6601 samples,30x10s circles.
lambda=[18;16.5],Ks=[23;10] constant,matching trial2; adaptation/common true.
freezeHysteresisWeights loggedfalse->true at121.11s;first PAC sampletrue121.15s.
Weights held exactly[1.5302619 1.3560675 1.3215875] afterward,zero range.
History-state postfreeze ranges[1.32867 1.32974 1.33001];other parameter
updates continue,adaptationActive98.04% overpostfreeze samples. Confirms intended
hardware-recorded switch behavior; does not resolve prior standalone Coder error.
Valid-sample history forceRMS lap12[1.3779 1.3461],lap30[1.3867 1.3436]N:
history force amplitude stabilizes with weights; history dynamics remain active.
Mean springs lap10[1042.01 126.511 660.584],lap20[1019.73 129.45 661.5],
lap30[1030.16 129.289 661.794]N/m: earlier downward means largely plateau,
although cycle-periodic changes persist. Damping/bias estimates still evolve.
Trial3 lap10thetaMAE.54238,phi.15918,length4.425/5.158mm;lap30theta.51956,
phi.14982,length4.281/4.900mm. Trial2 continuouslyadapting lap10.52496/.15411/
4.391/4.967;lap30.41577/.10194/3.520/4.026. Improvement lap10->30 in theta
is~.023rad frozen vs~.109rad adapting; phi~.009vs~.052rad. Frozen trial already
has somewhat worse prefreeze tracking,so separate-run differences are not a
perfect causal/repeatability-controlled estimate. Tuning vectors still absent
from saved names; payload/preconditioning/trim equivalence not established.
Evidence supports a useful contribution of continuing history adaptation under
these conditions,not that increasing h is physical hysteresis identification.
Joint history growth/stiffness decline is consistent with compensation tradeoffs;
freezing limits redistribution and does not cure the remaining directional lag.
No basis yet for permanently freezing weights or aggressive leakage/increasedgamma.
58 separate one-sample dropouts(57duringcircles),allmodelinvalid/status4/tipinvalid;
no faults or circle limiting. Circlecommon1.6-1.63499bar;commandpressuremin
[.914123 .8 .8],max[2.51173 2.45825 2.39412]bar;nominalresidualpeak9.41e-15N.
Freeze does not remove invalid-tip events; retain validity masking for forces.
Reproduction tmp/auditPacFreezeTrial3.m,report tmp/trial3_freeze_audit.txt,
figure tmp/trial3_freeze_comparison.png. No controller/model edits or model runs.

## 2026-10-16 - Hysteresis-only freeze implemented (offline)
User requested paste-ready code to freeze only hysteresis-weight adaptation.
pacAdaptiveTunableBlock appends Boolean input23 freezeHysteresisWeights;
true masks effective cfg.gamma(6:8) tozero before update,without changing the
configured input gamma. False/default for older MATLAB calls retains behavior.
No state resets; z history still evolves; h*z force still evaluated; remaining
parameter adaptation retains original validity/allocation/supervisor gates.
enableAdaptation=false still disables all estimate updates. No output changes.
validatePacHysteresisFreeze passed synthetic wrapper checks: learned nonzero
weights held exactly,history/force evolving,other adaptation active,resumption,
global disable. validatePacTunable/Adaptive/Pressure also passed. These checks
do not establish physical sampled-loop stability/passivity. No SLX edits,
hardware execution,commit or push. TUNABLE_BLOCK.md documents Boolean Step120s
for freeze afterten circles starting20s; log freeze input and tuning values.
Existing invalid-sample adaptiveForce diagnostic issue intentionally unchanged.
Additional standalone MATLAB Coder checks (GenCodeOnly and explicit -c) emitted
C sources but failed at reported build stage with 'Index exceeds array bounds'
and 'C compiler produced errors'. Cause not isolated; code-generation/build
success is NOT established by the passing MATLAB behavioral checks. No Simulink
model was compiled or executed. Retained diagnostics under tmp/pacFreeze_codegen
and tmp/pacFreeze_source for follow-up if integration encounters the same issue.

## 2026-10-16 - Trial2 constant-gain 30-cycle audit
Read-only adaptCommonP/1rad_10s_tr2.mat: 6601 samples,0-330 s,30 circles
20-320 s at10 s/circle. Actual constant lambda=[18;16.5],Ks=[23;10];
adaptation/common switches true. User reports logging tuning block outputs,
but adaptationGamma/parameterScale/parameterMaxRate and other tuning vectors
are absent from logsout; out contains NDI,dynamic_NDI_meas_len,
kin_NDI_meas_angles,kin_NDI_meas_len,logsout,simout,tout. simout is an unnamed
27-channel force/diagnostic timeseries, not a named ten-vector tuning dataset.
Do not substitute current source defaults: observed k2/k3 rates reach250,
whereas current pacTuningBlock.m source caps are200.

Exactly48 separate one-sample (.05 s grid) all-zero damping/force-evaluation
events;47 during circles,one at11.05 s. All48 match modelValid=false,status4,
configurationValid=false,orientationValid=false and sensorValid=[true false].
Base orientation valid throughout. Both logged frame pairs remain finite,
advancing and matched at these events: not repeated/mismatched-frame rejection.
First six events raw1-4 finite/raw5-8 NaN; invalid tip channel established by
sensorValid at all48. Measurement age=.05 s at drops vs~0 at good samples;
this is below current source .075 s threshold,not evidence of stale-age failure.
Inertia/Coriolis/damping/static/feedback/history/bias/request outputs all zero
at48/48 events because pacAdaptiveForce returns initialized outputs on invalid
input. Pressure commands exactly unchanged from preceding sample at all48.
No faults; no circle limiting; valid-circle nominal residual peak9.67e-15 N.
Stored parameter/history changes from each invalid sample to next sample are
exactly zero: no reset/update on those evaluations. Physical origin of missing
tip data (visibility/device/driver etc.) not determined by these logs.

IMPORTANT diagnostic defect: adaptiveForce_N is computed unconditionally in
the tunable wrapper after force evaluation. At invalid samples Kh=Dh=0 and
fh=fb=0, so it evaluates -Knom*q-Dnom*v,not valid adaptive compensation.
It is nonzero at48/48 drops (peak abs[11.2 15.55] N). These spikes must be
masked by modelValid; equivalent force/pressure is held,not zero. No fix applied.

Lap2/10/20/30 mean k coefficients:
[1123.1 144.9 700.18]/[1011.9 138.16 697.23]/[923.65 125.26 630.63]/
[846.78 122.15 567.26]. k2 touches50; final after return589.122/84.815/480.498.
Same laps mean h=[.23555 .17851 .19261]/[1.2415 1.1499 1.1183]/
[2.3935 2.2634 2.1867]/[3.4565 3.2953 3.1896] N; zRMS stays~.51-.55.
History forceRMS[.20277 .18676] ->[3.2046 3.0601] N. Thus increasing history
force comes from increasing weights,not unbounded history states.
Mean negative history-feature/error products remain positive in20-120,
120-220,220-320 s: [.02635 .02440 .02321], [.02331 .02251 .02179],
[.02149 .02085 .02033]. Current gradient law therefore keeps raising weights;
no leakage or convergence demonstrated. Late spring-feature mean gradients
are negative,consistent with declining mean k. This shows the update mechanism,
not physical parameter identification. History features correlate with measured
muscle velocities(~.64-.66); missing dynamics and compensation tradeoffs remain
plausible explanations,not established sources. Raw state feature scaling and
different units preclude directly comparing gradient magnitudes across groups.
No inference that all stiffness FORCE components decrease: adaptive stiffness
correction RMS lap2[4.458 3.308] vs lap30[1.877 4.257] N because configuration
and nominal subtraction matter. Hysteresis and stiffness channels jointly change.
Lap2thetaMAE .55129,phi .20651,lengthRMSE4.283/5.626mm;lap30 .41577/.10194/
3.520/4.026mm. Gradual improvement coexists with continued parameter drift.
Recommend valid-masked diagnostics, verified tuning logging, then a matched
history-weight freeze comparison; leakage/residual dynamics need separate design.
Offline reproduction tmp/auditPacTrial2.m; report tmp/trial2_audit.txt;
figure tmp/trial2_diagnostics.png. No controller/model changes or hardware runs.

## 2026-10-16 - Hysteresis-force growth in 10 s circles
Read-only 1rad_10s_tr1.mat separates history-state amplitude from weight drift.
History z is already approximately periodic by lap2 and remains within about
+/- .666; lap2 RMS [.5428 .5642 .5602], lap30 [.5754 .5677 .5650].
Mean adaptive history weights h1/h2/h3: lap2 [.2181 .1945 .1670], lap10
[1.2098 1.2487 1.1281], lap20 [2.2759 2.3718 2.1736], lap30
[3.2669 3.3975 3.1322] N. Hysteresis-force RMS grows from [.2096 .1676] N
on lap2 to [3.2423 2.9769] N on lap30; final-lap extrema approximately
[-4.4165 -4.2385] to [4.4578 4.1915] N. Growth comes predominantly from
adaptive weights, not hidden history-state accumulation. Current update uses
projected/rate-capped -scale.^2.*gamma.*Y'*s with no leakage; persistent
cycle-averaged gradient drives weights upward. Missing lag/force-map dynamics
or parameter tradeoffs are possible explanations, not identified causes.
No evidence of weight convergence within this recording; not proof of physical
hysteresis growth. Live gain changes confound interpretation. Increasing history
gamma may accelerate this drift; leakage/regularization/freeze comparisons are
future design options requiring controller analysis, not implemented fixes.

## 2026-10-16 - 10 s circle run with live feedback tuning
Offline 1rad_10s_tr1.mat: 0-330 s, 6601 PAC samples at .05 s; reference
phi=1, thirty 10 s circles over 20-320 s. Both adaptation/common switches
record true. Payload condition not established by filename/logs.
Actual lambda: [21;15] initially; second component changes to 21 at81.60,
18 at97.45,10 at172.10,17 at188.50,16 at211.75 s. Actual Ks: [25;20]
initially; second component 30 at105.75,28 at198.30,25 at238.20,20 at249.05,
18 at268.05,30 at298.75 s. Thus lap improvement is confounded by live tuning.
Lap2 thetaMAE .47249, phiMAE .18991, reduced lengthRMSE3.487/5.408 mm;
lap10 .40913/.11003/3.491/3.971; lap15 .37870/.09886/3.238/3.756;
lap28 .35557/.09903/2.903/3.750; lap30 .32598/.08923/2.693/3.474.
Metrics use valid model/orientation and observable theta; theta wrapped.
Theta signed mean equals MAE in every lap: persistent directional lag.
Lap30 equivalent lag .519 s at .62832 rad/s is not identified transport delay.
No latched faults; 55 zero damping-output samples all match modelValid=false
and status4; retained damping estimates stay >=40. Consistent with screenshot
duration/dropouts, not parameter resets. Circle valid98-100%, adaptation96-100%.
Only lap17 has limiting (.5%, residual peak .29728 N); other circle residuals
roundoff. Common1.5925-1.6903 bar; pressure min[.86471 .8 .8],
max[2.5416 2.4994 2.3901] bar.
Peak parameter rates [194.002 200 200 .037625 .040522 .050174 .040396
.042859 .13539 .139149]; k2/k3 reach current source caps, others below.
Gamma/rate-cap inputs absent, so exact experimental tuning not established.
Final parameters [690.71 97.702 460.83 42.999 43.015 3.3349 3.4664 3.1965
-.48048 .034948]; k2 reaches50 during run. Estimates are compensation values.
Compared with 20 s run final lap: theta .326 vs .288, phi .089 vs .064,
length2.693/3.474 vs2.536/2.717 mm. Different gains, speed and number of laps
preclude attribution to speed or adaptation alone. No hardware/model executed.
Reproduction: tmp/auditPac10s.m reads MAT only and prints metrics/gain changes.

Results below summarize existing recorded evidence; migration did not rerun
MATLAB tests or operate hardware. Distinguish calculation checks, model smoke
tests and measured behavior. Full settings/provenance remain in linked local
READMEs, original reports and dated log entries.

## NDI / estimator evidence — 2026-10-01

Straight calibration: 100/100 accepted pairs, 2.30 mm XYZ reference error;
separate preflight 98/100 pairs, median interval 0.0502 s. These checks do not
establish global orientation/strain accuracy after mounting changes.

Six 30 s sine trials, 601 samples each, ten invalid orientation samples among
3606. Filtered difference versus four-state CV Kalman used 0.05 mm measurement
SD and 0.003 m/s² per-step acceleration SD. At tau=0.12 s the mean prediction
proxy was about 0.106 versus 0.110 mm; no accuracy winner established.

| Filter tau (s) | Low-command velocity RMS (mm/s) | Prediction RMS (mm) | Velocity-change RMS (mm/s) |
| ---: | ---: | ---: | ---: |
| 0.12 | 0.0633 | 0.1061 | 0.6135 |
| 0.20 | 0.0534 | 0.1027 | 0.4064 |
| 0.30 | 0.0464 | 0.1011 | 0.2868 |
| 0.50 | 0.0381 | 0.1009 | 0.1827 |

Status: offline comparison completed; true velocity and controller phase lag
unvalidated. [Method and saved outputs](../estimator_comparison/README.md).
Source recording hashes and per-trial baseline metrics remain in the October 1 log.

## Early feedback implementation — 2026-10-01

Synthetic checks passed missing-frame hold, sine derivatives, error sign,
nonnegative/bounded pressures, force reconstruction, saturation and invalid-frame
zero command. Hardware-free copy with NI/NDI removed compiled and ran 0.5 s,
pressure 0–0.8002 bar. Status: narrow offline smoke pass; precedes later routing
and dashboard revisions and does not certify the current model or plant.

## Inverse and stiffness evidence — backfilled 2026-10-05

- Inverse recording `simulationresulta_20261001_1936.mat`: 2659 samples,
  167.836 s physical time, 99.85% valid, 7.48% limited, three qualifying holds,
  two stationary tails, one target maintained under 1 mm / two-second rule.
  Status: analysis completed, partial tracking success. Zero-rate holds cannot
  validate dynamic terms. [Report](../inverDynamic_test/results/report.txt),
  [method](../inverDynamic_test/README.md).
- Stiffness sweep `Comparision_results.mat`: grid 1350 down to 200 N/m;
  best tested 650, reduced-coordinate RMS 3.55231 mm over 600 valid samples;
  nominal baseline replay discrepancy <=4.61e-7 m. Status: replay consistency
  pass and conditional same-recording fit, no independent physical stiffness
  identification. [Report](../FWD_INV_COmparision/stiffness_sweep_results/report.txt).
- P1 static fit: means about 4132 N/m at phi=0.09425 to 953.333 at 2.05594 rad;
  low-angle SDs 1843 and 2597 N/m. Zero deadzone, pooled loading/unloading,
  command pressure assumed. User reported poor unloading agreement using phi
  alone. Status: empirical interpolation available, general law unvalidated.
  [Method and outputs](../planar_stiffness_test/README.md).

## Existing tracking / plotting — 2026-10-04–05

`testing_20261004_1311.mat`, 95–125 s: desired/measured phi=1.768/1.049 rad,
mean reduced errors [+9.37,-4.01] mm, requested force [16.01,0.39] N,
command [0,2.006,0.829] bar, limiting false. The inspected inverse used K=250,
not later K=650. Status: substantial steady error; modeled force consistency
does not demonstrate delivered force accuracy.

`plotCircleTracking.m` produced seven figure groups offline for three saved circle
recordings and an out-object input. Status: plotting execution pass, no tracking
performance claim; orientation-derived versus measured XYZ may differ.

## Common-pressure 300 s — 2026-10-05

6001/6001 valid; physical interval 0.0508 s, span 304.58 s. Final-five-second
tails at common means 1.2/1.6/2.0 bar gave neutral U=0.1065/0.1702/0.1958 and
V=0.0717/0.1444/0.1523 rad. Measured differential U gains
0.3096/0.4649/0.3897 rad/bar versus modeled 0.9641. Largest tail U range
0.00477 rad. Status: ideal equal-pressure neutrality disagrees with measured
pose; software commands matched regulator-path logs, chamber sensing absent.
Apparent K=2199/1464/1747 N/m is a conditional diagnostic.
[Report](../common_pressure_testing/results/test_300s/report.txt).

## Neutral-trim 387 s / validation 199 s — 2026-10-05

387 s: 7741/7741 valid; reference U,V=(0.03051,0.01501) rad.
Local Jacobian condition numbers 1.57/1.39/1.47 support distinct command
directions at the three means. Proposed triplets [1.328,1.292,0.981],
[1.766,1.766,1.268], [2.236,2.172,1.592] bar extrapolate beyond tested probes.
Status: candidate estimation only. [Report](../common_pressure_testing/results/trim_test_387s/report.txt).

199 s: 3980/3981 valid. The mean-1.6 candidate overshot reference with U/V
residuals (-0.10136,-0.15204) and (-0.09697,-0.14795) rad on two approaches.
Intermediate [1.7,1.7,1.4] residual (+0.01795,-0.01143), XYZ norm 6.57 mm
versus candidate ~26 mm. Status: extrapolated candidate failed neutral target;
repeat results do not prove equilibrium. [Report](../common_pressure_testing/results/trim_test_199s/report.txt).

## Manual trial-4 cycling — 2026-10-05

635.25 s, 12706/12706 valid, 17 on/17 zero episodes, four triplets.
Best of four was [1.6656,1.6492,1.4852] bar. Five later repeats:
U=0.02293 (SD 0.00103), V=0.00690 (SD 0.00247), magnitude=0.02404
(SD 0.00107) rad. Relative bend change 0.0047–0.0080 rad and XYZ shift
3.8–4.3 mm. Status: within-run local candidate, not exact straightness;
first approaches differ and short holds may still drift.
Source: `common_pressure_testing/manual_test_results.mat`; dated log entry.

## Local calibration 395 s — 2026-10-06

7900/7901 valid. +/-0.04 bar s/a probes around trial 4 gave approximate
`J=[[0.397,0.069],[-0.079,0.490]]` rad/bar after neighboring-trim subtraction.
Five identical-command trim holds spanned U/V=0.0096/0.0138 rad and differed
from manual-run means; tail slopes up to ~0.00035 rad/s. a+ achieved
(U,V)=(0.0067,-0.0007), but at 237.30 s an unexplained ~0.0063 rad
orientation-vector step occurred without a corresponding large XYZ/command
change. Status: useful local directions, equilibrium/repeatability incomplete.
Source: `common_pressure_testing/local_calibration_395s.mat`; dated log entry.

## Bounded feedback numerical checks / first recording — 2026-10-06

Pure-source checks at 0.05 s passed nominal bounds, zero windows, normal
0.8 bar/s slew and reference phases; extreme-acceleration case reached
[1.6656,1.6892,1.4452] bar, mean 1.6, limited=true. Status: nominal numerical
pass only; does not cover the subsequently observed invalid-sample jump.

Physical recording: 240 s, 4800/4801 valid; Kp=[1100;1100], Kd=[100;100].

| Window (s) | Desired U,V (rad) | Measured U,V (rad) | Mean length-error norm (mm) |
| --- | --- | --- | ---: |
| 115–125 | 0.0100,0 | 0.0054,-0.0065 | 0.112 |
| 185–195 | 0,0.0100 | 0.0016,-0.0032 | 0.211 |

Relative U response +0.0016, V response +0.0036 rad; small final slopes
~<=0.00006 rad/s support persistent residuals. Mean pressure stayed 1.6;
ds range -0.00749..+0.00288, da +0.00390..+0.00973 bar; no limiting in
3300 feedback samples. Initial trim drift confounds claimed zero-target improvement.
At 14.40 s invalid NDI produced an immediate zero command, ~33.3 bar/s;
recovery 0.8 bar/s. Status: intended response direction, incomplete tracking;
fault-path slew requirement violated. No post-Kill logs in this MAT.
Source: `common_pressure_testing/feedback_testing_results.mat`; latest log entry.

V2 source-only fault injection on October 6 passed: a single 0.05 s invalid
frame held the prior pressure, valid input resumed without a pressure step,
a sustained invalid interval slewed toward zero at no more than 0.8 bar/s,
and a latched fault did not automatically re-pressurize. Executed with MATLAB
R2025a using `testCommonPressureFeedbackInverseBlockV2.m`. Status: offline
source check only; no Simulink compile or NI/NDI run. The saved model still
contained V1 at the time of that source check; later saved-model inspection
finds V2 in the main inverse chart.

## V2 physical recording — 2026-10-06

`common_pressure_testing/feedback_testing_results_V2.mat` spans 240 s with
4800/4801 valid NDI samples; the one invalid sample was at 12.65 s. Both gain
logs stayed at Kp=[1100;1100], Kd=[100;100], so the proposed gain increase
was not part of this test. The desired U/V steps remained 0.01 rad.

| Hold (s) | V1 measured U,V (rad) | V2 measured U,V (rad) | V1/V2 length-error norm (mm) |
| --- | --- | --- | ---: |
| U, 115–125 | 0.00543,-0.00649 | 0.00415,-0.00760 | 0.112/0.133 |
| V, 185–195 | 0.00156,-0.00324 | -0.00040,-0.00555 | 0.211/0.248 |

The V2 zero-target 80–90 s pose was U,V=(0.00168,-0.00555), versus
(0.00379,-0.00530) in V1. Relative to that local baseline, the U step moved
U by +0.00247 rad in V2 versus +0.00164 in V1. Relative to the 150–160 s
zero-target hold, the V step moved V by +0.00221 in V2 versus +0.00360 in V1.
Both runs retained offset; different trim poses/history prevent attributing
the larger V2 absolute errors to V2 fault handling.

V2 `logsout` has 48 elements versus 56 in V1. It omits
`pressureCommand_bar`, `pressureLimited`, `inverseModelValid`,
`tauAvailable_N`, `tauDesired_N`, and the three inverse force outputs.
No post-Kill pressure command is logged. Thus the MAT
cannot establish whether pressure held through the invalid 12.65 s frame or
whether the local correction bounds were active. Read-only inspection of
`FBL_commonP.slx`, saved before the MAT file, found V2 in the main
`MATLAB Function1` chart and V1 in `MATLAB Function8`. This supports the
intended source version but is not a measured pressure-path check.

## Higher-gain V2 recording near straight — 2026-10-07

`common_pressure_testing/feedback_testing_results_V2_newGains.mat` spans
240 s at 0.05 s logging intervals. It records Kp=[1400;1400] and
Kd=[120;120], constant, and the same 0.01 rad U/V references. NDI and inverse
model were valid for 4800/4801 samples. At the single invalid sample, 22.00 s,
controller and logged post-Kill NI commands held exactly
[1.6656,1.64918,1.48522] bar from the prior sample and resumed at the same
value on the next valid sample. All three logged NI-command channels matched
controller pressure throughout. The Kill signal was logged enabled at its two
event samples. Actual chamber pressure was not measured.

| Hold (s) | Desired U,V (rad) | Measured U,V (rad) | Mean reduced length-error norm (mm) | Mean controller P1,P2,P3 (bar) |
| --- | --- | --- | ---: | --- |
| Zero, 80–90 | 0,0 | 0.00102,0.00108 | 0.020 | 1.6758,1.6459,1.4783 |
| U, 115–125 | 0.0100,0 | 0.00422,0.00449 | 0.090 | 1.6680,1.6502,1.4818 |
| Zero, 150–160 | 0,0 | 0.00143,0.00152 | 0.028 | 1.6759,1.6460,1.4782 |
| V, 185–195 | 0,0.0100 | 0.00157,0.00167 | 0.134 | 1.6759,1.6389,1.4853 |
| Zero, 210–220 | 0,0 | 0.00139,0.00148 | 0.027 | 1.6759,1.6460,1.4782 |

Relative to each preceding zero hold, U advanced +0.00320 rad for the U
reference, while V advanced only +0.00015 rad for the V reference. Absolute
errors are below both earlier recordings, but the V step was not independently
tracked. Feedback pressure maintained a 1.6000 bar common mean. Per-channel
command ranges during 55–220 s were [0.00920,0.01147,0.00857] bar, maximum
slew was [0.00803,0.01223,0.00709] bar/s, and pressureLimited was false for
all 3300 feedback samples. Overall maximum pressure rate, including the
start/end ramp, was 0.8 bar/s. NDI tip XYZ per-sample changes during feedback
had median 0.035 mm, 99th percentile 0.100 mm and maximum 0.239 mm. No tip
jump or command spike is evident in this recorded interval; distributed local
buckling cannot be assessed from tip sensing.

The near-straight interpretation is limited: thetaObservable=false during
3298/3300 feedback samples, yet inverseModelValid=true. Median feedback phi
was 0.00270 rad; its maximum 0.00897 rad was near the calibrated 0.5-degree
(0.00873 rad) observability threshold. During all five steady holds in the
table, thetaObservable was false and theta was held at 0.8166 rad by the NDI
orientation estimator. The logged measured lengths are inferred from phi and
this held theta, so U/V cannot resolve a change in bend direction there.
The stronger gains produced smooth bounded commands near straight in this
single run; they do not establish directional tracking or general stability
at the singular pose. Source: `ndiSupportFnc.orientationFromPose` and
`AuroraReader_s2_Position_only` for the held-angle behavior.

## Common-pressure circle recordings — 2026-10-07

Read-only analysis of `common_pressure_testing/FB_circle_phi_0.2rad.mat` and
`FB_circle_phi_1rad.mat` (both 240 s). The second filename is misleading:
its logged desired phi is **0.1 rad**, not 1 rad. The first run ramped phi
over 5–20 s; the second over 5–10 s. Both held theta=0 through 35 s, then
requested one smooth revolution over 35–155 s and a phi return over
155–175 s. The saved `FBL_commonP_circle.slx` inspected after the runs has
one-input `thetaCircle(t)` with `des_phi=0.1` internally; it does not yet
contain the later proposed `des_phi`/`circle_time_s` input interface. This
saved source should not be assumed to be the exact source for the earlier
0.2 rad recording. The user subsequently chose internal editable phi,
per-lap duration and circle-count variables in the one-input reference;
that revised code is also not installed in the saved model.

| Recording / logged target | Mean measured phi, 30–35 s | Desired l2,l3 (mm) | Measured l2,l3 (mm) | Mean controller P1,P2,P3 (bar) |
| --- | ---: | --- | --- | --- |
| 0.2 rad | 0.1213 rad | 1.30, 1.30 | 1.1625, 0.3410 | 1.53485, 1.72173, 1.54342 |
| Named 1 rad; logged 0.1 rad | 0.0613 rad | 0.65, 0.65 | 0.6463, 0.0808 | 1.60849, 1.68573, 1.50578 |

The corresponding l3 errors are 0.9590 and 0.5692 mm, whereas l2 errors
are only 0.1375 and 0.0037 mm. Across the circle window (35–155 s), mean
measured phi is 0.1156 rad versus 0.2 desired, and 0.0509 rad versus 0.1
desired. Logged theta and tip XYZ show substantial rotation/motion, so the
failure is bend-magnitude/directional tracking rather than no motion.

All three logged post-Kill NI-path pressure-equivalent signals equal
`pressureCommand_bar` at every recorded sample in both files (maximum
difference zero). `inverseModelValid` and NDI `orientationValid` are true
throughout the 10–220 s analysis window; `pressureLimited` is false there.
At 32 s in the 0.2 rad run, logged acceleration correction is about
[1.04,1.72] N and gravity/stiffness force about [1.88,1.30] N. The
controller therefore supplies nonzero correction while a persistent error
remains. The result is consistent with finite PD error under static
model/actuator mismatch, particularly in l3; it does not identify the
specific cause. Post-Kill NI-path agreement verifies the software command
route, not chamber pressure, which was not sensed. No model, controller,
hardware or analysis script was changed for this inspection.

The later `analyzeCommonPressureCircle.m` utility was run offline against
both MAT files and reproduced these hold values, the logged 0.1 rad target
in the second file, and zero controller-to-NI-path command discrepancy. It
prints a per-lap report for future recordings; no hardware/model validation
is implied by that offline check.
The expanded report also reads sparse Kp/Kd change-event samples correctly:
both existing runs report Kp=[1400,1400], Kd=[120,120]. It prints separate
validity counts, pressure decomposition and force summaries. This was
checked against both files without modifying them.

## One-radian circle trials 3 and 4 — 2026-10-07

`common_pressure_testing/FB_circle_phi_1rad_trial3.mat` logs one 25.01 s
circle over about 9.99–35 s; `FB_circle_phi_1rad_trial4.mat` logs five 30 s
circles over 10–160 s. Both actually request phi=1.0000 rad and log
Kp=[3000,3000], Kd=[120,120]. The 5–10 s window previously labelled
"pre-circle hold" by the report is a phi **ramp** (mean desired phi 0.495),
not a settled fixed-theta 1 rad hold. Trial 4's 20 s phi return starts at
160 s and finishes before the 220 s feedback cutoff.

| Recording | Valid circle samples | Mean measured phi (rad) | Phi MAE (rad) | Theta wrapped MAE (rad) | Measured turns | `pressureLimited` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Trial 3, one circle | 499/500 | 0.7693 | 0.2379 | 0.3026 | 0.976 | 468/500 |
| Trial 4, five circles | 2990/3000 | 0.7425 | 0.2575 | 0.2655 | 4.977 | 2894/3000 |

Trial 4 mean phi by lap is 0.7540, 0.7434, 0.7407, 0.7386, 0.7360 rad:
the persistent bend-magnitude deficit is repeatable across laps with a small
downward drift. Orientation-binned means range from about 0.926 rad at
desired theta 45–90 degrees to 0.566 rad at 180–225 degrees. This supports
direction-dependent tracking; it does not by itself identify gravity,
actuator response, or stiffness as the cause. In trial 4, l2/l3 circle RMSE
is 2.6058/3.5583 mm. Actual phi reaches up to 1.0099 rad in trial 4 and
1.1161 rad in trial 3 at some instants, so the recorded deficit is not an
absolute inability to attain 1 rad at every orientation.

Read-only inspection of raw MAT traces shows at least one channel at the
1.2 bar lower rail in 463/500 trial-3 and 2881/3000 trial-4 circle samples;
no logged channel reaches the 3 bar upper rail. Trial-4 lower-rail counts by
P1/P2/P3 are 1422/675/811. Thus lower-rail clipping explains nearly all
`pressureLimited` samples, although that flag also covers correction clamps
and slew. The saved inverse source currently has +/-1 bar s/a correction
bounds and 2 bar/s slew. Reconstructing its unclipped target from logged
force and initial static anchor matches every valid trial-4 circle pressure
command to numerical precision. On that reconstruction, symmetric `ds`
hits its +/-1 bar clamp in 1835/2990 valid trial-4 samples; antisymmetric
`da` never hits its clamp. The raw target falls below 1.2 bar in 2873/2990
valid samples, but never exceeds 3 bar. A uniform upward pressure shift
could preserve that raw differential within 1.2–3 bar in 2583/2990 valid
samples; this is an *offline allocation feasibility calculation*, not a
validated common-pressure setpoint or physical-response prediction. In the
remaining samples, the requested pressure span can exceed the 1.8 bar
available window (maximum reconstructed span about 2.10 bar).

Trial 4 has ten isolated invalid NDI/model frames during the circles and
trial 3 has one; logged controller pressure holds through inspected invalid
frames. Controller and post-Kill NI-path pressure-equivalent logs match
exactly in both runs. Actual chamber pressure was not measured, and the
known local buckle is not captured by tip-only NDI sensing. No Simulink,
controller or analysis code was modified for this result inspection.

## Five-second one-radian circle trial 6 — 2026-10-07

`FB_circle_phi_1rad_trial6_good_5s_0.5bar.mat` logs five 5 s circles from
10–35 s, desired phi=1 rad, Kp=[3000,3000], Kd=[22,22]. The user reports
lowering the active command floor to 0.5 bar. Valid/model/NDI orientation
samples are 489/500. Mean measured phi=0.9883 rad, phi MAE=0.0485 rad,
range 0.8864–1.1545 rad, and five-lap mean phi=0.9918, 0.9902, 0.9873,
0.9865, 0.9854 rad. Measured theta turns=4.989, but wrapped theta
MAE=0.3173 rad; l2/l3 RMSE=2.7649/3.0615 mm. The phi improvement relative
to 30 s trial 4 is substantial, while directional and individual-length
tracking are still imperfect.

`pressureLimited` is true in 65/500 circle samples. Per-channel pressure
ranges are [0.8057,0.7521,0.6052] to [2.6656,2.5299,2.2885] bar.
At least one command is below 0.8 bar in 93/500 samples (P2 20, P3 73),
but none reaches 0.5 bar. The observed per-channel command slew reaches
2.0 bar/s. Controller and post-Kill NI-path commands match exactly;
chamber pressure is not sensed. This run also changes circle time and Kd,
so it does not isolate the effect of the lower bound. The analyzer's
5–10 s "pre-circle hold" window is a phi ramp, not a steady hold.

Retrospective lag analysis of trial 6 finds a 0.24 s best shift between
desired and measured theta; per-lap fits are 0.25/0.24/0.23/0.23/0.23 s.
Comparing desired theta at `t-0.24 s` lowers wrapped theta MAE from
0.3173 to 0.0692 rad. The analogous offline length comparison lowers
l2/l3 RMSE from 2.765/3.062 to 0.820/1.110 mm. This quantifies an
effective phase offset, not a proven physical transport delay or a
modified controller result. Commands reach approximately 2.0 bar/s in
65/499 adjacent circle intervals.

The later `FB_circle_phi_1rad_trial6_good_10s_0.5bar.mat` logs five 10 s
circles at phi target 1 rad, Kp=3000, Kd=22. It has 990/1000 valid circle
samples, measured mean phi 0.9634 rad, phi MAE 0.0485 rad, theta MAE
0.2713 rad, 4.978 measured turns, l2/l3 RMSE 2.3119/2.5840 mm, and
pressureLimited 5/1000. A pure time shift fits 0.40 s here versus 0.24 s
at 5 s/circle. An offline linear fit of theta error against desired
angular speed gives intercept/slope 0.196 rad/0.121 s for the 10 s run,
versus 0.150 rad/0.133 s for the 5 s run. Thus the signed theta error
contains both an approximately speed-independent offset and a speed-related
part; no unique physical cause is identified. The saved model's current
`des_phi=1.5` setting was written after this 1 rad MAT and should not be
mistaken for its logged target.

## Documentation migration checks — 2026-10-06

Historical log copied before editing with matching SHA-256
`31A09742DA0F9EDEF769C9733C263897B4CA066C44B0778180429C8B5A8C1EA7`.
Read-only FBL_commonP archive inspection finds both supplied feedback scripts;
this updates the earlier old-script observation, without verifying wiring,
compilation or recording equivalence. No experiment results were regenerated.

## Dynamic Ki switch during 1.5 rad circles — 2026-10-08

`StaticTest/changeKi_dynamically.mat` logs five 10 s circles at desired
phi=1.5 rad with Kp=[3000,3000] and Kd=[42,42]. Ki changes from
[700,700] to [0,0] at 43.31 s, inside lap 3. The base remains nearly
vertical throughout. For a phase-matched comparison, complete laps 1-2
(20-40 s, Ki on) and laps 4-5 (50-70 s, Ki off) were used; valid NDI/model
counts were 389/400 and 395/400. Lengths are orientation-inferred. The
angle error is desired minus measured, with theta wrapped to [-pi,pi].

| Metric | Ki 700, 20-40 s | Ki 0, 50-70 s |
| --- | ---: | ---: |
| Phi signed error | +0.010 rad | +0.089 rad |
| Phi MAE | 0.0917 rad | 0.1080 rad |
| Phi error standard deviation | 0.111 rad | 0.083 rad |
| Wrapped theta MAE | 0.2208 rad | 0.2194 rad |
| l2 MAE | 2.391 mm | 2.391 mm |
| l3 MAE | 3.030 mm | 2.976 mm |
| l3 signed error | +0.478 mm | +1.002 mm |
| `pressureLimited` | 1623/2000 (81%) | 115/2000 (6%) |

Ki=700 substantially reduces persistent phi underbend but does not
materially improve theta or absolute length tracking. Phi variation and
pressure limiting increase with Ki on. Nominal pressure rails are touched
for about 14.5% of command samples with Ki on and 5% with Ki off; the
limited flag also includes internal correction and slew limits, so it is
not equivalent to rail saturation. The integral state was not logged. The
post-circle zero-reference hold occurs only with Ki off and is not a
matched comparison. This single on-to-off trial does not establish a
general stability margin or isolate slow drift. No model or hardware was
run by Codex.

## Integral-run fault at 123 s — 2026-10-08

`integral_term/plot_data_110s.mat` records through 140.05 s and identifies
model `FBL_commonP_circle_BRot`. The desired phi remains 1.5 rad and the
theta reference continues past 123 s. The tip 0B validity briefly drops
at 123.00-123.05 s, recovers, then drops again for 123.30-123.70 s;
base 0A and base-orientation validity remain good. InverseModelValid is
last true at 123.29 s. The controller holds [1.891,0.646,2.262] bar until
123.39 s, begins a 2 bar/s ramp toward zero at 123.40 s, reaches all-zero
at 124.50 s, and remains latched after NDI validity returns at 123.75 s.
Controller and NI-path command logs match (max difference zero); the
logged Kill state is 1. This validates the sustained-invalid command-path
latch and recovery lockout for this user-recorded run. It does not measure
actual chamber pressure or identify why tip tracking failed.

The saved SLX StopTime is 250 s, while the inverse and integral controller
time bounds are 220 s. No 120 s MATLAB Function time gate was found;
the apparent stop at 123 s was the fault response.

## Two-run Ki angle comparison script — 2026-10-08

`integral_term/compareKiAngleErrors.m` was validated offline with
`without_i_1rad.mat` and `withi_1rad.mat`. It generated the requested
theta, phi, and angle-error figures on their shared 0-80 s interval.
The logged desired references match. Kp differs (1800 versus 3000),
while Ki differs (0 versus 300) and Kd is 42 in both runs, so these
recordings cannot isolate the effect of Ki.

## 2026-10-08 - Nominal passivity controller numerical verification

`passivity_adaptive_control/validatePacFixed.m` passed MATLAB R2025a
hardware-free checks: nine-point +/-8 mm inertia positivity/symmetry,
Mdot-2C skew symmetry (1.9626e-17 maximum residual), finite-difference dM
agreement (2.5032e-10 relative error), force summation, matched-plant
closed-loop/energy identity (1.7101e-15 N balance residual), invalid-input
and base-validity gating, and gravity-direction response. Four-second
matched-model ode15s tracking ended with 3.8647e-7 m error norm.
These checks do not cover model mismatch, pressure limits/deadzone/slew,
filtered/discrete measurement dynamics, hardware, or SLX/codegen compilation.

## 2026-10-08 - PAC pressure allocation and fault supervisor

`validatePacPressure` passed MATLAB R2025a offline tests of trim-force bias,
all-active unconstrained allocation (maximum residual 2.2888e-15 N), pressure
rails/correction clipping, nominal activation crossing, 2.0 bar/s slew,
stationary-anchor condition, changed gravity without anchor recapture,
brief-invalid hold/recovery, sustained-invalid/stale-data latch, no spontaneous
fault recovery, disable/re-enable, Clock faults and combined-wrapper execution.
Stage-1 mechanical checks also passed unchanged. Code Analyzer reported zero
messages. MATLAB Coder C source generation passed for fixed input dimensions
and logical validity/enable flags. No binary build, SLX/model execution,
hardware access, chamber-pressure measurement or sampled closed-loop physical
validation was performed. Nominal below-threshold force estimates and the
constant trim calibration require experimental evaluation.


## 2026-10-08 - PAC.slx saved wiring review (no execution)

Read-only ZIP/XML inspection of the user-created PAC.slx confirms measured
[q;dq] (4), desired [qd;dqd;ddqd] (6), gravity_arm, orientationValid_control,
baseOrientationValid, conditional measurement timestamp, Clock and boolean
true enable are connected to the eight pacCommonPressureBlock inputs in order.
The sensor uses AuroraReader_s2_BaseRotation, COM12/921600, 0.05 s, the existing
20260930 calibration (file exists), and reference gravity [0;0;-9.81].
Pressure connects directly to the retained Kill switches, 0-3 bar saturation,
10/9 voltage conversion and NI outputs. Kill has no saved Value override
(default Constant value 1, which selects pressure through these switches).

Outstanding integration findings: controller and validMeasurementTime have no
saved explicit discrete rate and receive continuous Clock; configure both
Discrete 0.05 s to avoid continuous persistent-state execution and false
timestamp freshness between polls. Solver is VariableStepAuto. Active circle
reference is phi=1.5 rad, 20 laps at 5 s/lap; phi ramp begins at 15 s and
circles at 20 s, overlapping the earliest 15-20 s feedback blend. Hold the
reference through trim/engagement (status=3) before starting its ramp.
The XYZ scope desired input has no source. Controller diagnostics are named
and terminated, but no DataLogging enable entries were found in saved system
XML; naming and global logsout enable do not select individual signals.
Enable signal logging explicitly for controller outputs, inputs, validity,
gravity/timestamp and post-Kill pressures. Ensure root and isolated PAC folder
are on MATLAB path. No SLX compile/run or hardware was attempted; no model
or controller source changed. This remains the fixed-parameter passivity
baseline, with adaptation pending.


## 2026-10-09 - First PAC circle recording: feedback authority too low

Offline MATLAB read of PAC_circle/1rad_10s_1cycle.mat (78 logged signals,
0-61.1 s). Supervisor trim starts 10 s, anchor captures 15 s, blend finishes
20 s; one invalid hold 30.60-30.65 s, no latched fault. Circle 25-35 s:
phiOrientation mean 0.016283 rad (range 0.011401-0.019246), desired phi 1 rad.
All pressure commands remain above activation; limited false throughout circle.
Controller pressure equals post-Kill NI-path commands; ranges P1 1.6614-1.6794,
P2 1.6306-1.6817, P3 1.4452-1.4974 bar. Valid model 199/200 circle samples.
Feedback force components bounded approximately +/-0.295 N; inertial terms
approximately +/-0.028 N. Valid allocation residuals are roundoff (~4e-15 N).
Therefore no evidence of command gating or saturation as cause of poor bending.

Source nominal gains lambda=2 /s and Ks=10 N s/m give restoring coefficient
Ks*lambda=20 N/m (plus M*lambda velocity feedback), much weaker than nominal
650*[2 1;1 2] N/m stiffness. At qd=[6.5;6.5] mm, q near zero, static feedback
is about [0.13;0.13] N versus nominal stiffness force [12.675;12.675] N at the
desired pose. The compensation uses measured q, so near straight it supplies
little stiffness compensation; finite force-feedback gains must drive initial
error reduction. Initial gains were offline matched-model starting values, not
experiment-tuned gains. Increase/tune restoring authority with damping and
allocation limits considered; adaptation is still absent. No source/model edits
or hardware/model execution were performed.


## 2026-10-09 - Reusable PAC static analysis and first static test

Added PAC_circle/analyzePacStatic.m and local README at user request. Offline
only: automatic longest nonzero constant-reference hold, latter-half settled
metrics, inferred Ks/lambda from recorded forces/sliding signal, engagement/
fault/validity/pressure rails/freshness, post-Kill command comparison, force
components/residual, tracking plots, text report and compact analysis MAT.
Explicit hold window handles repeat trials. Gains are identified under settled
desired-rate assumptions rather than read from current controller source.

Validated on static_phi_test1.mat: 0-99.25 s; plateau 30.05-60 s; settled
45.025-60 s. Desired/measured phi .2/.01539 rad; inferred Ks=[10 10],
lambda=[2 2]; feedback means [.02254 .02947] N. Mean length errors
[1.1277 1.4707] mm. Anchor 15 s, full control 20 s; no latched fault or
limiting in hold; one invalid hold sample. NI-path command difference zero;
residual roundoff, timestamp age 0-.05 s. Weak restoring authority persists.
No controller or SLX edits; no hardware/model execution.


## 2026-10-11 - PAC static gain-pair comparison (user-supplied reports)

Compared lambda=400,Ks=10 with lambda=200,Ks=20; both Ks*lambda=4000 N/m.
Reported settled phi .15444 vs .15214 rad, errors .04556 vs .04786 rad.
Length errors [.21436 .38031] vs [.23269 .39172] mm. New run has no hold
limiting/fault/invalidity and controller-to-NI command difference zero.
Settled inertial-force peaks fall [.11041 .085709] to [.021285 .016562] N
(about fivefold); lower measured motion also contributes, so this is not a
controlled proof of gain effect. Gravity and captured forceBias differ between
runs. Static error remains near 24 percent, supporting a persistent mismatch
rather than allocation clipping. Keep 200/20 as provisional baseline; repeat
static holds/directions to check consistency before adaptive compensation.
Evidence is pasted reports/screenshots, not independently reloaded raw MAT.
No controller or SLX edits.


## 2026-10-11 - Direction-dependent PAC static performance

User-supplied reports/screenshots for static_phi_lambda_200_damp_20_theta_Piby2.mat
and static_phi_lambda_200_damp_20_theta_Pi.mat; raw MAT not independently
reloaded for this comparison. All use phi=.2 rad, lambda=200 /s, Ks=20 N s/m,
Ks*lambda=4000 N/m, hold 30.05-60 s and settled latter half.

| Desired theta | Settled phi | Phi MAE rad | Theta wrapped MAE rad | Mean l2/l3 error mm |
| --- | --- | --- | --- | --- |
| 0 (previous run) | .15214 | .04786 | .04643 | .23269 / .39172 |
| pi/2 | .19065 | .00935 | .16734 | -.34171 / -.071021 |
| pi | .16713 | .03287 | .03166 | -.15466 / -.27375 |

Both new holds have full engagement/anchor, 600/600 valid samples, zero limited
samples, all commands above activation and zero controller/NI command mismatch.
Gravity is similar to theta=0 trial. pi/2 feedback mean [-1.3668 -.28411] N;
pi feedback mean [-.61858 -1.0949] N. Static allocation residuals are roundoff;
pi/2 ramp/return shows oscillations and a transient residual outside the hold.
Good pi/2 phi alone hides a ~9.59 degree direction error and unequal length
errors. pi theta error ~1.81 degrees. Captured forceBias differs by trial.

Interpretation: persistent mismatch is direction dependent, not explained by
hold pressure clipping. Mechanics/stiffness/coupling, pressure-force mapping,
trim bias and orientation geometry remain candidate causes; no unique stiffness
identification from these tests. Keep 200/20 provisional baseline; formulate
adaptation in both length coordinates. Distinguish additive force offset from
stiffness/coupling mismatch using repeats/multiple amplitudes before selecting
adaptive parameters. No adaptation implementation authorized by this discussion.

Prior gain discussion: lambda converts length error to reference velocity (1/s),
Ks converts sliding velocity error to force (N s/m); direct feedback has
position coefficient Ks*lambda and velocity coefficient Ks. Exposing gains as
2x1 block inputs is possible but was discussed only, not implemented by Codex.
Constant gains per run preferred initially; abrupt lambda changes alter
reference dynamics and constant-gain stability arguments do not automatically
apply. Physical params.D remains distinct from controller Ks.


## 2026-10-11 - Combined PAC static sweep implementation

Added PAC_circle/pacStaticSweepReferenceBlock.m: 9 direction/amplitude trials
plus theta=0/.2 repeat, start25 s, 5/30/5/10 s rise/hold/fall/zero-rest, end525 s,
Stop Time535 s. Readiness must reach status3 before start. Zero-reference rest
does not vent pressure or recapture trim. User pastes block; no SLX edits.
Added analyzePacStaticSweep.m: detects >=10 s logged constant nonzero holds,
reuses per-hold static analyzer without intermediate artifact overwrite, prints
compact tracking/gains/force/pressure/quality reports and saves summary CSV, TXT,
MAT. Exploratory unconstrained affine residual fit is not physical stiffness
identification. analyzePacStatic adds optional fourth writeArtifacts flag.
Validated analyzer on existing single-hold recording and all ten reference
plateau IDs/amplitudes in MATLAB offline. Full multi-hold physical validation
is pending; no hardware or Simulink run performed.


## 2026-10-11 - Ten-hold PAC sweep: repeat changes and affine fit limits

Evidence: user-pasted static_phi_lambda_200_damp_20_sweep.mat sweep report and
summary screenshot; raw MAT not reloaded. All settled holds valid/engaged, no
limiting, above activation; lambda200/Ks20. Theta0 phiActual .07259/.14517/.22800
for desired .1/.2/.3; pi/2 .09539/.18283/.27894; pi .04626/.12473/.21352.
Same theta0/.2 repeated at hold10 gives .10945 vs .14517 at hold2 (24.6%
decrease); feedback [2.2642 2.4454] vs [1.0893 1.7723] N, despite fixed
forceBias [-.69562 -2.9912] N. Changes demonstrate time/history dependence
or another changing condition; cannot isolate hysteresis, creep, air supply,
orientation drift or temperature from report. Direction and trial order confound.
pi/2 theta MAE .54541/.27384/.20097 rad at increasing phi, so phi magnitude
alone overstates accuracy. pi/.1 theta omitted by observable-angle threshold.

Affine fit deltaK=[1159.7 -22.415;1241 987.44] N/m, b=[-.045984;-.45184] N,
training RMSE [.68587 .65032] N, predictor condition2.338. Full rank does not
remove history confounding; strongly nonsymmetric fit is not a conservative
stiffness correction. Do not insert into EoM. Next repeat same theta0/.2
several times with unchanged anchor and track zero-reference baselines; separate
repeatability/history from fixed configuration dependence before selecting
stiffness versus adaptive bias parameters. No code/model changes.


## 2026-10-11 - Adaptive PAC force/pressure path implemented

User authorized implementation in passivity_adaptive_control only. Added
pacAdaptiveConfig/InitialState/History/Force/Update/CommonPressureBlock,
validatePacAdaptive, validatePacAdaptiveTracking, replayPacAdaptive,
analyzePacAdaptive and ADAPTIVE_INTEGRATION.md. Existing nominal sources/SLX
files preserved. New wrapper has previous 8 inputs plus 2x1 lambda, 2x1 Ks,
boolean enableAdaptation; previous17 outputs plus Khat,Dhat,hysteresis force,
bias,10parameter estimates,3history states,active/bound flags,adaptive force.

Khat is spring-structured SPD: [k1+k3 k3;k3 k2+k3], positive independent
springs/nonnegative shared coupling; Dhat diagonal positive. Parameter bounds,
rate caps/deadband and scale-adjusted gradient in config. Three normalized
Bouc-Wen feature states use physical muscle velocities T*dq, with fixed internal
coefficients. Their nonnegative force weights and two signed bias coefficients
adapt. Does NOT adapt alpha/beta/gamma internally, nor M/C/G. Paper online
adaptation is K/D, not inertia/hysteresis-internal coefficients.

Adaptive controller uses Dhat*v instead of old D*dq, giving a changed nominal
error response (initial (Dhat+Ks)*lambda=12000 N/m at lambda200/Ks20).
Adaptive on/off comparison must use same new wrapper. Updates apply NEXT
sample, once per fresh timestamp, only valid status3/no limiting/residual
<=.02 N; dt gaps>.075 s skip updates. Freeze retains estimates/history;
restart simulation resets. Physical Kill and pressure supervisor unchanged.
Projection bounds parameters, not physical closed-loop stability. Fixed-feature
history, saturation, sampling/rate caps/deadband mean ideal asymptotic proof
cannot be claimed for hardware. Static damping/bias/pressure mapping confound.

Verification: algebra/regressor/force sums, SPD projections, descent sign,
ideal unprojected energy identity, history reversal, startup/fresh-only update,
invalid hold/fault checks pass in MATLAB. Ideal continuous force mismatch
simulation at lambda200/Ks20 (30s) improves final10s RMSE from
[.10960 .12827] mm to [.05988 .07452] mm. This is a chosen representable
plant, not physical validation. Replay of existing 10701-sample sweep finite,
7633 update samples, zero fault samples, pressure maxima [3 2.1240 2.3221]
bar; record-fixed replay is not closed-loop prediction. C source generation
passed with fixed dimensions/logical flags; no binary/SLX build or hardware
run. Initial generation commands failed due literal temporary output path;
corrected functional codegen call passed and generated workspace dirs removed.
Parameter estimates are bounded effective coefficients, not unique physical
identification. Integration/model wiring and experiment checks remain pending.

Adaptive analyzer report and plots additionally verified with synthetic
Simulink signal Dataset logs generated by offline wrapper calls.


## 2026-10-11 - First adaptive-wrapper physical tests: oscillatory baseline

Read three saved PAC_phi_0.2 recordings offline (statis, enable_inMiddle,
enable_fromBegining). During45-60s desired.2, measured phi means
.4226/.4300/.4298 rad, ranges .058-.694/.033-.696/.052-.696. Limited
fractions .897/.934/.910. Dhat*v force peaks near60N, feedback near30N,
request peaks near85N; allocation residual peaks near90N. Learning disabled
already oscillates, so adaptive gain tuning cannot explain baseline behavior.
Middle enable yields zero updates; beginning enable yields112 updates mostly
20-26s and tiny changes, then frozen. Gates intentionally inhibit learning
when limited. Limiting includes correction clipping/slew; pressure below3bar
does not mean force request realizable.

Important implementation/tuning correction: Dhat*v adds Dhat*lambda position
feedback, giving initial12000N/m at lambda200/Ks20, vs old direct4000N/m.
Offline continuous-force simulation did not validate this setting against
50ms sensor/filter/pressure delays. Nominal gain carryover unsuitable.
Proposed next new-wrapper fixed-gain baseline lambda50/Ks40 gives initial
(Dhat+Ks)*lambda=4000N/m and less lambda acceleration amplification.
Check learning-disabled response before adaptation; do not remove freeze gate
or increase adaptation rates/pressure limits to address oscillatory baseline.
No code/SLX edits made for this diagnosis.


## 2026-10-12 - Reduced-gain adaptive comparison

Offline read PAC_phi_0.2_static_reducedGains_engagedMiddle.mat and engagedFalse.mat.
Gains lambda50/Ks40. During45-60s, desired phi.2, adaptive/fixed actual
.16658/.15206 rad; MAE .03342/.04794 (~30.3% reduction). Mean length errors
[.14369 .29258]/[.22465 .40123] mm (~36%/~27% reduction). Feedback means
[.28708 .58477]/[.44933 .80235] N. No hold limiting in either; full run16
limited samples each, zero faults. Adaptive phi range .15499-.17547 versus
fixed .14910-.15659; learning continues, not fully converged. User reports751
updates; final spring coefficients [677.45 699.57 810.52], damping40.012/40.024,
small history weights, bias[-.00927 .12572] after return. Final values after
return differ from hold values; not physical parameter identification.
Keep gains/rates unchanged initially; longer hold or full sweep now useful
for directional/history tracking evaluation. No code/SLX changes.


## 2026-10-12 - Adaptive 10-hold sweep: good .3 rad, limited transfer

User supplied PAC_sweep.mat reports/screenshots, 10701 samples0-535s,
lambda50/Ks40,9601 adaptive updates,20 limited samples total,no faults.
All reported hold latter halves valid/engaged/no limiting/above activation.
At desiredphi.3 theta0/pi2/pi actual .30369/.29674/.29744; phiMAE
.00445/.00326/.00497 rad; thetaMAE .01817/.02628/.01833 rad. At .2 first
theta0/pi2/pi actual .18019/.19137/.18679. Finaltheta0/.2 actual .13643
vs initial .18019 (24.3% decrease), error .06357 vs .01980. This is
limited transfer/readaptation/history performance, not global convergence.
At pi2/.1 thetaMAE .46708rad despite good phi. Effective spring/history/bias
estimates bounded; D estimates only40.232/40.336 versus initial40, so no
evidence of uniquely identified physical damping. Final parameters are
compensation estimates; all nominal-baseline comparisons across gain/form
changes are confounded. Proposed application-relevant next amplitude .5rad
with unchanged gains/rates,10s ramp,60s hold,include return to earlier .3/.2
to evaluate transfer. Keep constraints and observe direction/length errors.
No controller/SLX/source changes.

## 2026-10-12 - High-phi adaptive sweep: unloading onset and slew-limited oscillation

Offline raw-data review of passivity_adaptive_control/PAC_sweep_high_phis_unstable.mat,
10701 samples, 0-535 s. Reference amplitudes .5/.8/1.5 rad. The onset at110 s
is the theta0/.8 unloading ramp110-115, not the 1.5 command (starts125).
During95-109 phi .804-.812 (mean.808), all orientations valid, no limiting,
force allocation residual roundoff. During110-115 adaptation active94%;
spring estimates [568.1 804.9 824.8] at110 become [467.7 682.6 397.9] at115;
k3 falls52%. K remains positive definite (eigenvalues678/2344 to561/1385 N/m).
First limited113.950 s, first invalid orientation115.700 s, first below
activation128.500 s. During115-125 desiredphi0, actual .0186-.462 rad,
SD.1008 rad, approximate dominant frequency1.9 Hz; limited and 2bar/s slew
hit86.5%, valid92%, adaptation active5.5%. Pressures minima[1.637 1.037 .9104]
bar: above .8 activation during this first oscillatory zero interval.
Feedback force swings about -12.7..11.9 N; requests -20.2..19 N and allocation
residual about +/-20 N. Before onset feedback magnitude was below.311 N.
Later pi2/.5 hold195-210 meanphi.515 conceals range.114-.889, limited94%.
Estimates remain bounded, damping changes small, zero latched faults. Later
activation residual can freeze learning even limited=false; these flags are
not equivalent. User-reported1537 updates/7232 limited samples.
Leading hypothesis: unloading triggers delayed force-feedback/pressure dynamics
interacting with slew limiting; rapid transient adaptation may contribute.
Timing excludes invalid orientation or activation crossing as initial trigger,
but does not prove adaptation causal or distinguish pneumatic/sensor delays.
Next diagnostic priority is transient plots/time-window metrics and a matched
frozen-estimate unloading comparison; consider slower10-20 s returns and learning
gating during return/recovery after review. Do not increase slew or remove gates
on this evidence alone. No controller or SLX changes, no hardware run.

## 2026-10-12 - High-phi repeat trial confirms unloading failure

Offline comparison PAC_sweep_high_phis_unstable_trial2.mat (0-193.7 s)
with trial1. Before110, trial2 phi .80527-.81167, SD.00102, no limiting.
Return starts110 s; first limited113.800 versus113.950 trial1; first
invalid orientation114.700 versus115.700; below activation128.300 versus128.500.
Oscillation grows approximately113.6-115 s; no mathematically unique onset
sample is asserted. During115-125 phi .01599-.45883, SD.09137, limited88%,
valid90.5%, residual peak20.347 N; trial1 corresponding phi .01865-.46166,
SD.10079, limited86.5%. Strong repeatability ties failure to .8 unloading,
not random startup or initial deadzone crossing. Slew limiting remains a
constraint/possible amplifier, not proven original cause. Adaptation-versus
sampled feedback/pneumatic delay causality remains unresolved. Hold4 meanphi
.5248 conceals MAE.2455. No controller/SLX changes or hardware execution.

## 2026-10-12 - Freeze-adaptation recordings
Offline changing_enabling_adaptive_tr1/tr2.mat reviewed. Both adapt: last active
104.75/102.90 s, last parameter change104.76/102.91 s. Trial1 is NOT an
adaptation-off-from-start baseline (1106 updates and changed parameters).
Trial2 quiet105-110 phi .8077-.8099, no limiting; frozen parameters then
unloading110-115 followed115-125 phi .0606-.4638, SD.0848, limited87%.
Thus continued parameter updates during unloading are not necessary for failure.
This does not exclude learned model/history dynamics as contributors.
Both logs .01 s vs earlier .05, measurement timestamp unique increments.05 s.
Possible multirate supervisor/slew effect: initial .5 hold limited~62-65% but
residual only~.005 N, so raw limited counts are not comparable across rates.
Need verify controller execution/sample holds before retuning. Trial1 already
oscillates95-105 phi .8324-1.045, whereas trial2 quiet. Trial2 1.5 hold mean
near desired conceals135-150 range1.438-1.625 and88% limiting. No source/SLX edits.

## 2026-10-12 - Next isolation experiment and verified constraints
Source review: pressure rails .2-3 bar, s/a correction +/-2.5 bar, slew2bar/s,
nominal activation .8 bar (minimum below activation), fixed mean1.6. No timeout.
Invalid grace.10s, max sample age.075s. Allocation assumes linear force above
activation; crossing activation yields residual without necessarily limited=true.
At static target tau includes Dhat*v-Ks*s = (Dhat+Ks)*lambda*error-Ks*dq.
Current D40/Ks40/lambda50 gives4000N/m restoring and40Ns/m explicit measured
velocity coefficient, plus M*a/C*v. Frozen adaptation does not freeze history z.
Hysteresis is plausible but not established. Next verify controller sample time
.05s consistent with fresh sensor data; then short identical .8 ramp/hold/return,
frozen estimates before return, compare Ks40 versus20 with lambda50 unchanged.
This reduces both feedback stiffness and velocity authority; cannot isolate pure
velocity causality. Retain limits; log rate/rail/correction reasons separately
when implementing diagnostics. No code changes in this discussion.

## 2026-10-13 - PAC .05s execution confirmed; severe wall-time slowdown
User sample-time overlay confirms PAC outputs .05s, NI path .01s. Pacing off;
user reports .05 simulation seconds per wall second. Profile attempt aborts in
AuroraReader_s2_BaseRotation.setupImpl preflight:100/100 valid paired,88 overruns,
median interval1.3434s, reason polling slower than configured period. Screenshot
PAC call count0: profiler cannot assess PAC runtime from this run. Reader source
preflight occurs before stepImpl; physical polling paced via tic/toc in stepImpl.
Thus this failure establishes slow acquisition/setup under profiling, not PAC
computation cost. Profiler overhead versus serial/host delay unresolved. Next
check preflight without profiler, preferably standalone zero-pressure acquisition;
retain .05 period/thresholds and do not bypass preflight. No code/model edits.

## 2026-10-13 - Restart resolves NDI slowdown; new frozen high-phi onset
User restarted MATLAB and acquisition now works; PAC inherited rate retained.
Offline static_phi_test1_unstable.mat: last parameter change102.710s, despite
reported disable110s (updates may naturally gate earlier).105-110 phi.8088-.8108
quiet no limiting.115-125 desiredzero actualphi.1398-.2195,SD.0185,limited2.5%,
valid100%: residual bend rather than previous large return oscillation.
Next1.5 ramp125-130: first limiting127.600, invalid128.000, activation crossing
128.150.128-130 limited77.5%, residualpeak20.48N;130-140 phi1.3395-1.5985,
limited85.5%,allAboveActivation0. New onset tied to larger loading ramp, not
continued parameter updates. Deadzone-map/allocation interaction plausible,
not proven initial cause (limiting/invalid precede crossing). Theta nearphi0
ill-conditioned; length/component tracking required to establish true lateral
oscillation. Prior pure unloading mechanism not universal across runs. Next
isolate amplitude at .8 vs1.0 with identical frozen setup and rates, avoiding
simultaneous gain/ramp changes; exact source limits remain .2-3bar/slew2.
No code/model changes.

## 2026-10-13 - Commit preparation and offline validation
User requested committing pending updates on codex/passivity-adaptive-controller.
Included PAC implementation/models, experiment recordings/reports, plotting updates
and project documentation. validatePacFixed, validatePacPressure, validatePacAdaptive
and validatePacAdaptiveTracking passed offline. No hardware models executed; these
checks do not establish physical sampled-loop stability. No remote push requested.

## 2026-10-13 - 1.0-rad target does not prevent frozen-loop oscillation
Offline static_phi_test1_unstable_tr2.mat reviewed. Last parameter change99.160s.
105-110 quietphi.8054-.8078, no limiting.110-115 unload limited15%, residual
peak10.657N;115-125 desiredzero phi.0618-.4539,SD.0892,limited88.5%,valid91%.
All commands above nominal activation in these intervals (min .9215bar during
115-125). Failure already present before third target1.0 starts125, so reduced
third amplitude cannot explain onset.130-140 phi.9784-1.2321,limited86%,
allAboveActivation100%,residual17.89N. Continued updates and nominal activation
crossing are not required. Previous run had much quieter return despite similar
frozen estimates; calibration bias/history/physical state differ. Leading issue
remains sampled feedback/pressure dynamics interacting constraints, not proven
hysteresis. Next restore verified explicit .05 PAC execution after MATLAB restart
and confirm wall time; then compare Ks20 vs40,lambda50,identical reference/freeze.
No source or model changes.

## 2026-10-13 - Trial3 .05s recording: timing alone insufficient
Offline static_phi_test1_unstable_tr3.mat,3668samples0-183.35 implies.05s
logging; report lambda50/Ks40 unchanged. Last parameter change102.65s.
105-110 quietphi.8088-.8102,no limiting.115-125 desiredzero phi.1069-.4612,
SD.0866,limited85%,valid90.5%,aboveActivation100%,residualpeak18.30N.
130-140 actualphi.9731-1.2578,limited90%,aboveActivation100%. Thirdtarget1.0
holdMAE.1653 vs trial2.1603; failure persists at .05 logging, controller rate
requires compiled verification. Quiet static micro-limiting flags disappear,
but actual transient force deficits persist. Ks20 comparison still not done.
Next change Ks40 to20 only with verified.05 timing and identical reference/freeze;
not simply more adaptive updates or removing pressure constraints. No code edits.

## 2026-10-13 - Ks20 removes observed oscillation, residual tracking errors remain
Offline static_phi_test1_20damping_tr4.mat:lambda50/Ks20,.05s logs,last parameter
change103.70s. Compared tr3Ks40,115-125 phi now.2033-.2573,SD.0119 vs.0866,
no limiting/invalidity and residual roundoff vs85% limiting/18.30N residual.
Return remains biased ~.2rad for desired0; suppression is not zero tracking error.
130-140 phi.8898-.9429,SD.0103,no limiting. Report1rad holdMAE.04165 vs.16528
tr3 (~75% reduction);thetaMAE.00186 vs.19117. Only16limited fullrun, allholds0.
pi2/.5 holdthetaMAE.51148,lengtherrors[-1.9865 -1.5063]mm despitephiMAE.04893:
directional compensation remains poor with frozen parameters. Source force law
Dhat*v-Ks*s changes both velocity coefficient (40->20Ns/m) and static corrective
coefficient (Dhat+Ks)*lambda ~4000->3000N/m; improved physical response does not
prove pure damping or identify hysteresis. Supports excessive sampled feedback
and slew interaction. Hold Ks20/lambda50 as provisional baseline; repeat same
run before adaptation-on full directional sweep or faster trajectories. No code
changes; results not a physical stability proof.

## 2026-10-13 - Full adaptive Ks20 sweep: accurate holds, recovered transient
static_phi_test1_20damping_full_tr5.mat,.05 logs0-535s,lambda50/Ks20,
7051updates288limited,no faults. All latter-half holdvalid100%,limiting0;
phiMAE .00208-.02109 and thetaMAE .00177-.03823. Directionpi2/.5 now
phiMAE.00208/theta.03758 versus frozen trial4 .04893/.51148 (separate-run,
not isolated proof). Repeattheta0/.5MAE.01139 vs.00867 initial.
Raw260-265 pi2/.8 unloading,265-275zero,275-2801radloading.245-260quiet.
265-275phi.0091-.3703,SD.0892,limited79.5%,valid95%,aboveactivation100%,
residualpeak12.059N.275-280limited79%,residual15.583N.280-285limited2%,
295-310none,phi.9774-.9900. Recovery after next nonzero target is not proof
zero-reference oscillation self-decays; reference changed275s.
k3 reaches0 lower projection bound during run, does not imply physical muscle
stiffness0; K can remain SPD via k1,k2. Nominal activation residual~.045-.048N
in1rad holds exceeds.02 adaptive gate and may freeze learning without limited.
Hold gains/rates; next identical fullrepeat or isolated extendedzerohold at265
before next ramp tests persistence/recovery. No faster circles or gain increase
until transient resolved. No code/model edits.

## 2026-10-13 - Raw k3/velocity/force analysis of full tr5
k3 first reaches0 at171.550s; stays0 throughout245-310 including quiet
245-260 and recovered295-310. Thus0 is not sufficient cause or sudden260trigger.
K eigenvalues at260 [462.25 1057.2],270 [285.48 826.74],280[225.31 809.65],
295[418.33 1029.7] N/m remainpositive. k2 drops462->225 overtransition.
245-260 velocityRMS[.01229 .01460]mm/s vs265-275[6.536 9.491]mm/s;
velocitypeaks~15/21mm/s.265-275 dampingForce(Dhat*v) min[-6.66 -9.503],
max[8.212 9.127]N; feedbackmin[-3.446 -5.069]max[4.223 4.881]; inertia
peaks<.784N; staticpeaks<3.54N. Force requestmin[-9.027 -14.98],
max[7.702 11.32].Large D*v is not directlyphysicalvelocity damping:
v=dqd+lambda*(qd-q),hence zero-referenceD*v=-D*lambda*q; sensitivity to
position error (reference velocity), feedback alsousesmeasuredvelocity.
Dominant force swings are damping-compensation plusfeedback, not inertia or
lossSPD. Measurements filtered/inferred cannot distinguish actual oscillatory
motion from derivative amplification/delay alone; length/pressure phase review
needed. Negativek3 unsupported; next isolate reference/velocityforce sensitivity
via gains/lambda comparison or simulation, do not blindly alterterms andclaim
passivity. No controller/model edits.

## 2026-10-13 - Lambda25/Ks20: oscillation suppressed, adaptation settling slower
Read full pasted report and raw static_phi_test_20d_25Lmabda.mat0-512.85,
17limited6503updatesno faults; all evaluated holdsvalid100%,limiting0.
Reference changed vslambda50:rise/fall1s,42s trialspacing vs5s/50s before,
calibrationbias also different; cannot attribute every improvement solelylambda.
pi2/.8 unload224-225,zero225-235 phi.1162-.2113 SD.0158,no limiting,
compared previouszero265-275SD.0892/79.5%limiting.1radpi2 hold251-266
phi.9503-.9626,MAE.04293; allocationresidual.0407N exceeds.02 adaptationgate.
Theta0/1rad residualpeak.8635N,P1mean.7403bar implies nominal activationmapping
limits force fidelity despite limitedfalse. Repeattheta0/.5 hold404-434 actual
starts.051-.085 during404-409,419-434range.174-.508 andMAE.17078, not just
fullysettled offset: directionalreadaptation/transfer poor. Last updates508.15.
Lambda reduction lowers(D+Ks)*lambda approx3000->1500N/m; measuredvelocity
coefficientKs20 unchanged. Conservative baseline25/20; next matchedreference
comparison lambda30 or35 only, not higher adaptation gains simultaneously.
Settling report must distinguish convergence to target vsdrift/model mismatch.
No source/model changes.

## 2026-10-13 - Lambda40/Ks20 improved tradeoff, brief return transient
Read report/raw static_phi_test_20d_40Lambda.mat,.05logs0-484.65,same1s ramps/
42s trialspacing aslambda25.58limited6036updatesnofaults;all evaluatedholds
valid100%,limiting0.225-235zeroafterpi2/.8 return phi.014-.2092,SD.0213,
limited2%,valid99%,residualpeak.465N: brief transient rather than sustained
lambda50 returnoscillation.1radtheta0MAE.0150 vslambda25.0401;
pi2MAE.0253 vs.0429. Repeattheta0/.5MAE.0183 vs.1708;
419-434phi.496-.5256,SD.0059 vslambda25.174-.508,SD.108.
Fullresidualpeak2.424N vs15.583lambda50full, .8635lambda25.
Initialstates/calibrationdifferent;singletrialnotstabilitymargin proof.
Provisional workingbaseline40/20 withsame.05timing;repeatbeforefastertrajectories,
retainbounds andrates. Activationmapresidual persists1radholds(.269/.055N),
canblockadaptationdespite nolimiting. No code/model edits.

## 2026-10-13 - Adaptive circle1rad30s per lap
Raw PAC_circle/AdaptiveController/1rad_60s.mat reviewed; filename60s but logged
reference theta speed.20944rad/s confirms30s period,7complete laps20-230 and
partial8 through235.7. phi1rad. Lap1MAEphi.0453 theta.2348(lengthRMSE2.079/
2.272mm) tolap7phi.0305 theta.1954(RMSE1.888/1.800mm). Meanphi~.984;
valid>99%,laplimited~.7-2.5%,no faults perreport. Allaboveactivation only65.7-
73.9%: nominal mapallocation errors persistwithoutpressureclamping and canfreeze
learning. k2floor50 for1023samples (~51.15s),k3zero52samples;coefficients cycle
withdirection andhistoryweightsdrift upward,not physicalparameter convergence.
Ks inferred20; staticlambda inference is invalid oncircle (nonzerodesiredrates),
do not report approximate40.8 as verifiedinput. AngleMAE.195rad~11.2deg with
phierror~.03; direction dominates lengtherror. Effective lag estimate fromMAE
alone not uniquelyidentified. No source/model edits. Keep30s baseline; inspect
signedphaseerror/request-vsallocatedforcebytheta before speed increase or
adaptationrate tuning. Need matchedPAC-vsFBL referencecomparison later.

## 2026-10-13 - Repeat10circle1rad30s: persistent phase error and allocation mismatch
Found PAC_circle/AdaptiveController/1rad_30s.mat matching6720samples0-335.95.
10complete circles20-320,period30s. Lap2thetaMAE.2890 tolap10.2752rad;
phiMAE.0487->.0579,lengthRMSE2.784/2.895->2.756/2.601mm; littleconvergence.
Lap3signedtheta.28584 equalsMAE:dominantlylag;equivalent1.36s at.20944rad/s,
not uniquelypuretransportdelay. Aboveactivationfractionlap2.355->lap10.278,
residualpeaks up to[3.647 5.561]N withalmostnolimitedflags. Lap3residualgate
passes~32.9% samples. k2floor50 on1555samples(~77.75s),k3zero47.
Adaptationupdates2261whilehistoryweights/damping/biascontinuechanging;notmodel
convergence. Allocationactive-map mismatchprioritybeforegainincrease: existing
fixedmeanlinearallocation allowscommandsbelowactivation thenforceevaluationclips,
nominalresidualcanblocklearning. Needpressure-vsangle/signedphaseandfeasibility
review;activationmapassumednotchambermeasurement. Endthetaerrorwhilephifalls
excludeasunobservable. No code changes.

## 2026-10-13 - Activation-consistent fixed-mean PAC allocator
User authorized pacAllocatePressure change. Exact activation-map inverse now
uses requested+forceBias to derive effective pressure differences. Preserves1.6
bar mean; when effective muscle pressure iszero, its inactive commandedpressure
can decrease below.8 to satisfyfixedmean. Candidate accepted only withinexisting
.2-3bar rails ands/a +/-2.5trimcorrectionclamps. Infeasiblecandidate retains bounded
legacyfallback,nowlimited=true fornonzeroallocationresidual; no hiddenmeanincrease
or relaxedconstraints. Supervisor2bar/slew remains andmayproduce transientresidual.
Nominal activationthreshold.8 remains assumption; doesnot identifyphysicalpressure.
validatePacPressure includes feasibleinactiveexactforceand infeasibleflag checks;
validatePacAdaptive passes. No SLX edited byagent. Repeat1rad30s baseline,
lambda40/Ks20 andrecordallocationresidual/limits/updates beforehighergainorspeed.

## 2026-10-13 - Fixed-mean activation inverse circle did not improve tracking
Raw1rad_30s_tr1.mat newallocator:lap10thetaMAE.3094 vsprior.2752,phiMAE.0888
vs.0579,RMSE3.141/2.821 vs2.756/2.601mm. Newlimited semantics includes
forceinfeasibility:cannotcompare3184vs25 asphysicalclampingcounts. Lap10limited
56.2%,activationabove20.2%. Offlinecandidateconstraints50-320 valid5384samples:
minimumviolations2244,maximum1486,scorrection496,acorrection0(overlapcounts).
Therefore fixedmean1.6 cannot exactlydeliver manyrequestswithinallbounds;
activation consistency improvement didnotresolvephysicaltracking. Needbounded
variablecommonpressure allocation/closestfeasibleforce solver to avoid legacy
fallback mismatch. Do notrelax rails or assumehighercommonalwayshelps. Proposed
next allocator changes separateauthorization;keepgains/rates forcomparison.
No controller edits in this review.

## 2026-10-16 - Documentation clarification and commit preparation
Pressure nullspace doesnot require3mechanicalcoordinates:retain2lengthcoordinates
and3pressureinputs;optional pneumaticpressurestates modelvalvedynamics. Fixedmean
allocation contributesforce mismatch;variablemeanproposal remainsunimplemented.
Wrappedtheta error blocksnippet supplied inchat atan2(sin(des-meas),cos(des-meas));
not added asstandalonefile;theta nearstraight poorlyobservable.
User requests allpendingupdates committed. Includesworkspaceanalyzers,PACplotting,
activationallocator,models,setup scripts,recordings/reportsandlogs. Offline
validatePacPressure/validatePacAdaptive passed; no hardwaremodel executed.
Allocator circle limitation retained/documented,not presented asperformancefix.

## 2026-10-16 - Switchable variable-common allocation and external PAC tuning
User authorized commonpressure adjustment andexternal tuning ports. New
pacAdaptiveTunableBlock keepsfirst11inputs/26outputs andaddsBoolean commonenable,
10fixedsize tuningvectors plusmeancommand/meanshift outputs. pacTunableDefaults
andTUNABLE_BLOCK.md providewiring anddefaults. Oldwrapper remainsunchanged.
pacAllocateVariableCommon projects2Dnominalforce ontoallactivepressure/correction
polygon usingboundededge/vertexenumeration; exactfeasibleforces,closestEuclidean
infeasibleforces,meannearestnominal. Keepsrails/slew/correctioncaps; no gradient
commonstate, noidentificationofphysicalpreload. OFF usescurrentfixedmeanallocator;
ON avoidsnominaldeadzone, canchangephysicalstiffness. Supervisortrim/anchor/fault
behaviorunchanged; switchedtargets passsamecommandslew. Update/historyaccept
optionalcfgdefaultingoldsettings. Externalgamma/scales/ratecaps/bounds/initial,
error/residualgates/historylaw/pressuremap/supervisortiming; mechanicalmodelfixed.
Structural/calibration settingsbefore-runonly;initialreadonfreshpersistentstate.
OfflinevalidatePacPressure/Adaptivepassed;validatePacTunable testsprojectionforce
bounds/feasiblevertices,OFFequivalence,custominitialvalue,actualmeanchangeunder
nonzerorequestand commandslew. No SLX/hardware runs. Oldfixedmeanexperimental
failurepreservedindocs. Provisionalcomparisonlambda40/Ks20,30scircles,sameinit.

## 2026-10-16 - Variable-common circle1rad20s verified
Rawexperimentations/adaptCommonP/1rad_20s_tr1.mat0-240 reviewed,10circles20-220,
period20s. Actual loggedlambda[28 28],Ks[5 5] constant,not prior40/20.
EnableAdaptation/commonenable loggedtrue at0/.02only(constants); no later sampled
coverage,butobservedadaptation/commonmotionconfirmsbehavior. Meancommandmatches
commonPressure exactly;range1.5971-1.9109bar;pressuresmin[.8 .802 .8],
max[2.5246 2.8033 2.3981]bar. Nearlyall validcircle samplesforce residualroundoff,
exceptlap7(140-160)peak4.344/.3663N,2.2%limited;restlaplimiting0.
Adaptiveactive97-99.5% eachlap. Lap1phiMAE.1117 theta.4230(lengthRMS3.048/4.249mm),
lap2theta.4359,lap6.3459,lap8.2813,lap10.2884;lap10phiMAE.0642,
lengthRMS2.536/2.717mm. Allthetasignedmeans=MAE:directionlags. Lap10equivalent
lag.918s at.31416rad/s,notpuretransportdelayproof. Lap7disturbancerecovers.
Allocator achievesnominalforceconsistency/avoidsthreshold betterthanfixedmean,
notactualforcevalidation;changedgains/speedprecludecontrolledtrackingcomparison.
Next matchedON/OFF withsame28/5,20s orreviewlaterbettergainrecording.
No controller/model edits.

## 2026-10-16 - 100g circle comparison
Offline adaptCommonP/1rad_20s_tr1_100g.mat userlabels100g. Same logged20s10circle
phi1reference,lambda28/Ks5,adapt/commonenabled. Payloadattachmenttiming/typeand
whethermodelpayloadchangednotlogged/verified; do notinfermassidentification.
No faults;lapvalid99-99.8%,activesteps98-99.5%;almostnolimits. Common1.6-1.851bar,
pressure.8-2.793bar. Nominalresidualmostlyroundoff;briefpeaks.1485Nlap4,.101lap6.
Lap2thetaMAE.4375 ->lap10.3112,phiMAE.0736->.0636,lengthRMS3.527/4.070 ->
2.842/2.872mm. Unloadedlap10theta.2884,phi.0642,length2.536/2.717:
loadedtheta~8%worse,phiaboutsame,length~12%/6%worse,notrepeatabilitystatistics.
Finalsprings579.45/136.85/578.8,damping43.579/44.363,history3.774/3.507/3.587,
bias[-.67175 .77705]. Parametercompensationdoesnotidentify100g; commonpressure
alone cannotcompensategravitybending becauseidealnullspace. Evidenceofmaintained
boundedtrackinginonepayloadrun,notproofadaptationbenefitwithoutmatchedOFFtest.
No code/model edits.

## 2026-10-16 - Estimate-output dropout and adaptation-rate audit
Read-only MATLAB inspection of adaptCommonP/1rad_20s_tr1.mat: 4801 samples,
0-240 s at 0.05 s. All 31 all-zero damping-matrix samples coincide with
modelValid=false and status=4 (invalid-data hold). Stored parameterEstimate
damping remains >=40; these are force-evaluation output zeros, not estimate
resets. pacAdaptiveForce initializes its diagnostic matrices/forces to zero
and returns early on invalid inputs. Underlying invalidity cause not isolated.
User screenshot extends to ~330 s and cannot be assumed to be this recording.
Observed peak absolute parameter rates in standard order:
[200 200 200 .0995549 .0849612 .060864 .0680641 .0639981 .2 .2].
Spring and bias rates reach current source maxRate values; damping/history
rates remain below current caps. Gamma/maxRate tuning inputs are not logged,
so exact run-specific values and clipping fractions are not established.
To judge gamma, inspect gates, bounds, unclipped gradients versus rate caps,
then use matched repeated gamma comparisons; residual tracking error alone
does not establish insufficient gamma. No code/model/hardware execution.

## 2026-10-19 - First RL open-loop30s trial validated
forward_test_RL_30s.mat read offline. Requested/NI-pathP1identical,exactwaveform
1.5*(1-cos(.2*t))bar(range0-3),P2=P3=0. Period31.416simulationseconds,
30scontains.955cycle. pressureDesired_bar nowdistinctinverseprediction;
pressureUsed_bar max2.2 ismodelmax(p-.8,0),NOTphysicalcommandclip.
601NDIsamples20Hz,600valid,invalid8.30s;pairedframes601,norepeats,basevalidall.
10overruns;hostelapsed31.484536vs30simulationseconds,hostdtmedian.05629,
min.021985max.347293. Needhost/simulationmappingandrecordedcommandtimes;
NIwritewalltimestampsnotpresent. Do notassumeuniformphysical50msormarkall
signalswithsamewallclock. Onlyphysicallograte20Hzdespite100Hzrepeatedlogs.
moduleTipXYZValid600;recomputedhelpermatchesloggedXYZexactlyonvalidsamples,
invalidXYZNaN. moduleendpointbounds[-.12238,-.001387,.09710]to
[-.002866,.014492,.17554]m;phiOrientationmax1.88401rad. Rawpositions/errors,
Rbase,gravity,flags,host/poll/frame,commands/volts present. Rawquaternionsabsent.
ConstantP2/P3/killlogonly2points0..18.86;notmissingacquisitionpersebutexportmust
handleconstantlifetimesexplicitlyorpersistmetadata. Othermodel/sensorname
collisionsandunnamedsignalsremain;exportbyprovenance. Datasetcontainsoutonly,
noembeddedtrialmanifest. Goodplumbingtrial,notcompleteRLtrainingbatch. No model
orsourcechanges, nohardwareexecution.

## 2026-10-19 - Open-loop RL pressure source blocks
User requests pressureexcitationblocks forreal-arm complete-dynamicslearning.
Created RL_data_collection/ standalone steppedramp,singlemusclesine,andcombined
common/s/a/rotatingpressureblocks. One sampledt input,internaleditableparameters,
outputsrequest3barvector/trial/phase/channel/cycle/done/endtime; nofeedbackor
activation subtraction. SineP1default0..3bar,.2rad/s,3cyclesx2trials,initial5/rest10,
total213.496s. StepP1levels1/2/3x2,rise5/hold10/fall5/rest10,total185s.
Combined4modes,.2rad/s2cycles,bias1.5,diff.5/commonpeak2,total336.327s.
KeepdownstreamKill/saturation,logfinalNIcommands andsensorhosttimestamps.
Datasetactionrawcommands,notpressureUsed ormodelpredictions. Existingtrialusable;
nochange tosimulatedforwarddeadzone withoutseparaterequest. READMEcoverswiring,
phaseids,metadata,stopdurationandconstantholdsemantics. Offlinebounds/slew/end/
invalidtime/sineequivalence checks passed. No SLX/hardware run.

## 2026-10-19 - Sequential pressure defaults and circle/helix references
User requests sequentialP1/P2/P3 andnew circlephi.2/.8/1/1.2/1.5,20speriod,
wait3s,fivecyclesandreturnzero;helix0->1.5loading/unloading3-4trials.
PressureSteps/Sine defaultsnowchannelSelection0(allseparately),durations545/
630.487s. AddedrlCircleSweepReferenceBlock620s(defaultphiRise/fall5s,zero10),
rlHelixReferenceBlock525s(4trials,3turnsload+3unload,zero10,alternateCWCCW).
Outputsdesiredtheta/phi+labels;NOTpressurecommands/actualpaths. Openlooprealization
requiresfeedforwardinversepath,notcontroller,notjustcomparisonportrewiring.
No per-laptheta pauses;phismooth,theta speedswitchesatphaseboundaries.
READMEwiringandscheduleexplicit. Offlinepressure/referencecontinuitybounds/speed/
endtests passed; no SLX edits orhardware execution.

## 2026-10-19 - RL CSV/report export and first full sine package
User requests reports,experimentinformation andCSVexportsfornewrecordingsin
RL_data_collection. AddedexportRlTrainingData.m:individualMATorfolderbatch,
per-signalnativeCSVwithtime/provenance,NDIsensoralignedobservations,NIcommands,
separatemodelpredictions,trial/phaseintervalsandtrialsummary,dictionary/inventory,
JSONmanifest,qualityreportandhostgapCSV. Ambiguousnameserrorunlessphysicalpath
resolves;constanttailsonlyexplicitlynoted;nonconstantnopastextrapolation.
OriginalMATunchanged. Optionalmetadata suppliedbyuser;missingmetadata flagged.
Physicaltransitiongapdefault.10s,invalid/gaptransitionsflagged;noautomatictrain/test
splitsorfakechamberpressure. DATASET_HANDOFF.md describessemanticsandlimits.
ActualSineWave.mat636simseconds:12721sensorrows12710valid,allpairedframesfresh,
17overruns12hostgaps>.10s,max4.499s;hostelapsed648.3145s. Sixlabelledtrials
P1/P1/P2/P2/P3/P3,threecycleseach,rawcommand0..3bar. Measuredphi maximaP1~1.95,
P2~1.96,P3~2.51rad:notnewsaferanges. Reviewedpackageexports/SineWave_reviewed,
ZIPexports/SineWave_RL_export.zip. Earlierinterimexportsalsopreserved.
CSV/JSONroundtrip,masks/transitiongaps,dimensionsandinitial30strialsparseconstant
checks passed. No SLX/hardware executed. NewdatareportsdonotcertifyRLcoverage.

## 2026-10-19 - Export mappings verified for all five experiment types
SavedcircleSweep/helix/combinePressure/StepWave/SineWave.mat inspectedoffline.
PhysicalthetaOrientation/phiOrientation fromNDIMATLABSystem andfinalNIpressure
signals commonacrossall. Exporterenhanced: requestedtheta_from_block/phi_from_block
becometheta_reference_rad/phi_reference_rad separatefrommodeldesangles. Whenno
rawpressure-sourcevector andanglereferencespresent,pressureDesired_bar supplies
pre-valvefeedforwardrequest,whileNIcommandsalwaysauthoritativeaction. Model
predictiontable retainsinversecalculationprovenance. All5summaryexports produced;
verifiedmandatorypressure/measuredcolumns andcircle/helixreferencecolumns.
Sensorrows/valid:Sine12721/12710,Step11001/10992,circle12501/12474,
combined6841/6832,helix10601/10575. No MAT orSLXmodified,nohardware executed.

## 2026-10-26 - Signed stiffness bound assertion
Fixed strict k1 lower-bound check (> -100) to >= -100. Accepted minimum bounds are [-100;-500;-100]; damping and hysteresis checks remain unchanged. Extended validatePacSignedK2 to verify all three signed stiffness boundaries, outward projection and inward recovery. MATLAB offline validation passed. No tuning-vector changes, model edits or hardware execution. This validates implementation only, not closed-loop stability with signed stiffness.

## 2026-10-26 - Framework and experiment checkpoint
User requested a repository checkpoint. Includes the signed stiffness helper/test, framework handoff plan, current user-saved PAC model, v10 recording/reports/audit scripts and newly supplied 1rad_4s_gainTesting_v11.mat. v11 is preserved as unreviewed data; no performance or configuration conclusions have been established for it. Signed-bound offline validation passed earlier in this session; no hardware model was run. Generated RL exports, duplicate presentation ZIP, model backups, temporary code-generation/PDF files and unrelated paper drafts remain local. No remote push requested.
## 2026-10-28 - Offline continuous-tuning supervisor prototype

User requested unattended between-batch tuning. Added host-only
pacTuneOnlineConfig/pacTuneOnlineStep and ONLINE_WORKFLOW.md. Synthetic
validatePacTuneOnline passed in MATLAB R2025a: exact candidate acknowledgement,
hold/settle/run/score progression, identical-gain repeat-to-target, stale-score
rejection, invalid coverage, limiting/fault stop, bounded local BO proposal,
global budget and apply timeout. No model was loaded or run, gains applied,
or SLX edited. Sandbox MATLAB startup failed with filesystem inconsistency;
the same offline test succeeded outside the sandbox.

This does not validate hardware integration, trajectory synchronization,
independent watchdog behavior, transport, score collection or hardware safety.
The prototype has no such adapters yet. Preserved-state comparisons
remain history-dependent, and repeat passes are operational checks rather than
independent reset confirmations. Scale-5 tuning requires a new campaign with
explicit current gamma and bounds; no historical scale-1 campaign was reopened.

## 2026-10-28 - Continuous tuning blocks, worker and new model copy

User authorized implementing the Simulink gain/trajectory integration and
provided baseline gammaK/D/H=.0515/.0005999/1.49999. Created
PAC_autotuning_online.slx from saved PAC_autotuning_manual_v2.slx without editing
the source. Saved source inspected from ZIP and model loaded/saved only; it was
not updated, compiled or run. No hardware model was executed.

Added fixed-size pacOnlineRuntimeStep, standalone generated block library,
separate-process pacOnlineWorker/WorkerStep, nonblocking localhost bridge,
configuration/snapshot creation, settings packer, zero-command stop confirmation
and reproducible new-copy installer. Runtime handles actual reference phases,
gain acknowledgement, independent metric accumulation, terminal ramps, heartbeat,
settings/gamma mismatch, faults, coverage and tracking limits. See
autotuning/ONLINE_SIMULINK_WORKFLOW.md for detailed contract and launch.

Offline verification:

- validatePacOnlineRuntime passed synthetic cycles, repeated-target finish,
  bounded local BO, duplicates, heartbeat, bounds/configuration/limit/clock stop.
  Online score matched independent pacTuneFixedVelocity offline reconstruction
  to <1e-10 in the synthetic fixture. Toy convergence fixture uses relaxed
  targets; real session retains 2mm/.01m/s/.025m/s^2 targets.
- Standalone library compiled/ran in new synthetic Simulink harnesses, both
  network disabled and localhost-enabled, without device blocks.
- Two real separate MATLAB processes exchanged synthetic UDP data and completed
  two batches with target-stop handshake; worker persisted aligned 121-column
  packet audit and checkpoints. Latest fixture under tmp/pac_online_transport_20261028_v2.
- A separate full-model test clone had BOTH complete NDI and valve subsystems
  removed and replaced by mock sources/terminators before compilation/run.
  This exposed and fixed measured vector width (four, not six), numerical
  enable-port types, library links and settings-feedback delay. The device-free
  full wiring compiled and ran. Actual arm model remains uncompiled/unrun.
- Test warning identified the saved dashboard pressure control's obsolete
  Kill:Value binding. The new copy explicitly binds to its actual pressure-enable
  constant, defaults Off and requires operator On before Run. Original untouched.

Scope limitations: no hardware timing/safety validation, deployable real-time
target, contextual GP, speed schedule or automatic rollback. Normal simulation
extrinsic transport only. Runtime watchdog is independent of worker, not of a
stalled Simulink/PC. UDP packet audit is not a substitute for normal logsout.
Adaptive state is preserved, so target repeats are operational rather than
independent reset confirmations. New session uses scale-5 settings, not old
campaign observations. No commit/push performed.

## 2026-10-28 - Dedicated baseline-only smoke configuration

Prepared PAC_autotuning_online_smoke.slx and separate
autotuning/online_sessions/smoke_phi1_T5_scale5_r01 from the implemented session.
Three circles at 5 s, phi1, final two scored. GammaK/D/H fixed to user baseline
.0515/.0005999/1.49999 with equal lower/upper bounds, budget one, no gain search.
Model created/saved only, not compiled/run. The same integration functions have
already passed device-free clone tests. Worker pure-state check verifies one
completed observation leads to budget stop without an optimizer proposal or
gain change. Operator-run monitored smoke and saved-log review remain pending.
Expected normal runtime reason23; missing tracking target alone is not a smoke
failure or proof of instability. Preserve raw logs including aborted trials.

## 2026-10-28 - Desktop Real-Time pacing incompatibility and fresh smoke retry

User screenshot: NDI preflight passed, then simulation refused to start because
Simulink pacing cannot be combined with Desktop Real-Time. This was an installer
mistake: device-free tests intentionally disabled pacing, so did not exercise
that hardware-specific compatibility check. Installer changed to pacing off;
generated initial/r01 model files corrected without running/compiling hardware.

Separate worker screenshot showed zero packets/zero completed batches, followed
by expected overwrite rejection on restart. Old log preserved. Prepared new
smoke r02 model/session with the same fixed gains/three circles, new identity
and worker initial connection timeout thirty minutes rather than five. Active
runtime watchdog and post-connection telemetry timeout remain unchanged.
R02 startup and physical behavior still require operator verification.

## 2026-10-28 - Windows build filename limit on r02

User screenshot identified generated C source path >260 characters under nested
session cwd/slprj. Added model-process pacOnlineConfigureBuild using
Simulink.fileGenControl cache/code folders under short local TEMP/pac_bld/PID.
All generated model-open scripts call it. Process isolation avoids shared-build
collisions; project, old generated files and logs preserved. No hardware model
run or Windows registry/policy change. User can retry same live r02 session
after calling helper in model MATLAB if worker is still listening.

## 2026-10-28 - Physical r02 smoke: false validity abort, shutdown verified

Read saved out and worker MAT only. Startup status0/1/2/3 at0/10/15/20s,
hold20.55, gamma acknowledgement20.60, approach22.60, spin-up27.60,
constant-speed run29.60. Supervisor abort at34.00 reason4, only4.4s/0.88
circles; zero completed batches, no final scoring window or optimizer trial.
All gamma values constant .0515/.0005999/1.49999; scale/rate/upper for H5/10/50.
Tracking limiting0, latched PAC status5 faults0. Isolated50ms status4 samples
recovered individually. Pressure peaks[2.13823,2.42442,2.60675]bar; final command
and threeDAQ outputszero. Terminal phase10 at39.00, commands<.01 at39.30,
simulation stop40.30 after1s dwell. Command shutdown verified, physical pressure
not sensed. Partial early RMSE[5.2206,8.9064]mm is not settled/target score.

Root cause: current raw validity false33.90, true33.95, false34.00; delayed
PAC status/modelValid remained invalid33.95 from prior frame. Old grace used
controlled (mixed current and delayed signals), bridging isolated dropouts.
Floating-point 34-33.9 also exceeds.10 slightly. Changed runtime grace to raw
current measurement validity, retained delayed flags for conservative score/
hold checks and PAC's own fault logic, added1e-9 tolerance. Hold stage no longer
calls one invalid frame movement; fresh measured movement still aborts.

validatePacOnlineValidity passes recorded-prefix replay through34.20s without
false abort, sustained-invalid stop beyond exact grace boundary, invalid hold
frame and actual-movement cases. validatePacOnlineRuntime synthetic suite
passes. Recorded replay beyond old desired-trajectory divergence cannot prove
later tracking; no such claim made. Reports MAT/CSV/README saved under
smoke_phi1_T5_scale5_r02/review_v1. Fresh r03 model/session with identical fixed
gains and no search prepared; runtime/worker source snapshot archived there.
No hardware model run by assistant. R02 is not a completed smoke pass; repeat
monitored fixed-gain smoke before optimizing.

## 2026-10-28 - R03 no-motion startup: worker handshake missing

Verified out metadata is actual R03 model. Online phase0 ->9 at50.05 reason9
->10 at55.05; end56.30. Actual online phi/theta0, candidate0. Status3 from20s,
all measurement/base/model flags valid, post20 velocity peaks[.000187003,
.000134689]m/s below.0005. Reconstructed non-heartbeat readiness from20.55s
for590 phase0 samples. Therefore startup's heartbeat gate was missing, not
insufficient learning or stationarity. Final pressure commandszero. Existing
unused old reference logs can show circles despite actual online targetzero.

Read-only netstat: both localhost32151/32152 on MATLAB backend28480, R03 cache
under60752. This supports stale model socket retained in instance now running
worker, blocking the other model process. Worker lock present, no final log.
Previously bridge initialization catch silently returnedzeros, hiding bind
error. Added explicit warning/transport flag and runtime reason13; worker
startup clears old bridge only with local models stopped. No user's sockets
or processes terminated by assistant.

New pressure-free connection probe validates matching session/baseline and
port ownership before model Run. Worker probes do not activate telemetry timeout
or enter observations. Alternate-port32251/32252 synthetic regression deliberately
seeded leftover model socket in worker, verified startup cleanup, probes11s
apart, subsequent two-batch completion and persisted audit. Existing R02
validity regression still passes. Reports under R03/startup_review_v1; old
raw recording preserved. No hardware runs/gain changes performed.

## 2026-10-28 - R033 completed motion, insufficient coverage, rejection handshake fix

Read saved R033 out and1733worker packets only. Three5s circles29.6-44.6s,
gamma .0515/.0005999/1.49999 constant, H scale/rate/upper5/10/50. No tracking
limiting or latched PAC fault. Last2circle exploratory scores independently
recomputed to<1e-8: lengths[4.681777,6.780124]mm, velocity error.01210224m/s,
roughness.03081964m/s^2. Coverage.884422 excludes GP training acceptance.

Score200samples:192rawvalid,184controlled after delayed flags,176/199pairs.
Eight tip-validity dropouts withbase valid/overrun0. Raw-only pairs183/199 also
below95%; exact missing/stale/tracking cause unresolved, policy/sensor code
unchanged. Per-lap length errors[5.3172,8.3457],[4.4699,7.1445],[4.8845,6.3950]
mm; not settled-performance proof. One motion batch, zero accepted observations.

Worker rejected usingoldrequest1 while runtime requested2, causing30s extra
hold. Reason11 at81.6, phase10 at86.6, zero commands86.9, stop87.9, finalDAQzero.
Pressure peaks[2.6315,2.6191,2.6050]bar. Fixed rejection reply ID, explicit
reasons, terminal immutability and stop withoutstationary wait. Invalid motion
attempts now retained/checkpointed separately from accepted observations.
validatePacOnlineRejection verifies recorded failed score ->reply2, retained
attempt/no accepted observation, stable repeats, immediate reason24 transition.
Runtime and R02 validity regression suites pass. No hardware/model execution
during review. Detailed artifacts under R03/review_r033_v1. Existing session
cannot be reused; no new hardware trial queued or gains changed.

