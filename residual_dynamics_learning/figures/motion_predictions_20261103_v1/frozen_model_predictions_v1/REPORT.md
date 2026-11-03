# November3 frozen force-GRU motion predictions

Additional visualization/evaluation requested by the user. No fine-tuning or
normalizer change. GRU64/B seed113 remains the October19/27-trained research
checkpoint, .5s training horizon, model SHA256:
8150352dd47230a2a4094b240e1e886b641b8c22ebad0cbfe300f63e9f526c6a.

## Figures and recorded settings

circle_free_prediction, sine_free_prediction, helix_free_prediction compare
measured responses (black), force GRU64 (blue), and constant velocity (orange).
Corresponding *_one_prediction figures use refreshed measured q/v each valid
step. nominal_free_prediction retains the poorer nominal model separately.
All seven sets have editable FIG, requested-vector PDF and300dpi PNG.

Actual logged November3 conditions: circle10s per revolution, requested phi
levels[.2,.8,1,1.5,2.5]; sine period about12.575s, channels1/2/3;
helix10s/turn, requested peak phi2.5. November2 had20s circles, sine~31.4s,
helix20s/turn and peak phi1.5. These are different conditions, not paired
experiments isolating the effect of speed or date. Neither tests2s circles.

| Selected window | Trial | Simulation interval (s) | Physical duration (s) | Missing measured rows |
|---|---:|---:|---:|---:|
| circleSweep | 3 | 159.00-169.00 | 10.117 | 0 |
| SineWave channel1 | 1 | 6.60-19.20 | 12.745 | 0 |
| SineWave channel2 | 3 | 102.00-114.60 | 12.740 | 0 |
| SineWave channel3 | 5 | 197.40-210.00 | 12.752 | 0 |
| helix | 1 | 6.60-64.95 | 59.047 | 4 |

Chronological window selection by the SAME rules as the November2 plotting
code: first eligible circle at requested phi1, first eligible sine period per
channel, first loading/unloading helix remainder after observed warmup. No
window was selected by minimum model error. Helix starts1.6s after its trial
begins, to initialize observed history; both loading and unloading are included.

## Full-window free prediction

Length-vector RMSE=sqrt(mean(e_l2^2+e_l3^2)) in mm, real valid observations only.
All15 model/window paths completed without the inherited numerical guard.
Free predictions have ONE measured initialization and causal observed hidden
prefix; subsequent q/v are integrated predictions, with known recorded pressure
commands/time. No future measured q, state reset or interpolation during rollout.

| Window | GRU64 (mm) | CV (mm) | Nominal (mm) | Nov2 GRU example (mm) |
|---|---:|---:|---:|---:|
| circleSweep | 1.2636 | 1.9485 | 11.6252 | 0.9440 |
| SineWave channel1 | 1.7184 | 9.7674 | 6.7421 | 2.3355 |
| SineWave channel2 | 1.3792 | 15.7346 | 9.6997 | 1.7784 |
| SineWave channel3 | 1.6876 | 19.5617 | 7.6702 | 1.7958 |
| helix | 2.8787 | 13.1685 | 11.8083 | 0.7627 |

The GRU is closer than nominal/CV on these five long windows. Circle and
helix errors increase versus the selected Nov2 examples; sine examples differ
less. This is evidence of remaining prediction error under the newer conditions,
not proof that speed alone caused the difference. The visible low-response
region is consistent with modeling combined deadzone/friction/hysteresis/lag;
these plots do not uniquely identify a physical actuator deadzone.

## Measured-state one-step replay

New replay helper was first verified against EVERY saved Nov2 circle anchor:
exact same anchors and maximum difference<1e-12m. Current measured-input replay
uses the same observed-history initialization, valid segments, recovery warmup
and future1s eligibility guard as the previous frozen evaluation. It never
uses future motion as model input. First free steps independently match every
window initial one-step prediction<1e-10m. Per-window metrics are in scores.csv;
whole_record_one_step_scores.csv retains all eligible anchors, not only examples.
CV remains better one-step. Do not equate free trajectory accuracy with one-
step or closed-loop hardware tracking accuracy.

## Coverage and recommendation for continued training

The coverage audit uses the actual59,630 older TRAIN rows and unchanged input
normalizer; source is auditNov3ForceCoverage.m. Some circle/helix q/v now extend
beyond those observed training ranges: q2 minimum old-14.98mm versus new
circle-22.18mm/helix-20.73mm; v2 maximum old.01322m/s versus new~.0209m/s.
Per-feature outside-range fractions are about2-3.5% for affected circle/helix
coordinates/velocities. This motivates added coverage; it does not isolate the
cause of all errors. Large standardized values alone are not proof of being
outside training range. The tiny dt edge excursion for sine is negligible.

Recommend NEXT offline experiment: continued training on November3 TRAIN
trials MIXED WITH original October19/27 TRAIN data; retain November2 for
evaluation after freeze. Improvement is plausible, not guaranteed, particularly
for faster/larger motions. Do not fit November3 alone and assume old conditions
will be retained. Sequential training can forget old behavior. Replaying old
examples is a practical mitigation, not a guarantee.

Before fitting, declare whole-trial/condition holdouts within November3 and keep
original older internal holdouts. Start from frozen GRU64 weights with the
existing normalizer and lower learning rate, fixed architecture/history/loss.
Do not refit normalization under existing weights silently. The saved research
model has weights, not a full Adam-optimizer resume state; a new fine-tuning
run must declare optimizer restart. A versioned trainer/manifest is needed:
the existing old-only trainer deliberately rejects November3 and initializes
new weights; its whitelist must not be silently relaxed for historical runs.

Use three comparison arms at matched extra-update budgets:
1. Frozen model reference (already archived).
2. Additional old-TRAIN-only optimization, to control for simply training longer.
3. Old+November3 TRAIN replay, same starting weights/budget/objective.
Select on old+November3 internal holdouts only, then freeze and compare on
November2. November2 is previously inspected, not blind final confirmation;
collect a later independent repeat/day for that. Keep the imputation proposal
out of this first data-addition experiment so its effect is not confounded.

November2 could be used for training instead, but then it loses evaluation-only
status and its post-fit plots are training diagnostics. A different held-out
cohort is necessary. There is no advantage merely from choosing the earlier
date; the current Nov3 coverage and preserving Nov2 comparisons favor Nov3.

References: forgetting https://doi.org/10.1073/pnas.1611835114 ; replay
https://proceedings.neurips.cc/paper_files/paper/2019/hash/fa7cdfad1a5aaf8370ebeda47a1ff1c3-Abstract.html .

## Provenance and scope

All four helix missing truth rows remain NaN/unscored; known commands/time
continue to drive free prediction. One-step replay keeps original gap resets.
PCC XYZ is reconstructed from inferred lengths, not measured sensor XYZ.
Near-straight theta is unobservable; circular angle errors use valid directions.
All raw/prepared/export/model hashes pass; Nov2 figure/source archive unchanged.
No training, optimizer update, normalization refit, gap interpolation, SLX edit,
hardware execution, migration, Git commit or push was performed.
