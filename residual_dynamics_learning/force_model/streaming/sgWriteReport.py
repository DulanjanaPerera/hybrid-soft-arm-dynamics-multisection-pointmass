"""Summarize completed saved experiments; no fitting or artifact selection."""
from pathlib import Path
import csv
import statistics as st
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CAMP = HERE / 'campaign_20261103_v2'
EVAL = HERE / 'evaluation_20261103_v1'
INV = HERE / 'inverse_20261103_v1'

def read(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def num(row, field):
    return float(row[field])

def mean(rows, field):
    return st.mean(num(r, field) for r in rows)

def ms(values):
    return f'{st.mean(values):.4f} +/- {st.stdev(values):.4f}'

assert (EVAL / 'VERIFICATION.md').exists()
screen = read(CAMP / 'screen.csv')
final = read(CAMP / 'finalists.csv')
one = read(EVAL / 'one_step.csv')
short = read(EVAL / 'short_paths.csv')
baseline = read(EVAL / 'short_baselines.csv')
quality = read(EVAL / 'force_quality.csv')
inv = read(INV / 'inverse_cases.csv')
coverage = read(HERE / 'data_20261103_v1' / 'coverage.csv')
history = read(CAMP / 'history_comparison.csv')[0]
ids = list(dict.fromkeys(r['id'] for r in final))
seeds = [113, 211, 307]
oldmeans = {i: mean([r for r in final if r['id'] == i], 'oldScore_mm') for i in ids}
winner = min(ids, key=oldmeans.get)
representative = next(r for r in final if r['id'] == winner and int(r['seed']) == 113)
config = next(r for r in screen if r['id'] == winner)

def macro(i, seed, table, field, horizon=None):
    rows = [r for r in table if r['id'] == i and int(r['seed']) == seed]
    if table is one:
        rows = [r for r in rows if r['predictor'] == 'streaming_force']
    if horizon is not None:
        rows = [r for r in rows if abs(num(r, 'horizon_s') - horizon) < 1e-8]
    return mean(rows, field)

one_stats = {i: [macro(i, s, one, 'vectorRMSE_mm') for s in seeds] for i in ids}
path_stats = {i: {h: [macro(i, s, short, 'pathVectorRMSE_mm', h) for s in seeds]
                  for h in [.1, .25, .5, 1]} for i in ids}
nom_one = mean([r for r in one if r['predictor'] == 'nominal'], 'vectorRMSE_mm')
cv_one = mean([r for r in one if r['predictor'] == 'constant_velocity'], 'vectorRMSE_mm')
eligible = [r for r in inv if r['knownTargetFeasible'].lower() in ('1', 'true')]
latencies = sorted(num(r, 'solve_ms') for r in eligible)
p95 = latencies[round(.95 * (len(latencies)-1))]
if (INV/'summary.json').exists():
    p95=json.loads((INV/'summary.json').read_text())['latencyP95_ms']
force_error = max(num(r, 'forceError_N') for r in eligible)
bound = max(num(r, 'boundViolation_bar') for r in eligible)
singular = min(num(r, 'minimumSingular_NperBar') for r in eligible)
force_max = max(num(r, 'maxForce_N') for r in quality)
force_bound = max(num(r, 'boundFraction') for r in quality)
rate_p95 = max(num(r, 'rateP95_Nps') for r in quality)
failures = sum(int(r['failures']) for r in one if r['predictor'] == 'streaming_force')
failures += sum(int(r['failures']) for r in short)
anchors = sum(int(r['anchors']) for r in quality if r['id'] == winner and int(r['seed']) == 113)
chunks = sum(int(r['chunks']) for r in coverage if r['partition'] == 'train')
selection_chunks = sum(int(r['chunks']) for r in coverage if r['partition'] == 'selection')

lines = [
    '# Streaming generalized-force residual: completed experiment sequence', '',
    '2026-11-03. Eleven completed offline fits: six architecture/horizon screen fits,',
    'one burn-in comparison and four additional finalist-seed fits. October19/27',
    'only for training and selection. All finalist weights and the representative',
    'research checkpoint were frozen before this November2 replay. The previously',
    'inspected November2 data are evaluation only, not blind confirmation.', '',
    f'The selected configuration is **{winner}: {config["widths"]} units, ',
    f'{num(config,"horizon_s"):g}s maximum physical training rollout**. Selection uses',
    'the lowest mean older-holdout score across three seeds; seed113 was declared',
    'in advance as the representative, rather than picking the luckiest seed.',
    'Research checkpoint: `release_20261103_v1/force_model.mat` in this folder.',
    f'SHA256: `{representative["sha256"]}`. Previous parent release is unchanged.', '',
    '**Main limitation:** the force GRU substantially improves the nominal model,',
    'but constant velocity remains better over all tested horizons through1s.',
    'The numerical inverse works on self-generated targets; this does not yet',
    'justify physical pressure correction or controller deployment.', '',
    '## Contract and checks', '',
    'The GRU updates each layer once from current measured q, causal estimated v,',
    'current pressure command and previous physical dt. Its outputs are two',
    'generalized forces in N on the mechanical RHS; lengths follow from the same',
    'fixed nonlinear nominal dynamics and SDIRK integrator. Chronological .5s',
    'training chunks carry measured-input hidden states, detached between updates.',
    'Auxiliary .5/1s free prediction starts with a separate pre-current state copy,',
    'uses integrated q/v and recorded future commands, and never commits its state.',
    'It is a known-command prediction test, not unknown-command forecasting.', '',
    f'Original whole-trial splits and TRAIN-only normalizer preserved exactly. {chunks}',
    f'train chunks, {selection_chunks} selection chunks in 241 valid mature segments;',
    '128 common selection anchors across eight eligible recordings. Day2 StepWave',
    'has two train chunks and no holdout coverage; day2 combined has neither.',
    'All anchors obey actual poll-start time caps and original validity/gap guards.', '',
    'Numerical/state checks pass for 1x32,1x64,2x32: stable gates, finite gradients,',
    'gradient versus central differences, chronological state replay, masked carry,',
    'future-target independence and auxiliary-state isolation. Full anchor time/',
    'segment/partition audit passes. Saved one-step metrics independently recomputed',
    'within1e-8mm; old release, normalizer, source/data/model hashes unchanged.', '',
    '## Six-model screen', '',
    'Old-only score = .7 macro one-step vector RMSE + .3 macro .5s free-path vector',
    'RMSE. Vector means sqrt(e2^2+e3^2) before time RMS, then equal recording means.',
    'This is NOT a single-coordinate or closed-loop tracking error.', '',
    '| ID | GRU units by layer | Horizon (s) | Parameters | Old score (mm) |',
    '|---|---|---:|---:|---:|',
]
for r in screen:
    lines.append(f'| {r["id"]} | {r["widths"]} | {r["horizon_s"]} | {r["parameters"]} | {num(r,"oldScore_mm"):.4f} |')
lines += ['',
    'All fits used seed113, batch32, Adam .0003, gradient clip1 and maximum1200',
    'updates with the same early-stop policy, checkpoint cadence, data, normalizer,',
    'force bound and loss weights. Length terms separately average one-step and',
    'free-path errors (.5/.5), plus .001 force-rate and .0001 magnitude penalties.',
    'Two layers did not win this screen. This is a practical capacity comparison,',
    'not parameter-matched proof that depth is inferior. Learning curves were still',
    'improving near the budget cutoff; rankings are budget-specific.', '',
    f'Leader burn-in .5s: {num(history,"burnHalfScore_mm"):.5f}mm; 1s: ',
    f'{num(history,"burnOneScore_mm"):.5f}mm. The worse 1s variant was not retained.', '',
    '| Finalist | Seed113 old score | Seed211 | Seed307 | Mean +/- SD (mm) |',
    '|---|---:|---:|---:|---:|',
]
for i in ids:
    vals = [num(next(r for r in final if r['id']==i and int(r['seed'])==s), 'oldScore_mm') for s in seeds]
    lines.append(f'| {i} | ' + ' | '.join(f'{v:.4f}' for v in vals) + f' | {ms(vals)} |')
lines += ['', '## November2 frozen evaluation', '',
    f'{anchors:,} common one-step anchors per model/seed; 60 matched free windows',
    'per model/seed, each scored within .1/.25/.5/1s physical caps. Primary replay',
    'refreshes q/v each sample; free windows have no measured-state resets.',
    f'Finalist numerical guard failures: {failures}. No failed paths silently removed.', '',
    '| Predictor | One-step vector RMSE (mm) | .1s path | .25s path | .5s path | 1s path |',
    '|---|---:|---:|---:|---:|---:|',
]
for i in ids:
    lines.append(f'| {i}, three seeds mean +/- SD | {ms(one_stats[i])} | ' +
                 ' | '.join(ms(path_stats[i][h]) for h in [.1,.25,.5,1]) + ' |')
for label,value in [('nominal',nom_one),('constant_velocity',cv_one)]:
    vals = [mean([r for r in baseline if r['predictor']==label and abs(num(r,'horizon_s')-h)<1e-8], 'pathVectorRMSE_mm') for h in [.1,.25,.5,1]]
    lines.append(f'| {label} | {value:.4f} | ' + ' | '.join(f'{v:.4f}' for v in vals) + ' |')
lines += ['',
    'The released seed113 representative has November2 one-step macro RMSE',
    f'{macro(winner,113,one,"vectorRMSE_mm"):.4f}mm and .5/1s path RMSE',
    f'{macro(winner,113,short,"pathVectorRMSE_mm",.5):.4f}/',
    f'{macro(winner,113,short,"pathVectorRMSE_mm",1):.4f}mm. The table above reports',
    'the full seed distribution, not only that representative.', '',
    'Mean/spread are across three fixed seeds of equal-run macro scores, not',
    'confidence intervals or independent hardware repeats. Per-coordinate RMSE,',
    'p95 absolute error, waveform phase and measured-speed groups are in CSVs.',
    'Speed groups are norm(v)<=.005, .005-.015 and >.015m/s; these are measured',
    'motion groups, not new cycle-speed experiments.', '',
    f'Max residual magnitude across finalist replays {force_max:.4f}N; maximum',
    f'per-run occupancy within1% of the80N head bound {100*force_bound:.4f}%;',
    f'largest per-run force-rate p95 {rate_p95:.4f}N/s. Force quality and all spikes',
    'are retained without post-hoc filtering. Rates use adjacent samples within',
    'valid segments; reset boundaries remain explicit in saved segment labels.',
    'These are effective compensation',
    'forces, not measurements or uniquely identified physical parameters.', '',
    'The parent release used a different warm-up/free-prediction contract and',
    'anchor set. Its previous 2.151mm one-step result cannot isolate a width/depth',
    'effect against this streaming experiment. No post-November2 tuning occurred.', '',
    '## Offline pressure inversion', '',
    'Best mean OLD-holdout configuration, fixed seed113. A known logged pressure',
    'generates its own target through B*p+r(q,v,p,h_previous). The solve holds',
    'common pressure at the logged mean, with0..3bar bounds and0.3bar step cap',
    '(6bar/s at intended.05s). The0..3bar range is the offline RL command domain,',
    'not a changed PAC setting. A deliberately displaced starting guess avoids',
    'trivial solves on holds. Trial optimizer calls reuse immutable pre-current',
    'hidden states; only the selected pressure commits one new state.', '',
    f'{len(eligible)}/{len(inv)} logged targets satisfy those prescribed bounds;',
    'excluded targets are explicitly retained as infeasible, not counted successes.',
    f'Max force-equation residual {force_error:.5g}N; max bound violation {bound:.5g}bar.',
    f'Minimum observed differential pressure Jacobian singular value {singular:.5g}N/bar.',
    f'Solve+selected-GRU latency median {st.median(latencies):.3f}ms; p95',
    f'{p95:.3f}ms; {sum(v>50 for v in latencies)}/{len(latencies)} exceed50ms.',
    f'Maximum latency {max(latencies):.3f}ms includes the first optimizer call;',
    'optimizer cold-start cost was retained, not excluded to improve the timing.',
    'Exact MATLAB summary and per-case timings/exit flags/continuity gains are',
    'retained. Latency excludes acquisition, nominal-control computation and',
    'real-time scheduling; this is NOT qualification for the50ms hardware cycle.', '',
    'Self-consistent inversion cannot establish true pressure correction or',
    'tracking improvement. Commands are pressure proxies; q is inferred from',
    'orientation. Missing measured force/chamber pressure and current-command',
    'dependence/actuator lag remain material limits. No PAC/PDFB/SLX integration',
    'or arm execution was performed.', '',
    '## Retained artifacts', '',
    '- campaign_20261103_v2: predeclared protocol, all checkpoints, source hashes,',
    '  split coverage, screen/history/seed CSVs and all six frozen finalist weights.',
    '- release_20261103_v1: selected research checkpoint, training source and hash.',
    '- evaluation_20261103_v1: all predictions/forces/truth/commands/time in MAT,',
    '  per-run/phase/speed and baseline CSVs, source hashes, independent verification.',
    '- inverse_20261103_v1: targets/pressures/states, cases CSV, timing and protocol.',
    '- Three FIG/vectorPDF/300dpiPNG sets: streaming_comparison,',
    '  circle_streaming_example and circle_short_free_example. Example windows are',
    '  first chronologically selected, with prediction mode explicit in captions.',
    '- Failed campaign_v1 retained for reporting-error provenance; no model released.', '',
    'New November3 recordings were not pooled. File migration, hardware testing,',
    'code/model deployment, Git commit and push remain unperformed.', '',
]
report = '\n'.join(lines)
(HERE / 'REPORT.md').write_text(report, encoding='utf-8')
(HERE / 'release_20261103_v1' / 'REPORT.md').write_text(report, encoding='utf-8')

entry = f'''## 2026-11-03 - Streaming force-GRU sequence completed

Explicit sequence executed: six architecture/horizon fits, one burn-in comparison,
four additional seed fits. New observed-state TBPTT primary path and independent
short free-prediction state; no future target inputs. Old October19/27 whole-trial
splits and normalizer unchanged; all finalists frozen before November2 replay.
Selected {winner},64-unit single layer,{num(config,'horizon_s'):g}s horizon,.5s burn-in;
fixed representative seed113, chosen config by three-seed mean older score.
Two-layer32+32 did not win under the common1200-update budget; not a general
depth result. November2 macro one-step {ms(one_stats[winner])}mm; .5/1s path
{ms(path_stats[winner][.5])}/{ms(path_stats[winner][1])}mm. Matched CV one-step
{cv_one:.4f}mm. CV remains better at all tested horizons through1s. All numerical
guard failures {failures}; per-phase/coordinate/
force-rate/seed results retained. New mode differs from parent release, so do
not attribute all changes to width. No post-test tuning or new November3 pooling.
Offline bounded, pressure-dependent model inversion: {len(eligible)}/{len(inv)} known
targets feasible, max equation residual {force_error:.4g}N, median/p95 latency
{st.median(latencies):.2f}/{p95:.2f}ms excluding acquisition/control scheduling.
Trial hidden states remain fixed; one selected-state commit passes. This is
self-consistency, not physical force or closed-loop tracking qualification.
Code/state/gradient/horizon/split/hash checks and independent saved-metric
recompute pass. Prior release unchanged. Reports/weights/CSV/MAT/FIG/PDF/PNG:
residual_dynamics_learning/force_model/streaming/REPORT.md. No SLX/hardware,
migration, commit or push.

'''
for name in ['RESIDUAL_DYNAMICS_RESULTS.md', 'VALIDATION_RESULTS.md']:
    path = ROOT / 'project_docs' / name
    text = path.read_text(encoding='utf-8')
    title, rest = text.split('\n',1)
    marker=entry.splitlines()[0]
    if marker in rest:
        start=rest.index(marker)
        end=rest.find('\n## ',start+len(marker))
        rest=rest[:start]+(rest[end:] if end>=0 else '')
    path.write_text(title+'\n\n'+entry+rest.lstrip('\n'), encoding='utf-8')
path=ROOT/'project_docs'/'PROJECT_STATE.md'
text=path.read_text(encoding='utf-8')
if '2026-11-03 streaming comparison IN PROGRESS' in text:
    start=text.index('2026-11-03 streaming comparison IN PROGRESS')
    end=text.index('\n\n',start)+2
else:
    start=text.index(entry.splitlines()[0])
    end=text.index('2026-11-03 user proposes two GRU layers',start)
text=text[:start]+entry+text[end:]
text=text.replace('Current Milestone: First generalized-force GRU trained on October19/27 and frozen; November2 evaluation improves nominal prediction but short-horizon accuracy remains insufficient for controller integration.',
                  'Current Milestone: Streaming force-GRU architecture/horizon/history/seed sequence complete; research checkpoint frozen and evaluated, offline pressure inversion checked; no hardware/controller qualification.')
path.write_text(text,encoding='utf-8')
path=ROOT/'project_docs'/'FIGURE_INDEX.md'
text=path.read_text(encoding='utf-8')
needle='## Other retained sources\n\n'
figure_entry='''- 2026-11-03 streaming force-GRU sequence:
  residual_dynamics_learning/force_model/streaming/evaluation_20261103_v1.
  streaming_comparison, circle_streaming_example, circle_short_free_example
  each FIG/vectorPDF/300dpiPNG. Six-model old-only screen/learning curves,
  repeated-seed results, frozen Nov2 comparison versus nominal/CV; explicitly
  separate measured-state one-step from free integrated prediction. All source,
  predictions, force spikes, pressure/time, CSV/MAT and SHA256 retained.
  Prediction/model-inverse self-consistency only, no controller/physical-force claim.

'''
if figure_entry not in text:
    text=text.replace(needle,needle+figure_entry)
path.write_text(text,encoding='utf-8')
print(report)
