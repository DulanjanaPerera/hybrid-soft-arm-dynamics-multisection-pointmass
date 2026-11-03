"""Report frozen local experiments; never fit or select using evaluation data."""
from pathlib import Path
import csv
import hashlib
import shutil
import statistics as stats

ROOT = Path(__file__).resolve().parents[3]
SRC = Path(__file__).resolve().parent
CAMP = SRC / 'campaign_20261103_v1'


def read(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    return float(row[key])


def mean(rows, key):
    return stats.mean(number(row, key) for row in rows)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def update_doc(path, title, entry):
    text = path.read_text(encoding='utf-8-sig') if path.exists() else ''
    marker = f'## {title}\n'
    if marker in text:
        before, rest = text.split(marker, 1)
        next_section = rest.find('\n## ')
        text = before + (rest[next_section + 1:] if next_section >= 0 else '')
    if path.name == 'PROJECT_STATE.md':
        lines = text.splitlines()
        lines[0] = 'Last Updated: 2026-11-03'
        lines[1] = 'Current Milestone: GRU capacity comparison complete; current GRU64 release retained pending controller integration. See latest capacity results below.'
        text = '\n'.join(lines) + '\n'
        at = text.index('# Current project state') + len('# Current project state')
        text = text[:at] + '\n\n' + marker + '\n' + entry + '\n' + text[at:]
    else:
        text = marker + '\n' + entry + '\n\n' + text
    path.write_text(text, encoding='utf-8')


def main():
    assert (CAMP / 'VERIFIED.mat').exists()
    cases = read(CAMP / 'internal_results.csv')
    short = read(CAMP / 'all_seed_evaluation' / 'short_paths.csv')
    one = read(CAMP / 'evaluation' / 'one_step.csv')
    motion = read(CAMP / 'motion_comparison' / 'scores.csv')
    names = [('gru64', '1x64'), ('gru128', '1x128'), ('gru64x2', '2x64')]
    means = {key: mean([r for r in cases if r['arm'] == key], 'compositeScore_mm') for key, _ in names}
    winner = min(means, key=means.get)
    representative = next(r for r in cases if r['arm'] == winner and r['seed'] == '113')
    assert sha(CAMP / 'candidate_force_model.mat') == representative['sha256']
    report = [
        '# GRU capacity comparison: completed offline experiment', '',
        '2026-11-03. Nine architecture/seed cases: 1x64, 1x128, 2x64; seeds113/211/307.',
        'Training gradients use October19/27 TRAIN only. November3 whole-trial holdouts participate in checkpoint and architecture selection.',
        'November2 remains excluded from fitting/selection, but is previously inspected evaluation, not blind confirmation.', '',
        '## Matched protocol', '',
        '- Stage1: original streaming protocol, up to1200 Adam updates at3e-4; checkpoint by old holdouts.',
        '- Stage2:1200 extra old-only updates at1e-4 with fresh Adam from stage1 selected weights.',
        '- Same batch32, clip1, .5s measured-state TBPTT chunks/.5s auxiliary free path, .5s burn-in, losses, physics and integration.',
        '- Unchanged old TRAIN normalizer, old/Nov3 splits and208 internal anchors. No interpolation or new training data.',
        '- Stage2 checkpoint score=.5old+.5Nov3; each cohort=.7one-step+.3half-second path, equal recording means.',
        '- Same10% old-regression guard relative to EACH case\'s stage1 checkpoint; starting checkpoint eligible.',
        '- Architecture selected by three-seed mean composite; seed113 fixed in advance as representative, not luckiest seed.',
        '- Exact existing64-unit stage1 fits and seed113 continuation reused by hash: six new stage1 and eight new continuation fits.',
        '- Same update budget is not proof all architectures converged. Parameter counts differ; no parameter-matched depth claim.', '',
        '## Internal architecture comparison', '',
        '| Architecture | Parameters | Mean score (mm) | Sample SD | Seed113 | Seed211 | Seed307 |',
        '|---|---:|---:|---:|---:|---:|---:|',
    ]
    rows_out = []
    for key, label in names:
        rows = [r for r in cases if r['arm'] == key]
        values = [number(r, 'compositeScore_mm') for r in rows]
        row = f'| {label} | {rows[0]["parameters"]} | {stats.mean(values):.5f} | {stats.stdev(values):.5f} | ' + ' | '.join(f'{v:.5f}' for v in values) + ' |'
        report.append(row)
        rows_out.append({'architecture': label, 'parameters': rows[0]['parameters'], 'meanScore_mm': stats.mean(values), 'sampleSD_mm': stats.stdev(values)})
    report += ['', f'Predeclared winner: **{winner}**. Frozen representative candidate: candidate_force_model.mat.',
               f'SHA256 `{representative["sha256"]}`.',
               'The existing fine-tuned GRU64 research release remains unchanged. This candidate is an offline result, not deployed or controller-qualified.', '',
               '## All-seed short-path evaluation', '',
               'Scores below are per-seed recording-macro scores, then means/SD across three seeds. Short paths use the same common anchors.', '',
               '| Cohort | Architecture | One-step mean +/- SD (mm) | Half-second path mean +/- SD (mm) | Mean recording force-rate P95 (N/s) |',
               '|---|---|---:|---:|---:|']
    for cohort in ['old_internal', 'nov3_internal', 'nov2_external']:
        for key, label in names:
            rows = [r for r in short if r['cohort'] == cohort and r['variant'] == key and abs(number(r, 'horizon_s') - .5) < 1e-9]
            first = [mean([r for r in rows if r['seed'] == str(seed)], 'oneStepVectorRMSE_mm') for seed in [113, 211, 307]]
            free = [mean([r for r in rows if r['seed'] == str(seed)], 'pathVectorRMSE_mm') for seed in [113, 211, 307]]
            rate = mean(rows, 'forceRateP95_Nps')
            report.append(f'| {cohort} | {label} | {stats.mean(first):.5f} +/- {stats.stdev(first):.5f} | {stats.mean(free):.5f} +/- {stats.stdev(free):.5f} | {rate:.3f} |')
    report += ['', '## Whole-record one-step evaluation: fixed seed113', '',
               'All eligible anchors, equal means across five recordings; length-vector RMSE, not single-coordinate error.', '',
               '| Cohort | 1x64 | 1x128 | 2x64 |', '|---|---:|---:|---:|']
    for cohort in ['nov2_external', 'nov3_selection', 'nov3_train']:
        values = [mean([r for r in one if r['cohort'] == cohort and r['variant'] == key], 'lengthVectorRMSE_mm') for key, _ in names]
        report.append('| ' + cohort + ' | ' + ' | '.join(f'{v:.5f}' for v in values) + ' |')
    report += ['', 'Nov3 TRAIN is diagnostic only; these weights were not optimized on Nov3 TRAIN, but it is not the declared whole-trial validation subset.', '',
               '## Long recorded-input free predictions: fixed seed113', '',
               'One observed initialization plus causal observed recurrent prefix; thereafter predicted q/v only. Future recorded pressure/time known.',
               'Missing observed truth stays NaN/unscored. All Nov3 comparison windows are the existing held-out whole trials, including helix trial2.', '',
               '| Day | Motion | Channel | Trial | 1x64 | 1x128 | 2x64 |', '|---|---|---:|---:|---:|---:|---:|']
    for r in motion:
        if r['variant'] != 'gru64':
            continue
        matched = [next(x for x in motion if all(x[k] == r[k] for k in ['day', 'runId', 'channel', 'trial']) and x['variant'] == key) for key, _ in names]
        vals = [f'{number(x,"lengthVectorRMSE_mm"):.5f}' if x['complete'] in ['1', 'true', 'True'] else 'FAILED' for x in matched]
        report.append('| Nov' + r['day'] + ' | ' + r['runId'] + ' | ' + r['channel'] + ' | ' + r['trial'] + ' | ' + ' | '.join(vals) + ' |')
    circles = {}
    for day in ['2', '3']:
        circles[day] = {key: number(next(r for r in motion if r['day'] == day and r['runId'] == 'circleSweep' and r['variant'] == key), 'lengthVectorRMSE_mm') for key, _ in names}
    baseline_windows = [r for r in motion if r['variant'] == 'gru64']
    wins = sum(number(next(x for x in motion if x['variant'] == 'gru128' and all(x[k] == r[k] for k in ['day','runId','channel','trial'])), 'lengthVectorRMSE_mm') < number(r, 'lengthVectorRMSE_mm') for r in baseline_windows)
    tradeoff = (f'GRU128 improves {wins}/10 fixed long examples relative to GRU64, while others worsen; these correlated examples are not independent statistical trials. '
                f'Circle free RMSE worsens on Nov2 {circles["2"]["gru64"]:.3f}->{circles["2"]["gru128"]:.3f}mm and Nov3 {circles["3"]["gru64"]:.3f}->{circles["3"]["gru128"]:.3f}mm. '
                'It improves several sine/helix examples. The short-horizon winner is therefore a candidate for measured-state model augmentation, with mixed long-horizon accuracy. '
                'Retain the current GRU64 reference; do not infer improved physical tracking from short prediction scores.')
    report += ['', '**Interpretation:** ' + tradeoff, '', '## Bounded inverse self-consistency: fixed seed113', '',
               'Each model generates its own target force at known logged pressure, then solves from a displaced initial pressure guess.',
               'This checks numerical invertibility/state handling, not physical force accuracy, pressure correction quality or closed-loop stability.',
               'Offline command bounds0..3bar; slew6bar/s at assumed.05s; fixed mean logged pressure. Candidate evaluations keep immutable pre-current state; selected state committed once.', '',
               '| Architecture | Feasible cases | Max force-equation error (N) | Median (ms) | P95 (ms) | >50ms |',
               '|---|---:|---:|---:|---:|---:|']
    for key, label in names:
        rows = read(CAMP / f'inverse_{key}' / 'inverse_cases.csv')
        times = sorted(number(r, 'solve_ms') for r in rows)
        p95 = number(read(CAMP / f'inverse_{key}' / 'summary.csv')[0], 'latencyP95_ms')
        report.append(f'| {label} | {sum(r["knownTargetFeasible"] in ["1","true","True"] for r in rows)}/{len(rows)} | {max(number(r,"forceError_N") for r in rows):.6g} | {stats.median(times):.3f} | {p95:.3f} | {sum(t>50 for t in times)} |')
    report += ['', 'Timing is this workstation only, excludes acquisition/control/communication, includes any cold/JIT outliers and is not a real-time deadline guarantee.',
               'Force-rate measurements quantify predicted residual variation, not measured mechanical jerk or physical pressure smoothness.', '',
               '## Retained artifacts and limitations', '',
               'capacity_overview and six circle/sine/helix figures have editable FIG, vector PDF and300dpiPNG copies.',
               'Predictions, CSV metrics, trial-window metadata, checkpoints, source snapshots, hashes and numerical verification retained.',
               'See VERIFICATION.md. Existing models, raw data, normalizer, controllers, Simulink files and user edits preserved.',
               'Three optimization seeds do not represent three independent hardware sessions. Architecture ranking is budget/data/protocol-specific.',
               'A fresh unseen day/condition and matched PAC/PID FB controller tests remain necessary; no hardware, deployment, migration, commit or push performed.', '']
    (CAMP / 'REPORT.md').write_text('\n'.join(report), encoding='utf-8')
    with (CAMP / 'architecture_summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows_out[0]));writer.writeheader();writer.writerows(rows_out)
    card = '\n'.join(['# Capacity candidate model card', '', f'Architecture selected: {winner}; representative seed113.',
                      f'File: candidate_force_model.mat; SHA256 `{representative["sha256"]}`.',
                      'Outputs two generalized residual forces in N on RHS of the fixed nominal dynamics; not three pressures.',
                      'Use sgCell with cfg.widths and separate hidden states per layer. Existing normalizer unchanged.',
                      'Pressure-dependent residual requires bounded implicit inverse and immutable candidate-evaluation state.',
                      'Old TRAIN only; old/Nov3 selection; Nov2 evaluation only, all previously inspected.',
                      tradeoff,
                      'This is an offline research candidate. Existing released GRU64 retained. No controller/hardware qualification.', ''])
    (CAMP / 'MODEL_CARD.md').write_text(card, encoding='utf-8')
    verification = CAMP / 'VERIFICATION.md'
    text = verification.read_text(encoding='utf-8')
    qa = ('\n## Artifact review\n\nAll seven PNG figure bundles visually inspected: readable axes/legends, consistent architecture colors, no missing model traces. FIG/PDF copies retained.\n'
          'An overview-only CSV text-column import error was repaired using explicit string imports and a render-only retry. No fits, predictions, numerical metrics or checkpoint choices changed; initial failure and successful retry remain in the evaluation diary.\n')
    if '\n## Artifact review\n' in text:text=text.split('\n## Artifact review\n')[0]
    verification.write_text(text + qa, encoding='utf-8')
    title = '2026-11-03 - Three-seed GRU capacity comparison completed'
    table = '\n'.join(report[report.index('| Architecture | Parameters | Mean score (mm) | Sample SD | Seed113 | Seed211 | Seed307 |'):report.index('## All-seed short-path evaluation')-1])
    entry = ('User authorized quick offline1x64/1x128/2x64 comparison. Three seeds113/211/307, matched original1200-update stage plus1200 old-only continuation updates; unchanged physics/loss/masks/normalizer.\n'
             'Six new stage1 and eight new continuation fits; exact existing64-unit fits reused by hash. Old TRAIN gradients only; old/Nov3 holdouts select checkpoints and architecture; Nov2 evaluation excluded from selection.\n\n'
             + table + '\n\n'
             'Representative seed113 fixed in advance, not chosen by best seed. Current fine-tuned GRU64 release retained; separate frozen candidate archived.\n'
             + tradeoff + '\n'
             'All nine cases have common-anchor evaluation; three representatives have whole-record one-step, fixed long motion predictions, and208-case bounded numerical inverse checks.\n'
             'Report: residual_dynamics_learning/force_model/capacity_comparison/campaign_20261103_v1/REPORT.md.\n'
             'Figures FIG/PDF/PNG, predictions/checkpoints/source/hashes and verification retained. No physical force/closed-loop qualification; fresh confirmation and PAC/PID FB integration pending. No hardware/SLX/migration/Git changes.')
    for name in ['PROJECT_STATE.md', 'RESIDUAL_DYNAMICS_RESULTS.md', 'VALIDATION_RESULTS.md']:
        update_doc(ROOT / 'project_docs' / name, title, entry)
    figures = ('Seven capacity comparison figure bundles at residual_dynamics_learning/force_model/capacity_comparison/campaign_20261103_v1:\n'
               'capacity_overview plus motion_comparison/nov[2,3]_[circle,sine,helix]_before_after. Each FIG/PDF/300dpiPNG.\n'
               'Architecture legends1x64/1x128/2x64; fixed seed113. Long predictions are known-command free runs, not controller tracking.\n'
               'All Nov3 windows held-out, missing truth unfilled. MAT/CSV and source snapshots retained; all negative results preserved.')
    update_doc(ROOT / 'project_docs' / 'FIGURE_INDEX.md', title, figures)
    for name in ['training_console.txt', 'evaluation_console.txt']:
        if (SRC / name).exists():shutil.copy2(SRC / name, CAMP / name)
    shutil.copy2(Path(__file__), CAMP / 'evaluation_source' / Path(__file__).name)
    artifact = CAMP / 'artifact_hashes.csv'
    paths = sorted(p for p in CAMP.rglob('*') if p.is_file() and p != artifact)
    with artifact.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream);writer.writerow(['relativePath', 'bytes', 'sha256'])
        for path in paths:writer.writerow([str(path.relative_to(CAMP)), path.stat().st_size, sha(path)])
    print(f'Report/logs updated; {len(paths)} artifacts hashed. Winner {winner}. Current release preserved.')


if __name__ == '__main__':
    main()
