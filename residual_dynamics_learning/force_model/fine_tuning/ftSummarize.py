"""Write reports from frozen fine-tune artifacts; never select on external data."""
from pathlib import Path
import csv, statistics as st, json, hashlib

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CAMP=HERE/'campaign_20261103_v1'
def read(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def num(r,k):return float(r[k])
def avg(rows,k):return st.mean(num(r,k) for r in rows)
def gain(before,after):return 100*(1-after/before)
assert (CAMP/'VERIFICATION.md').exists()
R=read(CAMP/'internal_results.csv');selected=min(R,key=lambda r:num(r,'compositeScore_mm'));name=selected['arm']
O=read(CAMP/'evaluation/one_step.csv');P=read(CAMP/'evaluation/short_paths.csv')
L=read(CAMP/'motion_comparison/scores.csv');W=read(CAMP/'motion_comparison/windows.csv')
splits=read(CAMP/'data/nov3_trial_split.csv');coverage=read(CAMP/'data/coverage.csv')
inverse=json.loads((CAMP/'inverse/summary.json').read_text())
arms=['frozen','old_only','mixed_replay'];labels={'frozen':'Original frozen','old_only':'Old-data continuation','mixed_replay':'Old+Nov3 replay'}
def one(cohort,variant):return avg([r for r in O if r['cohort']==cohort and r['variant']==variant],'lengthVectorRMSE_mm')
def long(day,id,ch,variant):return num(next(r for r in L if r['day']==str(day) and r['runId']==id and r['channel']==str(ch) and r['variant']==variant),'lengthVectorRMSE_mm')
baseline=R[0];initial_hash=baseline['sha256'];model_hash=selected['sha256']
train_chunks=sum(int(r['chunks']) for r in coverage if r['group']=='2' and r['partition']=='train')
held_trials=sum(r['partition']=='selection' for r in splits)
lines=['# Fine-tuned generalized-force GRU64: completed offline comparison','',
       '2026-11-03. User authorized fine-tuning for eventual PAC and PID feedback',
       'model augmentation. Two weight continuations completed, with a frozen',
       'reference. No controller or hardware changes. November2 was read by this',
       'experiment only after all weights and checkpoint choices were frozen.','',
       f'**Selected research checkpoint: {labels[name]}**, using the predeclared',
       'equal old/Nov3 internal score. File: release/force_model.mat.',
       f'SHA256: `{model_hash}`.',
       'The selected weights were optimized on older TRAIN data; November3 internal',
       'holdouts participated in checkpoint selection. The independently fitted',
       'old+Nov3 replay candidate is retained in mixed_replay/force_model.mat.',
       'Selection preceded external November2 and long-window results; it was not',
       'changed afterward. Neither candidate replaces the earlier release files.','',
       '## Protocol and data roles','',
       '- Exact warm start from streaming/release_20261103_v1 GRU64/B seed113,',
       f'  SHA256 `{initial_hash}`. No random network reinitialization.',
       '- One recurrent layer64 units, same .5s training rollout and chunk, .5s',
       '  burn-in, same physical-time guards, nominal mechanics/integrator,80N',
       '  force bound and length/force regularization. Both branches of training',
       '  retain the checked measured-state/independent-free-copy contract.',
       '- Original October19/27 train/selection assignments and59,630-TRAIN-row',
       '  normalizer unchanged. Old128 selection anchors retained exactly.',
       f'- November3: allfive supplied recordings audited; {len(splits)-held_trials}',
       f'  complete TRAIN trials and {held_trials} whole held-out channel/trials.',
       f'  {train_chunks} new TRAIN chunks,80 new selection anchors. Deterministic',
       '  middle listed trial per channel, the same split rule as the older data.',
       '  Holdouts: circle3; sine1/3/5; helix2; steps3/9/15; combined2.',
       '- Adam restarted because the original released model lacks Adam state.',
       '  Learning rate1e-4 versus original3e-4; gradient clip1; batch32; seed113;',
       '  exactly1200 additional updates per arm, checkpoints every150. New',
       '  checkpoints retain Adam averages/squares, RNG, stream state and queues.',
       '- Old-only arm:32 old streams. Mixed arm:16 old +16 Nov3 streams; equal',
       '  stream counts, not a claim of equal frame counts with variable dt.',
       '- No interpolation, filtering, objective/architecture changes or Nov2',
       '  fitting. Nine Nov3 holdouts are internal previously inspected validation;',
       '  Nov2 is previously inspected external evaluation, not blind confirmation.','',
       'Selection: .5*old+.5*Nov3 cohort score. Within each cohort, equal recording',
       'means of .7 one-step + .3 half-second free-path vector RMSE. Older score',
       'must remain <=110% of the initial old score; the initial checkpoint is',
       'eligible so a worse fine-tune cannot be promoted by default. Both trained',
       'candidates improve BOTH initial cohort scores. They trade off which cohort',
       'improves more; only one optimization seed was run, not a significance test.','',
       '| Candidate | Extra updates selected | Old internal score (mm) | Nov3 internal score (mm) | Composite (mm) |',
       '|---|---:|---:|---:|---:|']
for r in R:lines.append(f'| {labels[r["arm"]]} | {r["extraUpdatesSelected"]} | {num(r,"oldScore_mm"):.5f} | {num(r,"nov3Score_mm"):.5f} | {num(r,"compositeScore_mm"):.5f} |')
lines+=['',f'Selected composite improves {gain(num(baseline,"compositeScore_mm"),num(selected,"compositeScore_mm")):.1f}%.',
        'Nov3 replay gives the largest gain on NEW internal score, while old-data',
        'continuation gives the larger old-domain gain and slightly lower composite.',
        'This shows that additional optimization mattered; extra data was not an',
        'automatic improvement of every metric under this fixed update budget.','',
        '## Whole-record measured-state one-step results','',
        'Equal-recording macro q2/q3 vector RMSE in mm. Same fixed initial seed113',
        'checkpoint for before/after, not the earlier three-seed mean0.563mm.','',
        '| Cohort | Original | Old-data continuation | Old+Nov3 replay |',
        '|---|---:|---:|---:|']
for cohort in ['nov2_external','nov3_selection','nov3_train']:
    lines.append('| '+cohort+' | '+' | '.join(f'{one(cohort,a):.5f}' for a in arms)+' |')
lines+=['',f'Selected Nov2 one-step gain: {gain(one("nov2_external","frozen"),one("nov2_external",name)):.1f}%.',
        'Nov3 TRAIN metrics are explicitly training diagnostics. Original CV remains',
        f'better one-step: Nov2 {one("nov2_external","constant_velocity"):.5f}mm and',
        f'Nov3 heldout {one("nov3_selection","constant_velocity"):.5f}mm. These are',
        'prediction errors with measured-state refresh, not hardware tracking errors.',
        'Per-coordinate RMSE/p95, numerical guards and force-bound occupancy are',
        'retained in evaluation CSV/MAT files. No failed predictions were discarded.','',
        '## Short free-prediction results','',
        'Each cell is a macro mean of recording-level path vector RMSE, mm.','',
        '| Cohort/horizon | Original | Old-data continuation | Old+Nov3 replay |',
        '|---|---:|---:|---:|']
for cohort in ['old_internal','nov3_internal','nov2_external']:
    for hz in [.5,1.]:
        vals=[avg([r for r in P if r['cohort']==cohort and r['variant']==a and abs(num(r,'horizon_s')-hz)<1e-8],'pathVectorRMSE_mm') for a in arms]
        lines.append(f'| {cohort}, {hz:g}s | '+' | '.join(f'{v:.5f}' for v in vals)+' |')
lines+=['','## Complete motion-window diagnostics','',
        'Free propagation has one measured initialization and observed recurrent',
        'prefix, then integrated predicted q/v and known recorded pressures/time.',
        'No future measured motion enters the prediction. Missing truth remains',
        'NaN/unscored. All30 candidate/window paths complete the numerical guard.',
        'These long diagnostic windows exceed the .5s training horizon.','',
        '**Every Nov3 window here belongs to a whole held-out trial.** Helix uses',
        'trial2 (76.6-134.95s) instead of the prior illustration of TRAIN trial1.',
        'Therefore compare its original4.0294mm to its updated values, not to the',
        'previous trial1 value2.8787mm. Nov2 windows are unchanged from prior plots.',
        'PCC paths are reconstructed from q, not independent measured sensor XYZ.','',
        '| Day/window | Original (mm) | Old-data continuation | Old+Nov3 replay |',
        '|---|---:|---:|---:|']
for w in W:
    day=int(w['day']);id=w['runId'];ch=int(w['channel']);label=f'Nov{day} {id}'+(f' P{ch}' if ch else '')+f' trial{w["trial"]}'
    lines.append('| '+label+' | '+' | '.join(f'{long(day,id,ch,a):.4f}' for a in arms)+' |')
lines+=['','Selected continuation improves the Nov2 circle and allthree sine windows',
        'and its helix example. The Nov3 circle and sineP1 improve; sineP2/P3 and',
        'held-out helix2 worsen slightly. Mixed replay improves new-day one-step',
        'prediction most but worsens several long rollouts, including both circles',
        'and sineP2. These regressions are retained, not hidden by an average score.',
        'Training/selection focused on measured-state and short free prediction;',
        'improvement does not extend automatically to all long trajectories.','',
        '## Offline inverse and verification','',
        f'Selected model solves {inverse["knownFeasible"]}/{inverse["cases"]} self-generated',
        f'feasible targets. Max force-equation residual {inverse["forceErrorMax_N"]:.5g}N;',
        f'max pressure-bound violation {inverse["maxBoundViolation_bar"]:.5g}bar.',
        f'Median/p95 observed inverse+selected-GRU time {inverse["latencyMedian_ms"]:.3f}/',
        f'{inverse["latencyP95_ms"]:.3f}ms; {round(inverse["fractionOver50ms"]*inverse["cases"])}/{inverse["cases"]} exceed50ms.',
        'Cold optimizer calls are retained. Identical network architecture and',
        'different machine/JIT load mean these timings are not a claimed training-',
        'induced speedup or qualification for the full50ms hardware cycle. Bounds',
        'are the offline RL0..3bar domain, not new PAC runtime settings. Trial',
        'pressures reuse fixed pre-current hidden state and commit once per sample.',
        'This is model self-consistency, not measured physical force accuracy.','',
        (CAMP/'VERIFICATION.md').read_text(encoding='utf-8'),
        '## Retained artifacts and downstream scope','',
        '- PREDECLARED_PROTOCOL.mat/PROTOCOL.md and data/: source/whole-trial masks,',
        ' 208 selection anchors, unchanged normalization, raw/prepared hashes.',
        '- old_only/ and mixed_replay/: all checkpoints, optimizer/stream states,',
        ' training curves, internal selection tables and frozen weight files.',
        '- release/: selected research weights and model card; previous releases',
        ' untouched. Candidate alternative remains in mixed_replay/force_model.mat.',
        '- evaluation/: all three models, Nov2 external, separate Nov3 TRAIN and',
        ' heldout results; one-step and .1/.25/.5/1s predictions/forces/CSV/MAT.',
        '- motion_comparison/: fixed windows, predictions/truth/commands/time,',
        ' scores and six before/after figure bundles. fine_tune_overview is the',
        ' seventh bundle, each editable FIG/PDF/300dpiPNG; all regressions retained.',
        '- inverse/: all208 cases, pressure/state/latency/Jacobian diagnostics.',
        '- training_source/evaluation_source and hashes, verification and diaries.','',
        'The downstream goal remains matched PAC and PID FB model-off/on physical',
        'comparisons. Controller gains/laws, SLX and hardware were not changed.',
        'Resolve residual definition versus PAC adaptation, pressure/common-pressure',
        'coverage and measured-velocity/runtime conventions before integration.',
        'No explicit payload estimator, physical parameter identification, passivity',
        'proof or closed-loop tracking improvement is established by this fit.',
        'No migration, Git commit or push. Future confirmation requires fresh data.','']
report='\n'.join(lines);(CAMP/'REPORT.md').write_text(report,encoding='utf-8')
card=f'''# Fine-tuned force-GRU64 research checkpoint

Selected arm: {name};1200 additional weight updates from original seed113.
Selection used old/Nov3 internal holdouts before Nov2 reading. New weights
were optimized on October19/27 TRAIN only; a separately fitted old+Nov3 replay
candidate is retained in ../mixed_replay/force_model.mat.

SHA256: {model_hash}

One GRU layer64, two force outputs in N on RHS of
Mddq+(C+D0)dq+K0q+bound+G=B*p+r_tau. Fixed nominal K0=650*[2 1;1 2], D0=40I,
same corrected bounds/physics. Input is q2,q3,v2,v3,p1,p2,p3,previous physical dt,
normalized with the original older-TRAIN mean/std. Recurrent state evolves;
weights stay fixed during inference. Use sgCell/sgFlow, not legacy forceResidualCell.
The active network width is cfg.widths; legacy cfg.hiddenUnits is not its runtime.

Nov2 one-step macro vector RMSE {one('nov2_external',name):.5f}mm versus original
{one('nov2_external','frozen'):.5f}mm. Full-circle example Nov2
{long(2,'circleSweep',0,'frozen'):.4f}->{long(2,'circleSweep',0,name):.4f}mm,
Nov3 heldout {long(3,'circleSweep',0,'frozen'):.4f}->{long(3,'circleSweep',0,name):.4f}mm.
Some new-day long windows regress; see ../REPORT.md and complete comparisons.

Normalized inputs, filters, nominal dynamics and trained force bound remain
unchanged. Current-pressure dependence requires a bounded inverse calculation.
PAC adaptation/residual consistency and physical tests remain pending. This is
a frozen research candidate, not an installed PAC/PID FB model or hardware release.
'''
(CAMP/'release/MODEL_CARD.md').write_text(card,encoding='utf-8')
(HERE/'README.md').write_text('''# Force-GRU fine-tuning

Completed campaign: campaign_20261103_v1/REPORT.md.
Two exact-weight continuations, old-only and balanced old/Nov3 replay. The
predeclared guarded equal-cohort score selects old-only, with replay retained
as an alternative. No old trainer/source/dataset/release was overwritten.

MATLAB entry points: ftRun(pwd), then after freeze ftEvaluate(pwd),
ftInverseCheck(pwd), ftMotionCompare(pwd), ftVerify(pwd). Fresh versioned output
folders required. ftMotionCompare(pwd,true) renders saved results only.
ftSummarize.py reports saved results; it never selects using Nov2 outcomes.

Freeze before external evaluation; declare every data role before fitting.
No hardware, controller/SLX edits, migration or Git actions are included.
''',encoding='utf-8')
entry=f'''## 2026-11-03 - Force GRU64 fine-tuning completed

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
campaign_20261103_v1/release/force_model.mat, SHA256 {model_hash}.
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

'''
docs=ROOT/'project_docs';marker=entry.splitlines()[0]
for filename in ['PROJECT_STATE.md','RESIDUAL_DYNAMICS_RESULTS.md','VALIDATION_RESULTS.md']:
    p=docs/filename;text=p.read_text(encoding='utf-8')
    if marker in text:
        start=text.index(marker);end=text.find('\n## ',start+len(marker))
        text=text[:start]+entry+(text[end+1:] if end>=0 else '')
    elif filename=='PROJECT_STATE.md':
        text=text.replace('# Current project state\n\n','# Current project state\n\n'+entry,1)
    else:
        title,rest=text.split('\n',1);text=title+'\n\n'+entry+rest.lstrip('\n')
    if filename=='PROJECT_STATE.md':
        lines=text.splitlines();lines[1]='Current Milestone: Force-GRU64 fine-tuning complete; guarded internal selection favors old-data continuation, with mixed Nov3 candidate retained. Prediction gains verified; PAC/PID FB integration remains pending.';text='\n'.join(lines)+'\n'
    p.write_text(text,encoding='utf-8')
p=docs/'FIGURE_INDEX.md';text=p.read_text(encoding='utf-8');needle='## Other retained sources\n\n'
figure_entry='''- 2026-11-03 force-GRU fine-tuning comparison:
  residual_dynamics_learning/force_model/fine_tuning/campaign_20261103_v1.
  fine_tune_overview plus motion_comparison/nov2_circle/sine/helix_before_after
  and nov3_circle/sine/helix_before_after. Seven FIG/PDF/300dpiPNG bundles;
  full trajectories/forces/truth/masks/pressure/time, score/windows CSVs,
  checkpoints/optimizer states, source/model hashes retained. Original, old-data
  continuation and old+Nov3 replay, including regressions. Nov3 motion windows
  all heldout; helix uses trial2, unlike prior TRAIN-trial1 illustration. Frozen
  weights, measured initialization then no state resets; no physical tracking claim.

'''
if figure_entry not in text:text=text.replace(needle,needle+figure_entry,1)
p.write_text(text,encoding='utf-8')
p=docs/'RESIDUAL_CONTROLLER_COMPARISON_PLAN.md';text=p.read_text(encoding='utf-8')
text=text.replace('No continuation fit has begun. Freeze/qualify that model before deciding which\ncheckpoint to integrate. The shared inverse and PAC consistency work above\nremain pending. All arm experiments are operator-run; no autonomous hardware.',
'''Fine-tuning is now complete: fine_tuning/campaign_20261103_v1/REPORT.md under
force_model. The equal-cohort internal score selects old-data continuation;
old+Nov3 replay is retained. Old-only gives larger Nov2/circle gains, while
replay gives lower Nov3 one-step error but worsens several long windows. No
controller integration follows automatically from this checkpoint selection.
Shared inverse/PAC consistency, operating-region and physical qualification
remain pending. All arm experiments are operator-run; no autonomous hardware.''')
p.write_text(text,encoding='utf-8')
p=docs/'RESIDUAL_DYNAMICS_LEARNING_PLAN.md';text=p.read_text(encoding='utf-8')
learning_entry='''\n\n## 2026-11-03 - Fine-tuning execution completed

The authorized comparison is complete in force_model/fine_tuning/
campaign_20261103_v1/REPORT.md under residual_dynamics_learning. Exact weight
warm starts, fixed loss/physics/normalizer, fresh Adam1e-4,1200extra updates per
arm. Original split and nine Nov3 whole-trial holdouts preserved. Predeclared
guarded equal-old/new internal score selects old-only continuation; mixed Nov3
candidate retained. Mixed reduces new one-step error more but regresses some
long free trajectories; selected continuation improves circles and allfive
Nov2 long examples, with some Nov3 long regressions. Keep allresults.

No Nov2-informed refit/reselection. No interpolation, capacity, preprocessing,
controller or hardware changes. This completes the prior proposed continuation
experiment; it does not implement PAC/PID FB integration. The next model study
could investigate the short-loss/long-rollout trade-off using development data,
with a fresh confirmation cohort for resulting changes. Do not treat improved
prediction RMSE or inverse self-consistency as proven physical tracking gains.
'''
if '## 2026-11-03 - Fine-tuning execution completed' not in text:text+=learning_entry
p.write_text(text,encoding='utf-8')
print('Selected model:',name,model_hash)
print('Nov2 one-step: %.6f -> %.6f mm'%(one('nov2_external','frozen'),one('nov2_external',name)))
print('Reports/model card/project logs written. Both candidate weights and all regressions retained.')
