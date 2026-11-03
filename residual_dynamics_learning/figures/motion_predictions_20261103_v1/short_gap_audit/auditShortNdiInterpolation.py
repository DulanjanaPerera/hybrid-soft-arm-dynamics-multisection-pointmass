"""Offline short-gap feasibility preview; never edits raw/normalized training data."""
from pathlib import Path
import csv, math, statistics, collections, hashlib, json

ROOT=Path(__file__).resolve().parent.parent
FOLDER=ROOT/'residual_dynamics_learning/figures/motion_predictions_20261103_v1'
OUT=FOLDER/'short_gap_audit'
OUT.mkdir(exist_ok=False)
MAX_SPAN=.12

def value(row,name):
    return float(row.get(name,'nan'))

def quantile(values,p):
    a=sorted(values);x=(len(a)-1)*p;i=int(x)
    return a[i]+(a[min(i+1,len(a)-1)]-a[i])*(x-i)

def write(name,fields,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def blocks(mask):
    start=None
    for i,flag in enumerate(mask+[False]):
        if flag and start is None:start=i
        elif not flag and start is not None:
            yield start,i-1
            start=None

summaries=[];gap_rows=[];test_rows=[];sources=[]
for name in ['circleSweep','SineWave','helix','StepWave','combinePressure']:
    file=FOLDER/'exports'/name/'observations.csv'
    r=list(csv.DictReader(file.open(newline='',encoding='utf-8-sig')))
    sources.append(dict(source=str(file),sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
    n=len(r);t=[value(x,'sensor_poll_start_host_s') for x in r]
    q=[[value(x,'length_inferred_m_2'),value(x,'length_inferred_m_3')] for x in r]
    original=[a[:] for a in q]
    channel=[value(x,'channelId') if 'channelId' in x else 0 for x in r]
    trial=[value(x,'trialId') for x in r];phase=[value(x,'phaseId') for x in r]
    good=[value(x,'orientation_valid')==1 and value(x,'base_orientation_valid')==1 and value(x,'fresh_sensor_frame')==1
          and all(math.isfinite(v) for v in a) for x,a in zip(r,q)]
    gaps=list(blocks([not x for x in good]));imputed=[False]*n;span_by_row={}
    for a,b in gaps:
        bracketed=a>0 and b<n-1 and good[a-1] and good[b+1]
        span=t[b+1]-t[a-1] if bracketed else float('nan')
        same=bracketed and trial[a-1]>0 and all(trial[j]==trial[a-1] and channel[j]==channel[a-1] for j in range(a-1,b+2))
        base_ok=all(value(r[j],'base_orientation_valid')==1 for j in range(a,b+1))
        input_ok=bracketed and all(0<t[j]-t[j-1]<=.1 for j in range(a,b+2))
        endpoint_flags=bracketed and all(value(r[j],s)==0 for j in [a-1,b+1]
                                       for s in ['ndi_tip_out_of_volume','ndi_tip_partial_out_of_volume','ndi_tip_sensor_broken'])
        candidate=bool(bracketed and same and base_ok and input_ok and endpoint_flags and b==a and 0<span<=MAX_SPAN)
        if candidate:
            f=(t[a]-t[a-1])/span;q[a]=[original[a-1][d]+f*(original[b+1][d]-original[a-1][d]) for d in [0,1]]
            imputed[a]=True;span_by_row[a]=span
        gap_rows.append(dict(runId=name,firstRow=a+1,lastRow=b+1,missingSamples=b-a+1,
                             startSim_s=value(r[a],'time_sim_s'),bracketingSpan_s=span,sameTrial=same,
                             baseValidInside=base_ok,oneFrameCandidate=candidate,
                             tipStatus=','.join(str(int(value(r[j],'ndi_tip_sensor_status'))) for j in range(a,b+1)),
                             tipRejectMask=','.join(str(int(value(r[j],'ndi_tip_reject_mask'))) for j in range(a,b+1))))
    preview=[]
    for i in range(n):
        preview.append(dict(sourceRow=i+1,timeSim_s=value(r[i],'time_sim_s'),physicalTime_s=t[i],
                            q2_raw_m=original[i][0],q3_raw_m=original[i][1],q2_preview_m=q[i][0],q3_preview_m=q[i][1],
                            q1_preview_m=-sum(q[i]),originalNdiValid=good[i],imputed=imputed[i],
                            bracketingSpan_s=span_by_row.get(i,float('nan')),trial=trial[i],channel=channel[i],phase=phase[i]))
    write(name+'_preview.csv',list(preview[0]),preview)
    # Hide REAL valid observations, reconstruct from surrounding REAL endpoints.
    # No truth exists for the actual dropped frames; this is a proxy test only.
    for missing in [1,2]:
        eligible=[]
        for a in range(1,n-missing):
            b=a+missing-1;rows=range(a-1,b+2)
            if not (all(good[j] for j in rows) and trial[a]>0):continue
            if not all(trial[j]==trial[a] and channel[j]==channel[a] and phase[j]==phase[a] for j in rows):continue
            if not all(0<t[j]-t[j-1]<=.1 for j in range(a,b+2)):continue
            eligible.append(a)
        selected=set()
        if eligible:
            cap=min(400,len(eligible));selected={eligible[round(i*(len(eligible)-1)/max(1,cap-1))] for i in range(cap)}
        # Include valid neighborhoods of actual losses so testing is not solely
        # uniform sampling of easy quiet holds. Samples can overlap: no CI claim.
        near={a+offset for a,b in gaps for offset in range(-10,11)}
        selected.update(near.intersection(eligible))
        e2=[];e3=[];vector=[]
        for a in sorted(selected):
            b=a+missing-1;duration=t[b+1]-t[a-1]
            for k in range(a,b+1):
                fraction=(t[k]-t[a-1])/duration
                pred=[original[a-1][d]+fraction*(original[b+1][d]-original[a-1][d]) for d in [0,1]]
                err=[1000*(pred[d]-original[k][d]) for d in [0,1]]
                test_rows.append(dict(runId=name,hiddenFrames=missing,sourceRow=k+1,phase=phase[k],trial=trial[k],
                                      channel=channel[k],bracketingSpan_s=duration,q2Error_mm=err[0],q3Error_mm=err[1],
                                      vectorError_mm=math.hypot(*err)))
                e2.append(err[0]);e3.append(err[1]);vector.append(math.hypot(*err))
        summaries.append(dict(runId=name,hiddenFrames=missing,windows=len(selected),scoredSamples=len(vector),
                              q2RMSE_mm=math.sqrt(statistics.mean(v*v for v in e2)),
                              q3RMSE_mm=math.sqrt(statistics.mean(v*v for v in e3)),
                              vectorRMSE_mm=math.sqrt(statistics.mean(v*v for v in vector)),
                              vectorP95_mm=quantile(vector,.95),maxVectorError_mm=max(vector)))

write('gap_candidates.csv',list(gap_rows[0]),gap_rows)
write('artificial_mask_errors.csv',list(test_rows[0]),test_rows)
write('reconstruction_summary.csv',list(summaries[0]),summaries)
write('sources.csv',list(sources[0]),sources)
summary=dict(actualEpisodes=len(gap_rows),actualMissingSamples=sum(x['missingSamples'] for x in gap_rows),
             oneFrameCandidates=sum(x['oneFrameCandidate'] for x in gap_rows),maximumBracketSpan_s=MAX_SPAN,
             actualMissingTruthKnown=False,usedForTrainingOrValidation=False,rawOrExporterChanged=False)
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
report=['# November3 short NDI gap audit','',
        'Offline feasibility preview only. Raw recordings, exports, validity flags,',
        'training masks, all frozen weights and November2 evaluation are unchanged.','',
        f'{summary["actualMissingSamples"]} invalid NDI samples in {summary["actualEpisodes"]} episodes.',
        f'{summary["oneFrameCandidates"]} one-frame episodes meet the conservative preview rule:',
        'one missing sample, valid clean endpoints, positive same channel/trial,',
        'base still valid, individual physical intervals <=.1s and bracketing',
        'valid-to-valid span <=.12s. Two-frame episodes remain missing.',
        'Preview interpolates q2/q3 on actual host timestamps; q1=-q2-q3.',
        'Pressure commands are never interpolated. No held-invalid value is',
        'silently reclassified as a real measurement. Each preview row has',
        'originalNdiValid and imputed flags; it has NOT entered any model pipeline.','',
        '## Reconstruction proxy test','',
        'Real valid blocks were artificially hidden, including neighborhoods of',
        'actual dropouts. Linear reconstruction uses bounding observations.',
        'Errors below compare with the hidden real measurements. They do not',
        'measure error at actual gaps, whose truth is unavailable. Overlapping',
        'test blocks are correlated; this is not an independent validation set.','',
        '| Recording | Hidden frames | Windows | Vector RMSE (mm) | p95 norm error (mm) | Max norm error (mm) |',
        '|---|---:|---:|---:|---:|---:|']
for r in summaries:
    report.append(f'| {r["runId"]} | {r["hiddenFrames"]} | {r["windows"]} | {r["vectorRMSE_mm"]:.4f} | {r["vectorP95_mm"]:.4f} | {r["maxVectorError_mm"]:.4f} |')
report += ['', '## Recommendation and next experiment','',
           'Short-gap interpolation is plausible as a separate offline training',
           'input experiment, not an unconditional replacement for validity.',
           'Keep the original measured-only target mask; do not train force or',
           'acceleration truth labels from filled lengths. Retain long-gap resets',
           'and time/pressure/calibration/trial boundaries. Preserve unchanged',
           'measured-only validation, and compare causal operation under simulated',
           'dropouts before accepting a new recurrent-state policy.','',
           'Two-sided interpolation uses the next valid sample, so it is not',
           'available at zero delay online. A deployed estimator must predict/hold',
           'causally, expose uncertainty and missingness, or explicitly accept delay.',
           'Masking only the filled targets does not remove future leakage if',
           'imputed inputs already advanced the validation hidden state.','',
           'An alternative is gap-aware physics rollout with loss only on real',
           'observations. Free predictions can already traverse a missing target',
           'using known commands/time without imputing measured q. Measurement-',
           'conditioned training needs an explicit tested missing-input/state',
           'policy; the current frozen reset convention has NOT been changed.','',
           'Do not interpolate wrapped theta near straightness. If raw rotations',
           'are later used, quaternion SLERP is preferable to blending Euler angles.',
           'Smoother appearance alone is not evidence of correct derivatives or dynamics.','',
           'References: MATLAB fillmissing supports physical SamplePoints, MaxGap',
           'and returned fill masks: https://www.mathworks.com/help/matlab/ref/fillmissing.html',
           'Quaternion SLERP: https://www.mathworks.com/help/uav/ref/quaternion.slerp.html','']
(OUT/'REPORT.md').write_text('\n'.join(report),encoding='utf-8')
print(json.dumps(summary,indent=2));print('\n'.join(report[report.index('## Reconstruction proxy test'):]))
