# Force-GRU fine-tuning

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
