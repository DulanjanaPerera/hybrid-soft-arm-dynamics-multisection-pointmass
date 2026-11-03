# Fine-tuned force-GRU64 research checkpoint

Selected arm: old_only;1200 additional weight updates from original seed113.
Selection used old/Nov3 internal holdouts before Nov2 reading. New weights
were optimized on October19/27 TRAIN only; a separately fitted old+Nov3 replay
candidate is retained in ../mixed_replay/force_model.mat.

SHA256: 1c38edc36ba6ff4a41da8f974a88ca1c5ea775655426da91124a35b255a9c870

One GRU layer64, two force outputs in N on RHS of
Mddq+(C+D0)dq+K0q+bound+G=B*p+r_tau. Fixed nominal K0=650*[2 1;1 2], D0=40I,
same corrected bounds/physics. Input is q2,q3,v2,v3,p1,p2,p3,previous physical dt,
normalized with the original older-TRAIN mean/std. Recurrent state evolves;
weights stay fixed during inference. Use sgCell/sgFlow, not legacy forceResidualCell.
The active network width is cfg.widths; legacy cfg.hiddenUnits is not its runtime.

Nov2 one-step macro vector RMSE 0.36143mm versus original
0.53588mm. Full-circle example Nov2
0.9440->0.6229mm,
Nov3 heldout 1.2636->1.0047mm.
Some new-day long windows regress; see ../REPORT.md and complete comparisons.

Normalized inputs, filters, nominal dynamics and trained force bound remain
unchanged. Current-pressure dependence requires a bounded inverse calculation.
PAC adaptation/residual consistency and physical tests remain pending. This is
a frozen research candidate, not an installed PAC/PID FB model or hardware release.
