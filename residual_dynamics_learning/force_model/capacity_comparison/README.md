# GRU capacity comparison

Offline comparison requested November3: 1x64, 1x128, 2x64, seeds113/211/307.
Each receives the same original1200-update old-only stage (Adam3e-4), followed
by1200 extra old-only updates (fresh Adam1e-4), preserving the existing
normalizer, physics, loss, masks, history and whole-trial splits. Exact matching
previous64-unit fits are reused by hash. No interpolation or new TRAIN data.

Stage1 checkpoint selection uses older holdouts. Stage2 selection uses equal
old/Nov3 cohort means, .7oneStep+.3halfSecond path per cohort, and the same
10% old regression guard relative to each stage1 checkpoint. Architecture
ranking uses mean composite across three seeds; representative seed113 fixed
in advance. Different parameter counts: this is a practical width/depth test,
not a parameter-matched proof. Same update budget does not guarantee convergence.

`acRun(pwd)` creates fresh campaign_20261103_v1 and freezes all nine fits before
Nov2 loading. `acFinish(pwd)` evaluates every seed on common short-path anchors,
the three fixed representatives on whole-record one-step and fixed long
motion windows, and all three on bounded offline inverse self-consistency.
This does not test force ground truth or closed-loop tracking. Nov2/Nov3 were
previously inspected; a fresh experiment is still needed for confirmation.

Current released model is preserved. A candidate is archived separately;
controller integration, hardware, SLX, migration and Git actions are excluded.
