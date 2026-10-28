# r03 fixed-gain smoke: validity-alignment correction

Use this folder's start_worker.m in the second MATLAB process and open_model.m
in the model process. Open helper sets a short build cache; pacing remains off.
Do not use an old loaded r02 model with the r03 worker. Enable operator pressure
before Run only when ready, and monitor with physical kill available.

Same three5s circles, phi1, gamma.0515/.0005999/1.49999, hysteresis scale/rate/
upper5/10/50, one-batch budget, no gain optimization. r02 record/log preserved.
Corrected raw-validity grace so isolated dropouts are not joined by delayed
PAC status. Real sustained invalidity, actual movement during application,
PAC faults, limiting and heartbeat watchdog remain guarded.

Expected normal terminal reason23 and completed batch1. Save complete out as
online_smoke_00g_phi1_T5_3cycles_r03.mat if unloaded (otherwise label actual
payload) in this folder, including aborted recordings. Review before search.
The old r02 run was aborted at34s after0.88circles; it was not a smoke pass.
