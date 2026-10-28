# Smoke retry r02: pacing corrected, old attempt preserved

The r01 arm model did not start: Desktop Real-Time rejects Simulink pacing.
The old worker exited with zero completed batches/zero packets and preserved
worker_log.mat. Its overwrite guard is intentional; do not delete that log.

Use PAC_autotuning_online_smoke_r02.slx and THIS folder's start_worker.m.
The new model has EnablePacing=off and a new session identity. Worker waits up
to 30 minutes for the first packet so preflight/compilation can finish. Do not
reuse a still-open r01 model with this worker: their session identities differ.
open_model.m opens only the correct r02 model and adds paths.

One fixed-gain batch: three 5 s circles, phi1, final two scored; gamma
[.0515;.0005999;1.49999], fixed equal bounds, budget one, no optimization.
Set pressure enable On before Run when ready; monitor with physical kill
available. Expected normal reason23 budget exhausted, followed by pressure
command slew-to-zero and automatic stop. No hardware model was run to test
this correction. Start manually after review.

Save full out in this folder as online_smoke_00g_phi1_T5_3cycles_r02.mat if
actually unloaded (otherwise label actual payload), and retain worker logs.
Review results before gain exploration. Any further attempt needs a fresh
session if this worker has already exited and written its log.
