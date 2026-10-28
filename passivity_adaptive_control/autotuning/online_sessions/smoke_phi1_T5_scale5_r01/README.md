# Operator-monitored smoke test only

Model: root PAC_autotuning_online_smoke.slx. Original and ten-candidate models
are preserved. One batch, three 5 s circles, phi 1 rad, no payload assumption
unchanged from the current campaign. Score last two circles. Gamma fixed at
[.0515;.0005999;1.49999]; lower and upper bounds equal baseline, budget one.
No optimization or gain exploration occurs. Hysteresis scale/rate/upper 5/10/50.

1. Run this folder's start_worker.m in a SECOND MATLAB process. It only listens
   on localhost. Do not run another online/synthetic worker on the same ports.
2. Run open_model.m in the model process. Review the setup and keep physical
   pressure kill available. Pressure enable defaults Off: set it On BEFORE Run
   when ready. Leaving it Off causes terminal shutdown, not an experiment.
3. Start the model manually and monitor the entire test. Expect startup/trim,
   bend and rotation ramps, three circles, withdrawal, command pressure slew
   to zero and automatic simulation stop, approximately one minute if holds
   are achieved promptly. Use physical kill if motion becomes unacceptable.
4. Save the full out.logsout as online_smoke_00g_phi1_T5_3cycles_r01.mat in this
   folder, provided it was actually unloaded; otherwise label actual payload.
   Keep worker_log.mat and batch_001.mat. Report any abort reason and do not
   restart/overwrite this session; a repeat needs a fresh session/configuration.

Expected normal finish: onlineReason=23 (single-batch budget exhausted).
This is intentional, even if accuracy targets are missed. Smoke acceptance is
correct plumbing/sequence, constant gamma, acceptable monitored motion, no
faults/limits, valid logging, and proper pressure-command shutdown. A short
three-circle test does not demonstrate settled adaptation or unattended safety.
Review its saved logs before enabling a multi-candidate search.
