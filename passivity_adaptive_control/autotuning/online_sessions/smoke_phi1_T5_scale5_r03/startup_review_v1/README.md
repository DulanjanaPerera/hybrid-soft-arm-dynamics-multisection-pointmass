# R03 startup timeout: no circle was commanded

Recorded model identity PAC_autotuning_online_smoke_r03 verified. Phase0 until
50.05s, terminal withdrawal9 with reason9 startup timeout, phase10 at55.05,
recording ends56.30s. Actual online desired phi/theta remained zero, candidate
ID0: no circle, gain application or optimization occurred. Constant baseline
gamma preserved. PAC normal status3 from20s until shutdown, no latched faults.
Velocity peaks20-50s[.000187003,.000134689]m/s, below.0005 stationary threshold;
measurement/base/model flags valid. Non-heartbeat startup gates ready20.55s,
590 startup samples eligible. Therefore no accepted fresh worker heartbeat was
the missing startup condition. Final pressure commands zero.

Transport evidence at inspection: netstat showed both model port32151 and
worker port32152 owned by MATLAB backend28480, while R03's generated cache was
under PID60752. Worker lock existed without a completed worker log. A stale
model bridge retained in the MATLAB instance now running the worker is the
supported port-conflict explanation; heartbeat was not logged directly.
No physical instability or gamma conclusion follows from this stationary run.

Worker startup now clears its old model bridge (only with models stopped).
Model bridge bind failure warns once and reports explicit negative transport
status, producing runtime reason13 rather than hiding the problem behind an
eventual startup timeout. Added pacOnlineConnectionCheck and this session's
check_connection.m: model-side pressure-free probe checks port, session and
baseline. Worker probe replies do not start the post-connection10s timeout or
record a trial. Old r03 model recording preserved; use a distinct filename for
retry. No new physical run or gain changes by assistant.

Regression tested on alternate32251/32252 ports, isolated from user's worker:
deliberately seeded old model socket in worker process, verified cleanup,
probe PASS, another PASS after11s with no model telemetry, synthetic two-batch
completion and persisted audit. Also existing recorded-validity regression
still passes. Probe is not hardware timing/trajectory validation.

Logs also contain unused old reference des_phi/des_theta at.01s in addition
to actual online reference at.05s. Actual reference is zero. Installer now
disables logging unused reference outputs for future generated copies; do not
confuse those unused traces with the PAC target.
