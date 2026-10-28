# Candidate05 rejected; next manual run is candidate06

Reviewed2026-10-28. Candidate05 matched its exact proposed tuning and logged
software initialization, with no live edits/protocol mismatch. Physical reset is
per standing user attestation. Tracking faults0, pressure limits0; coverage99.17%.

Last3circles110-140s: l2/l3 RMSE=[1.83602,2.06862]mm, vector2.76589mm,
velocity .00372471m/s (passes .01), roughness .03358536m/s^2 (fails .025).
Valid scoring coverage99.333%, pair98.667%. Length limits fail too. No parameter
rate-cap contacts; no stiffness/damping bound contact. Raising caps is unsupported.
Roughness first exceeds .025 on lap7, remains elevated, and reaches .037364 on
lap12 despite declining length error. Near-stop10.30%, phiPeak1.196 are secondary
diagnostics. Reject this candidate; do not repeat it as confirmation.

Candidate03 remains best observed (vector2.21092mm), not accepted/confirmed.

Follow-up on perceived quicker response: first-lap l2 RMSE is lower for05
(3.258mm vs03 3.967mm), but l3 is higher(6.219 vs5.668mm); combined first-lap
error7.021 vs6.918mm. A descriptive sine/cosine/DC least-squares fit at the10s
cycle frequency on valid samples gives early20-50s lag03 [.674,.844]s vs05
[.669,.931]s; late110-140s03 [.185,.226]s vs05 [.292,.321]s. This approximation
does not resolve individual catch-up bursts, but does not support overall faster
phase tracking for05. Crossing desired positions briefly differs from sustained
cycle error or convergence below both1mm thresholds.
Candidate05 vector error is about25.1% higher. Three gammas changed together and
this run occurred the following day, so individual gain causality or unchanged
physical conditions cannot be established solely by this comparison.

Next initial-design point, not a demonstrated improvement:

```matlab
k_g = 0.42192531104075343;
d_g = 0.0020323562614174683;
h_g = 2.5725609289380671;
```

Change only these grouped gamma sources before the same full-reset procedure.
All other gains/bounds/rates/history/switches/pressure settings stay fixed
(activation.2, slew6, lambda18/16.5, Ks15/10). No payload, phi1,12circles10s,
tracking20-140s,ramp140-145s,pressurizedhold145-150s then manual pressurekill.
Save out in the parent folder as candidate06_00g_phi1_T10_12cycles_r01.mat and
send it before another candidate. Preserve interrupted runs with notes.

Five observations counted,five remaining including confirmations. Current
candidate remains exploratory; gammaD is close to original .002 baseline, but
this alone does not guarantee smoothness. No gains applied, hardware/model run,
controller edits, commit or push by assistant. See artifacts/candidate05_r01 and
status_after_candidate05 for complete metrics and provenance.
