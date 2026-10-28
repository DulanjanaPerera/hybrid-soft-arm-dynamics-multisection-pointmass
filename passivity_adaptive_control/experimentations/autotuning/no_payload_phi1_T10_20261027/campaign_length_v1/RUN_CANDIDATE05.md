# Candidate04 reviewed; next manual run is candidate05

Candidate04 matched the archived proposal and fixed settings, with fresh logged
software initial state and no live tuning edits. Physical reset/calibration is
per standing user attestation. Trial imported as candidate04_r01, valid_infeasible.

Late110-140s: l2/l3 RMSE=[1.69217,1.94034]mm, vector2.57456mm,
velocity .00302789m/s, velocity-error-rate roughness .02072852m/s^2.
Coverage98.833%, pair97.667%. Both1mm length limits fail; other limits pass.
Tracking coverage99.125%, adaptation98.25%, pressure limits0 and faults0 during
tracking, no upper bound contact; rare early k2/k3 lower and k3 rate contacts.
Combined length error is about16.4% above candidate03 (2.21092mm), which remains
best observed. This is not an isolated causal test of any one gamma value.

Next initial-design exploration point, not a demonstrated improvement:

```matlab
k_g = 0.28570361215890433;
d_g = 0.0077863906974842481;
h_g = 3.5952162270591348;
```

Change only these grouped gamma sources. All bounds, rates, history, pressure
activation .2/slew6, lambda18/16.5, Ks15/10 and switches remain fixed. Full reset,
no payload, phi1, 12circles at10s, complete150s sequence (20startup,120tracking,
5ramp-down,5pressurizedhold then manualkill). Use exact proposal precision.

Save out in the parent experiment folder as:
candidate05_00g_phi1_T10_12cycles_r01.mat

Send the recording before another candidate. Keep any interrupted recording and
notes under a distinct name. Four observations counted, six remaining including
confirmations; no candidate meets all targets yet. No live gains applied, model
run, controller/SLX edit, commit or push by assistant.
