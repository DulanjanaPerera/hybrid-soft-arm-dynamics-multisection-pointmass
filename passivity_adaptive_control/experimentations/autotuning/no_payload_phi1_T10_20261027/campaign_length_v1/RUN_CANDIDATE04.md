# Candidate03 reviewed; next manual run is candidate04

Candidate03 exact tuning/software initialization matched the proposal, with no
live tuning changes, no faults or pressure limiting. 12 circles at10s, phi1.
Physical reset/calibration follows the standing user attestation, not independent
measurement. Original recording remains unchanged; manifest/report are in
artifacts/candidate03_r01. Late110-140s length RMSE=[1.42549,1.69001]mm,
vector2.21092mm; velocity .00314726m/s; error-rate roughness .0207654m/s^2;
coverage98.83%, valid-pair97.67%. Only the two length limits fail. Combined length
RMSE is ~9% below the mean of the two baseline scores, but this is one candidate
run and not proof of repeatable improvement or identified physical parameters.
Lap12 roughness .022242 remains below .025 but is higher than early laps.

Next initial-design point (not a demonstrated improvement):

```matlab
k_g = 0.58034718095284454;
d_g = 0.0037139224909066579;
h_g = 4.4053683649622428;
```

Change only those sources before a fresh full-reset run; all other tuning,
switches, pressure activation .2/slew6, lambda18/16.5 and Ks15/10 stay fixed.
No payload, phi1, 12 circles at10s. Tracking20-140s, ramp140-145s,
pressurized hold145-150s then manual pressure kill using the existing procedure.
Use full proposal precision so logs match candidates/candidate04/candidate.mat.

Save out in the parent experiment folder as:
candidate04_00g_phi1_T10_12cycles_r01.mat

Send the recording before another candidate. Preserve any interrupted run under
a distinct name with notes. Current status:3 observations counted,7 remaining
including any confirmations; no candidate has met every target. Initial design
is still exploring gamma settings; no claim of converged optimizer/arm behavior.
No live gains, controller/model changes, hardware runs, commit or push by assistant.
