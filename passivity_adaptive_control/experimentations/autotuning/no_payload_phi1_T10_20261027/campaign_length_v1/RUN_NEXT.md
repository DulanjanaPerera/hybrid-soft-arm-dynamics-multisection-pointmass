# Next operator-run trial: candidate03

Update: candidates03 through08 have been reviewed and recorded. See
REVIEW_AFTER_CANDIDATE08.md before another run; candidate09 is not queued.
The original candidate03 instructions below are retained as history only.

Prepared from two existing baseline screenings and user-agreed limits/bounds.
Campaign ready for manual intake; no gain was applied and no hardware/model run.
Review candidates/candidate03/review.md for the full exact proposed configuration.

Change only the three grouped gamma sources:

```matlab
k_g = 0.645028694505006;
d_g = 0.0010090123623767366;
h_g = 5.8362570107825071;
```

These values are an initial-design optimizer proposal, not demonstrated improved
gains. Preserve their precision so logged tuning matches the archived candidate.
Keep lambda=[18;16.5], Ks=[15;10], bias gamma=1, all bounds/rates/scales, history,
pressure settings (activation .2, slew6), supervisor and switches unchanged.

Fully depressurize/reset, freshly restart and repeat the same calibration/trim
and preconditioning procedure per the user's attestation. No payload; phi1;
10 s/circle, 12 circles. Tracking20-140s; ramp-down140-145s; pressurized hold145-150s.
Keep the existing manual pressure kill/stop method after the hold. Do not stitch
trials together just by toggling pressure: persistent controller/history/supervisor
state and the physical reset would need a separate verified reset implementation.

After completion save out to the parent experiment folder as:
`candidate03_00g_phi1_T10_12cycles_r01.mat`. Do not overwrite an existing file.
Retain aborted/faulted runs with an explicit suffix and notes; do not resume them
as the same fresh run. Send the saved path for offline import before another trial.

Acceptance over110-140s: q2/q3 RMSE<=1mm each; vector length velocity RMSE<=.01m/s;
velocity-error-rate roughness<=.025m/s^2; pressure limiting=0; minimum valid
coverage95%, no latched fault, fixed logged tuning. Angular limits disabled.

Budget interpretation is explicitly conservative:10 total recorded observations,
including2 existing baselines and any confirmations, leaving8 additional runs.
Two fresh confirmation repeats are required after a candidate meets the targets;
if budget runs out before that, stop and review rather than add runs automatically.
The current model/source snapshot is archived, but is not an embedded historical
model artifact from the baseline files; physical reset/calibration is user-attested.
