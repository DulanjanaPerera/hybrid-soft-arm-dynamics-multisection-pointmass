# Circle test4: the nominal inverse targets differ

2026-11-09. Offline inspection/replay of circlesweep_FWD_INV_test_4.mat only.
Model1.54, operator StopCommand at271.6s,5433 observer rows,5118 valid inverse
contexts. The MAT contains only out; it does not retain the option variables.
The current saved SLX was inspected as ZIP/XML without loading Simulink or
executing callbacks. It corroborates the desired q/v/a collection input and
the acceleration-only single-pass reference input. It is supporting evidence;
full numeric replay of the actual recording establishes executed behavior.

## Finding

Define F0(q,v,a)=M0(q)a+(C0(q,v)+D0)v+K0q+bound(q)+G0(q).
Let A be the minimum-nonnegative pressure allocation (max3bar), followed by
the shadow slew rule where applicable. The traces are:

- Actual pressureRequested: collection A(F0(qDesired,vDesired,aDesired)).
- Single-pass nominal: A(F0(qMeasured,vMeasured,aDesired)).
- Single-pass corrected: A(F0(qMeasured,vMeasured,aDesired)-rActual).
- rActual: one GRU call at actual q/v, actual NI command, physical dt and
  immutable actual pre-observation history. No desired q/v enters that call.

Thus the two nominal traces are not theoretically the same calculation.
The new shadow block has no desired-length or desired-velocity input. It is
an observed-state force/pressure reconstruction, not a desired-trajectory
inverse feedforward block. This was documented but was not the comparison
the user intended. Matching pressure conversion did not match force targets.

Recomputing F0 at desired q/v/a matches logged tauDesired within3.56e-15N.
Applying the collection map matches pressureRequested within4.45e-16bar.
Recomputing F0 at measured q/v plus desired a matches shadow target exactly.
Replaying current source matches both pressure outputs within1.12e-15bar.
Actual NI commands, pressureRequested, and inverseContext pressure agree
exactly on valid rows; commandAge=0. This establishes a target mismatch,
not a subtraction sign error or pressure-output wiring interchange.

## First circles,13.00-62.95 simulation seconds

969 valid contexts. Desired absolute peak lengths:[2.5999,2.5999]mm.
Measured absolute peaks:[0.3742,0.1796]mm; includes the observed static offset.
Length tracking RMSE:[1.7772,1.8064]mm. The small actual response is real in
the logged orientation-derived lengths; correct forward prediction is not
evidence that the desired trajectory was achieved.

| Peak pressure by channel (bar) | P1 | P2 | P3 |
|---|---:|---:|---:|
| Actual command / desired-state nominal inverse | .2423 | .2403 | .2423 |
| Current measured-state nominal shadow | .0000 | .0430 | .0325 |
| Current measured-state corrected shadow | .2584 | .2585 | .2522 |
| Desired target minus same frozen actual residual, offline only | .5228 | .4682 | .4835 |

The last row is a diagnostic counterfactual before slew, computed as
A(F0(qDesired,vDesired,aDesired)-rActual). It illustrates the missing force
target, not validated required pressure or a guaranteed tracking result.
Peaks occur at different times. Generalized-force corrections need not raise
every channel at every instant because allocation has a common-mode freedom.

Corrected vs actual pressure RMSE in these circles:[.00660,.00873,.00694]bar.
It is almost reconstructing the command that produced the observed response.
Corrected slew limiting:0/5118 valid contexts across the entire recording;
corrected max-pressure scaling:0. The closeness is therefore not caused by
the .3bar-per-tick clipping in this experiment. Nominal shadow slew limiting
occurs3026/5118 later in the recording, but0 times in the first70s.

## What forward accuracy does and does not establish

The forward GRU model takes actual current state/command and predicts the
next50ms response. Each prediction origin is reinitialized from measurements;
this is not an inverse-pressure predictor or autonomous desired tracking.
Full-record matched prediction length RMSE (mm):

| Predictor | l2 | l3 |
|---|---:|---:|
| Nominal |3.7507|3.8072|
| Residual + nominal |.2257|.2315|
| Constant velocity |.0547|.0581|

During the first70s, augmented RMSE is only[.0658,.0552]mm while CV gives
[.0050,.0064]mm. The arm barely moves there, making short-horizon prediction
easy. These results support better local prediction than nominal, not an
identified deadzone, accurate pressure sensitivity, or a certified inverse.

Because r depends on pressure, freezing rActual at a different hypothetical
pressure is approximate. Reevaluating the network at the offline desired
candidate pressure with actual q/v/history gives full-record model-balance
error norm median/P95/max[1.4778,2.2331,3.5686]N. This is a model consistency
diagnostic, not measured physical force error. No force or chamber-pressure
sensor data establish the physically required pressure here.

## Next implementation decision

For an apples-to-apples desired inverse comparison, supply a shared nominal
force demand, ideally the existing logged tauDesired_N2x1, to a separate
desired-demand subtraction block. Its nominal output maps that exact demand;
its corrected output maps demand-rActual. Keep the GRU's actual measured
state, actual command and causal history. Do not substitute desired states
into its recurrent observation stream. Apply/log the same pressure policy
and constraints to both branches and retain pre-limit outputs.

For PAC, the demand must instead be the complete PAC force demand including
reference shaping, feedback and adaptive correction. A collection feedforward
target is not a substitute for that controller. Candidate-pressure residual
consistency and physical tracking still need separate validation.

No runtime/controller/NN source or SLX change was made in this audit. No
model, hardware, acquisition, training, commit or push was run.

## Retained evidence

- target_comparison_first70.png/.pdf/.fig: desired/measured lengths, actual
  vs observed-state nominal/corrected pressure, offline desired-target case.
- summary.json,first_window_metrics.csv,trace.csv,analysis.mat.
- saved_model_readonly_wiring.json and input_sources_sha256.json.
- independent_csv_checks.json:5118 rows; frozen balance1.57e-13N,
  desired nominal balance1.29e-7N (within intentional1e-8bar zero threshold),
  no corrected slew difference and exact actual/requested command identity.
- Reproduce in a fresh directory with analyzeSinglePassTargetMismatch(root,folder).
