# V2 saved-model wiring review

2026-11-09. Inspected residual_FWD_INV_validation_v2.slx read-only as ZIP/XML;
model SHA retained in review.json. No model load/callback/compile/edit/run.

Correct saved connections:

- MATLAB Function8 corrected pressure -> Switch upper/true data input1.
- Original MATLAB Function3 pressureRequested -> Switch lower/false input3.
- MATLAB Function9 useResidual -> Switch control input2.
- Switch -> Unit Delay -> Pneumatic Valves; delay Ts0.05s,IC[0;0;0].
- Existing downstream Kill/saturation paths retained.
- Command feedback is assembled after Kill/saturation, sampled with its clock
  at0.01s, coherently rate transitioned to0.05s, split[3 1] and routedA8/A9
  to observer pressure and timestamp inputs6/7.

Required correction:

MATLAB Function9 input2(valid) currently branches From6/A5. This resolves to
NDI Sensor reading output7 baseOreintationValid, also observer input3. It is
not the desired-pressure block's readiness flag. The validity function source
is otherwise correct: useResidual=valid&&all(isfinite(p)).

Branch MATLAB Function8 output3(desiredPressureDiagnostic bus), add a Bus
Selector selecting **valid**, and connect that scalar Boolean to MATLAB
Function9 input2 in place of A5. The same bus can remain logged/terminated
on another branch. Keep Function9 input1(p) connected to corrected pressure.
Use Switch criterionu2>0 with threshold0, or u2~=0, for the Boolean selector.

The existing finite test already rejects cold/invalid NaN predictions under
the present desired-pressure source; the wrong base flag should nevertheless
be replaced by the actual readiness flag, not treated as equivalent readiness.

The command-path delay matches the recommended break of the reported
same-timestep loop and introduces50ms command delay. This static review does
not establish successful full-model compilation, timing or closed-loop
tracking/startup/dropout safety. Physical Kill/enable/stop retain priority.
No further block-code change is needed for this wiring correction.
