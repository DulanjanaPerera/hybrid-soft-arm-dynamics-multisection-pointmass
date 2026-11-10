# Nominal forward block dimension correction

2026-11-09. Operator compilation screenshot reports inconsistent concatenation
at armS_single_dynamics: [dq; M\(...)], followed by unresolved dX dimensions
in residual_FWD_INV_validation/MATLAB Function1.

Read-only saved SLX inspection confirms the forward block's pressureCommand,
X,dX,lengths,tau and pressureUsed symbols all use inherited size-1. The
Integrator initial condition is already[0;0;0;0]. The screenshot does not
identify which inherited operand first became inconsistent. An error at this
concatenation requires resolving the state/acceleration shapes; it does not
indicate a learned-model force magnitude problem.

Paste the entire **rlForwardModelBlock.m** into **MATLAB Function1**, the
nominal forward simulator that outputs dX. This keeps its existing3inputs/
7outputs and original constants/EoM/pressure clipping. It explicitly reshapes
state4x1 and pressure3x1 and constructs torque2x1/derivative4x1. No shared
armS_single_dynamics/core source is changed.

Set these sizes explicitly in Symbols/Edit Data (double, fixed size):

| Symbol | Size |
|---|---|
| pressureCommand_bar | [3 1] |
| X | [4 1] |
| stiff | 1 |
| dX | [4 1] |
| lengths_m (wire comp_lengths) | [3 1] |
| tau_N | [2 1] |
| pressureUsed_bar | [3 1] |
| theta_rad,phi_rad | 1 |
| thetaObservable | 1,boolean |

Also explicitly size desiredPressurePrediction_bar and
nominalPressurePrediction_bar in the desired-pressure block as[3 1]. Keep
desiredMotion[6 1], diagnostic scalar existing bus. If its pressure output
still arrives as a1-D vector at a block requiring a column, a Reshape to[3 1]
can make that interface explicit; do not transpose it to a row or change its
three channel order. Keep the Integrator state column initialization.

The derivative must be [dl2;dl3;ddl2;ddl3]. Each term stacked by the RHS is2x1.
Do not silence this error by setting dX to a matrix or by flattening a wrongly
dimensioned multi-column acceleration after it was computed.

Update Diagram manually after installation, without starting hardware. No
actual user-model compilation was performed during this fix. If dimensions
remain unresolved, inspect the remaining compile diagnostic and upstream port
dimensions rather than guessing a different physics parameter.

## Direct actuation is a separate integration task

The desired-pressure block was installed for observation and returns NaN/
valid=false during cold or invalid forward contexts. A valid prediction alone
does not provide pressure enable/startup, last-command retention, dropout
grace/fault handling or a physical Kill supervisor. Connecting it directly
to valves removes the original startup source while it is still warming up.
Preserve the existing protected command path; do not bypass validity or force
invalid predictions to be treated as valid. An explicit supervised transition
is required before a residual pressure prediction is used as a valve command.
The current change addresses numerical shape only, not that integration.

## Offline evidence

residual_controller/forward_shape_validation_20261109_v1/summary.json:
50states/pressures,100row/column variants, explicit output shapes pass.
Derivative difference from the unchanged original nominal RHS <=8.882e-15.
Includes inputs outside0..3bar to verify unchanged clipping. No SLX edited,
model/hardware/acquisition run, commit or push.
