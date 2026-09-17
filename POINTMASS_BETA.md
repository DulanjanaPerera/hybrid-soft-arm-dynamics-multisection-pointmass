# Point-mass translational beta interface

The active three-section point-mass model uses a `3 x 3` numeric `beta`:
row `n` belongs to section `n` (base to tip), and its columns are
`[beta_v1, beta_v2, beta_v3]`. The default is `ones(3,3)`. Beta is held
constant with respect to the six generalized coordinates during a run.
The mass position remains the selected backbone point `cog_xi(n)`; beta
does not change that point or the gravitational potential.

With upstream translational Jacobian `A`, upstream angular blocks applied
to the selected local position `B`, and local Jacobian `P = p_q`, the
section's translational contribution is

```text
M11 = m*(A'*A + A'*B + B'*A + beta_v1*B'*B)
M12 = m*(A'*P + beta_v2*B'*P)
M21 = M12'
M22 = m*beta_v3*(P'*P)
```

For section 1, `A` and `B` are empty and only `M22` contributes. `Mi_h`
applies the complete product rule to these same blocks; the section mass
multiplies both `M` and every derivative slice. `christoffelSymbol` uses
`dM(i,j,h) = partial M(i,j)/partial q(h)`.

The compiled runtime interface is now:

```matlab
beta = ones(3,3);
dX = armS_dynamics_N3_entry_mex_mex(t,X,L,r,cog_xi,mi,g, ...
    Kmin,D,tau,mu,lKbounds,beta);
```

Change `beta` between calls without rebuilding. The interpreted
`armS_dynamics_nume` reads `params.beta`, defaulting to `ones(N,3)` when
that field is absent. The interpreted entry and core also default to ones
if beta is omitted. The MEX requires the runtime beta argument. The two
`runDynamicSimulation_armS_nume*` scripts initialize and pass all ones.

These are the *translational* coefficients of the CoG paper's equations
(22) and (25). The paper's fitted values were not applied. Positive beta
entries alone do not guarantee that the coupled mass matrix is positive
definite; check Cholesky or eigenvalues over the intended operating range.

`armS_dynamics_nume_N3`, `armS_dynamics_nume_v2`, and
`armS_dynamics_recursive` are older alternate implementations and are not
called by the active simulation scripts or MEX entry. They are not a
tunable-beta interface.
