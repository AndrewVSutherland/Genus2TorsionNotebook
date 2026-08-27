# Boundary-component audit for `[3,18]` and `[6,12]`

## Scope and reproduction

This is a bounded local resolution of the boundary cases left open by
`contact9_318_finite_scout.sage` and
`contact6_m36_612_local_feasibility.sage`.  It does not run another rational
height search and it is not a global resolution of either parameter space.

The reproducible files are

```text
code/task4_boundary_component_audit.sage
data/task4_boundary_component_audit.txt
```

Run:

```bash
timeout 240 sage code/task4_boundary_component_audit.sage \
  > data/task4_boundary_component_audit.txt
```

The completed run takes about four seconds in the current environment.

## Contact 9 and `[3,18]` at 7

Put `t=eps*s`.  Since the family depends on the sign and parameter only
through `eps*s^9` and `s^2`, the two sign charts are the same `t`-chart.
The exact discriminant is factored by the script.  Modulo 7 the finite
special fibres are

```text
t=0,1: (x-1)^2*x^3
t=2:   (x+3)*(x^2+4*x+6)*(x^2+6*x+4)
t=3,5: (x+3)*(x+1)^2*(x^2+2*x+5)
t=4:   (x+1)^2*(x^3+4*x^2+6*x+5).
```

The good `t=2` fibre is the previously excluded group `[2,36]`.  The fibres
at `t=3,4,5` have one split ordinary node.  Their normalizations and Neron
identity groups give

```text
t=3,5: E(F_7)=[12], #J^0(F_7)=12*6=72, 3-part=3^2
t=4:   E(F_7)=[9],  #J^0(F_7)= 9*6=54, 3-part=3^3.
```

The first exceptional directions have node thicknesses `2` or `4` at
`t=3,4`, and `1` or `2` at `t=5`.  Thus their component groups add no
3-primary order.  This proves that `t=3` and `t=5` cannot contain the
required prime-to-7 subgroup `[3,9]`.

At `t=4`, the two points above the node differ by an order-3 class which is
3-divisible in `E(F_7)=[9]`.  Consequently the generalized Jacobian has the
needed rational 3-rank as well as 3-primary order `3^3`.  This chart is not
locally excluded.  In the original variables it is

```text
eps*s = 4 mod 7,
(s,eps) = (4,+1) or (3,-1).
```

### The pole and infinity

At the pole `t=-1`, put `d=t+1`, `x=d*z`, and multiply the equation by
`d^4`.  The normalized special polynomial is

```text
(1/64)*z^4 - 1/4,
```

which is squarefree modulo 7 and has genus-1 normalization of order 8.
The fifth branch point approaches infinity to order 9.  The resulting split
node has thickness 18.  Hence the identity torus contributes one factor of
3 and the component group has 3-part `9`; its total 3-primary capacity is
again `3^3`.  This is a capacity test, not proof that the full group is
`[3,9]`, so the pole is a secondary live chart.

At `t=infinity`, after setting `d=1/t` and scaling by `d^2`, the coefficient
orders are

```text
[1,1,1,1,0,2]
```

and the special polynomial is `x^4`.  The four-root cluster has Newton slope
`-1/4`, so a ramified extension of degree 4 is needed before it separates.
The `t=0`, `t=1`, and infinity charts are not decided by this first ordinary-
node analysis.  They should not be attempted before the simpler live `t=4`
chart.

## Contact 6 and `[6,12]` at 5

The script re-enumerates the full seven-variable intersection and extracts
all 42 points where the five defining equations have Jacobian rank 5.  Their
six boundary signatures occur with counts

```text
2  A0+contact_q_double+degree_drop+half_q0+v3eq1
8  A0+contact_q_double+degree_drop+v3eq1
4  contact_q_double+half_q0+singular_curve
6  contact_q_double+half_q0+singular_curve+v3eq1
6  contact_q_double+singular_curve
16 contact_q_double+singular_curve+v3eq1.
```

Each rank-five point has exactly 25 lifts modulo 25.  All 1,050 lifts are
enumerated and checked in the original five integer equations.  None is
transverse to every boundary initial form.

Only eight residue points even have a promising partial exit:

```text
(a,b,U,v)=(1,0,3,1),  L=1 or 4,
(q,r) in {(1,3),(2,0),(2,2),(3,4)}.
```

For each, 20 of the 25 lifts make the curve discriminant nonzero modulo 25.
They remain on both `U^2-4*v^2=0` and `v^3-1=0` modulo 25, so the marked
cubic contact is still degenerate.  The audit follows all 160 such lifts one
more digit.  Their 4,000 lifts modulo 125 all have exit set exactly

```text
{singular_curve};
```

none leaves either torsion-degeneracy divisor.  The exact substitution
`(a,b,U,v)=(1,0,-2,1)` also makes all three contact equations divisible by
`L^2-1`, explaining the persistent horizontal boundary branch.

## Decision

1. Continue `[3,18]` only on the ordinary-node chart `eps*s=4 mod 7` first.
2. Treat the pole `eps*s=-1 mod 7` as a secondary chart requiring an explicit
   component-group extension calculation.
3. Defer `t=0`, `t=1`, and infinity until those two simpler charts are closed.
4. Stop the current affine `[6,12]` boundary route.  Its only curve-smoothing
   class stays torsion-degenerate through modulus 125; a third blowup has no
   present evidence of reaching the intended open contact locus.

The Neron calculations only test the 3-primary subgroup needed by `[3,18]`.
Passing them is necessary, not sufficient for a rational specialization or
for geometric simplicity.
