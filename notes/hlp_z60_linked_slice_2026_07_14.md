# The linked `B5=q0` lane and the smallest transverse pullback

## Outcome

The tempting link at the split HLP cyclic-`[60]` point,

```text
B5 = q0 = x^2-1728/125,
```

does not provide the desired transverse family.  It cuts the marked moduli
to a one-dimensional lane, but that lane is tangent to the unique
Humbert-4 branch.  In a concrete one-dimensional chart its compatible
degree-two involution lifts through order `30` exactly.  Thus this is a
locally split/dead lane for the present search, not a shortcut to a simple
cyclic-`[60]` Jacobian.

The genuinely transverse construction remains the explicit rational base
line obtained by integrating the order-3 identity.  Above it, however, the
order-5 and order-4 data are finite algebraic pullbacks, not rational
functions supplied by the tangent calculation.

## The linked tangent calculation

Use the monic normalization of the three marked identities

```text
A^2-F*B^2 = k5*q5^5,
H^2-F     = k3*q3^3,
(q0*L)^2-F = k4*q0*u4^2.
```

At the seed, `B=q0`.  The two coefficient equations `B=q0` have independent
differentials, so they have the expected codimension two in the
three-dimensional marked moduli space.  Nevertheless the exact differential
ranks are

```text
rank(link equations)                    = 2,
rank(link equations + Humbert-4 normal) = 2.
```

Thus the Humbert normal is already in the span of the two link equations.
The linked locus has no first-order direction leaving the split surface.

For a completely concrete affine curve, additionally impose

```text
q0 = x^2-1728/125,
L(0) = -2775,
k4 = 46250000.
```

Eliminating `F` through the order-4 identity leaves `19` variables and `18`
equations.  Their Jacobian has rank `18` at the HLP point.  Hence this is a
smooth one-dimensional local marked branch, but its unique tangent pairs to
zero with

```text
N_H4 = (0,81125,0,904800,0,9916416,0).
```

## Exact formal split test

The script constructs the unique marked formal branch one coefficient at a
time.  Every step is a rational solve with a matrix no larger than `19x19`.
It then seeks an involution of the binary sextic in the form

```text
(X:Z) |-> (-X+b(s)Z : c(s)X+Z)
```

and solves

```text
F_s(-X+bZ,cX+Z) = lambda(s)*F_s(X,Z).
```

The seven invariance equations lift through order `30`.  Moreover, at every
computed order,

```text
b(s) = (-1728/125)*c(s),
lambda(s) = (1-(1728/125)*c(s)^2)^3.
```

These are precisely the multiplier relations for an involution preserving
the quadratic form `X^2-(1728/125)Z^2`.  Independently, the exact germ
calculation in `notes/hlp_z60_linked_germ_2026_07_14.md` proves more: the
smooth five-dimensional linked germ is the coordinate saturation of a
three-dimensional even-sextic subgerm.  It is therefore locally contained in
Humbert 4 to all formal orders.  The order-30 jet is a direct check of that
all-orders dimension/submersion argument, not the sole basis for it.

A Padé scout using the same `30`-jet found no rational expression of
numerator and denominator degree at most `4` for any nontrivial marked
coordinate (with four unused coefficients held back for validation).  So
there is also no small rational parametrization hiding in this chart.

## Why the transverse line is only a finite pullback

Integrating the nonzero order-3 tangent coordinates gives the exact rational
base family

```text
q3(t) = x^2-316/25 -(3/1620896000)*t*x,
H3(t) = 29600-2775*x^2
        -(137/202612)*t*x +(75/1620896)*t*x^3,
F(t)  = H3(t)^2-46250000*q3(t)^3.
```

It satisfies `F(0)=F_HLP` and `dF/dt|_0=x`, so it really is transverse to
Humbert-4.  The order-3 layer is rational by construction.  The remaining
normalized blocks have the following square sizes:

```text
order 5: 11 variables / 11 equations,
order 4:  7 variables /  7 equations.
```

Both auxiliary Jacobians are invertible at `t=0`.  Therefore their pullback
is a one-dimensional algebraic curve, etale over the `t`-line at the HLP
point, with a unique formal point over `Q[[t]]`.  This proves existence of a
formal lift, not rational functions in `t`.

The order-4 block makes the obstruction especially visible:

```text
F = q0*(q0*L^2-k4*u4^2).
```

Hence a rational specialization first needs a rational quadratic factor
`q0` of `F(t)`.  Quadratic factors of a generic sextic already form a finite
degree-`15` cover (the `15` unordered pairs of six roots), before imposing
the halving data.  The order-5 norm identity supplies a second finite cover.
A rational `t != 0` must make rational points appear in both covers at once.

Consequently the smallest honest transverse object now available is the
explicit `t`-line together with `18` finite auxiliary equations.  It is not
an explicit rational cyclic-`[60]` family, and treating the implicit-function
lift as one would be incorrect.

## Reproduction and memory

```text
python3 code/hlp_z60_linked_slice_tangent.py
python3 code/hlp_z60_linked_slice_formal.py
python3 code/hlp_z60_linked_slice_reconstruct.py
```

Measured peak resident memory was approximately `11.5 MB`, `12.0 MB`, and
`55.0 MB`, respectively.  No Groebner basis or large elimination was run.
