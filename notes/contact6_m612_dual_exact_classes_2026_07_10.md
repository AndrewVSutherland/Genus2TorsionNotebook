# Exact dual-class tests for the contact-6 `[6,12]` route

Date: 2026-07-10.

The existing discriminant covers

```text
R1 halved => DB*DC is a square,
R2 halved => DC is a square,
R3 halved => DB is a square
```

are necessary conditions only.  The new file

```text
code/contact6_m612_dual_class_exact.m
```

constructs both standard Richelot twist signs of

```text
Delta*R1*R2*R3
```

and directly calls `IsDivisibleBy(T,2)` on the three rational dual-kernel
classes.  It also computes the exact rational torsion and marks which twist
retains two independent rational `3`-directions.  This gives the missing
exact test for the mixed class `R1`, as well as exact confirmation after the
cheap `DB` and `DC` tests for `R2` and `R3`.

On the known geometrically simple contact-6 seed

```text
(a,b)=(133/39,-7/13),
```

the `+` twist has exact torsion `[6,6]`, retains the rational `[3,3]`, and
none of `R1,R2,R3` is divisible by `2`.  The `-` twist has torsion `[2,2]`
and does not retain `[3,3]`.  Thus the exact class test independently
recovers the previously known negative Richelot result and identifies the
correct twist without comparing models returned by
`RichelotIsogenousSurfaces`.

Reproduce with

```text
magma -b code/contact6_m612_dual_class_exact.m
```

For a search driver, load with `NoMain:=true` and
`run_known_seed:=false`, then call `ExactDualClassRecord(a,b,eps)` for
`eps=+1,-1`.  A genuine target on the distinguished dual must occur on the
twist with `has_33=true`, have at least one exact half, and finish with exact
torsion invariants `[6,12]`.

The companion derivation

```text
code/contact6_m612_dual_halving_equations.m
```

computes the two universal square-quartic equations for each dual class.  No
common factor or low-degree factor appears over `Q` on this raw chart.  For
`R1` the two irreducible equations have shapes

```text
(total degree, M-degree, N-degree, terms)
(12,6,3,225), (16,8,4,700).
```

For `R2` they have shapes `(14,5,3,269)` and `(19,8,4,701)`; for `R3`,
`(13,5,3,248)` and `(18,8,4,648)`.  This makes the intended workflow clear:
use the resultant-square covers and the small-prime boundary conditions as
cheap filters, then call the exact class test only on rational core points.
The universal direct cover is not the right first elimination problem.
