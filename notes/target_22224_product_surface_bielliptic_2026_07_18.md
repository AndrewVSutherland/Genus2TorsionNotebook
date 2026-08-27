# A rational full-cover surface and why it is not the simple target

Date: 2026-07-18

While ranking the `tor2228.txt` bank by the number of successive good primes
at which `3` divides the Jacobian order, the longest-lived rows all satisfied

```text
a*d = b*c.
```

After projective scaling, write

```text
(a,b,c,d) = (1,r,s,r*s).
```

The four squareclasses for the full `A(2,2,2,8)` cover reduce to

```text
R0 = r^2*s^2,
R1 = (1+r)(1+s)(1+r*s),
R2 = r^2*(1+r)(1+s)(r+s),
R3 = s^2*(1+r)(1+s)(r+s).
```

Put

```text
q = 2*t/(1-t^2),
s = (1-q^2*r)/(q^2-r).
```

Then

```text
(1+r*s)/(r+s) = q^2,
(1+r)(1+s)(r+s)
    = (1+q^2)*((1-r^2)/(q^2-r))^2,
1+q^2 = ((1+t^2)/(1-t^2))^2.
```

Thus all four radicands are rational squares.  This gives a rational
two-parameter surface on the full order-8 cover.

The three strongest near-`3` bank seeds lie on very small fibres:

| tuple | `r` | `t` | first killing good prime |
|---|---:|---:|---:|
| `(754,1716,4611,10494)` | `66/29` | `1/2` | `67` |
| `(65,185,247,703)` | `37/13` | `1/2` | `53` |
| `(14,116,259,2146)` | `58/7` | `2/3` | `61` |

These unusually long local streaks do not indicate a likely geometrically
simple record.  Let `k=r*s`.  For

```text
C: y^2=x(x+1)(x+r^2)(x+s^2)(x+k^2),
```

the branch set is invariant under `x -> k^2/x`, and

```text
iota(x,y) = (k^2/x, k^3*y/x^3)
```

is a rational non-hyperelliptic involution.  Therefore `C` is bielliptic and
its Jacobian is geometrically `(2,2)`-isogenous to a product of elliptic
curves.  No point on this surface can establish the desired geometrically
simple `[2,2,2,24]` record.

The surface remains useful as a source of forced `11`/`13` boundary templates
and as a transverse deformation centre.  Finite-fibre mask code is in
`code/target_22224_product_surface_masks.m`; corrected contact pullbacks for
the fibres `t=1/2,2/3` are families 5 and 6 in
`code/target_22224_a2228_curves_plus3_geometry.m`.
