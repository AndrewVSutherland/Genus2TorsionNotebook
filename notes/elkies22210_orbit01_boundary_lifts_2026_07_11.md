# Orbit `01` boundary lifts at 3 and 7

The affine unit-circle chart for the Elkies Clebsch--Klein Orbit-`01`
halving cover has no points modulo `3` or `7`.  This does **not** give a
local obstruction: the proper `(P1)^4` closure has boundary points, and
there are genuine smooth-open `Q_3` and `Q_7` branches above them.

The exact lift classifier is

```text
code/elkies22210_orbit01_boundary_lifts.py
```

## Projective equations

Write `t_i=X_i/Y_i`, for the four parameters indexed by `i=2,...,5`.
If `E_j` denotes the elementary symmetric function of the four `t_i^2`,
the affine equations are

```text
E4-E1-5=0,
E3-E2+10=0.
```

The script clears `product Y_i^2` and works in the unique integral chart of
each `P1`: either `Y_i=1`, or `X_i=1` with `Y_i=0 mod p`.  It enumerates all
base points and lifts the two exact multihomogeneous equations through
successive powers of `p` by their four-variable affine tangent system.

The counts reproduce the independent closure enumeration:

```text
p=3:  56 mod 3  ->  1512 mod 9  ->  13608 mod 27
p=7: 248 mod 7  -> 28616 mod 49 -> 1402184 mod 343.
```

All mod-`p` points are boundary points.  At `p=3` they have Jacobian rank
`1`; at `p=7`, 56 have rank `1` and 192 have rank `2`.

## Quantitative Hensel certificates

For a fixed pair of variables, let `d` be the valuation of the corresponding
`2 x 2` Jacobian determinant.  The quantitative multivariable Hensel lemma
says that

```text
min_i v_p(F_i(a)) > 2*d
```

produces an exact zero whose displacement from `a` has valuation at least
`min_i v_p(F_i(a))-d`.

The open factors checked by the script are

```text
X_i, Y_i, Y_i^2-X_i^2,
Y_i^2*X_j^2-X_i^2*Y_j^2,
Y_i^2*Y_j^2-X_i^2*X_j^2.
```

Their nonvanishing is exactly `r_i != 0`, `r_i^2 != 1`, and pairwise
distinct `r_i^2`.  Since `-1` is nonsquare at both primes, every
`X_i^2+Y_i^2` is automatically a unit in a primitive projective chart.

For `p=3`, take the chart with the first parameter at infinity and chart
variables

```text
(Y_2,X_3,X_4,X_5) = (9,7,10,22).
```

Then

```text
F(a) = (2319921,185938497),
v_3(F_1)=v_3(F_2)=5.
```

Using columns `0,1`, the Jacobian determinant is

```text
-28562661461376,    v_3(det)=2.
```

Thus the exact root is congruent to this approximation modulo `3^3`.
Every open factor at the approximation has valuation at most `2`, so all
remain nonzero at the exact root.  This proves

```text
H01(Q_3) has a smooth-open point.
```

For `p=7`, use the same projective chart and

```text
(Y_2,X_3,X_4,X_5) = (7,246,212,3).
```

Here

```text
F(a) = (24473311709,1068847960173),
v_7(F_1)=v_7(F_2)=3.
```

The determinant on columns `1,2` is the `7`-adic unit

```text
-49896096399360.
```

The exact root is congruent modulo `7^3`, while every open factor has
valuation at most `2`.  Hence

```text
H01(Q_7) has a smooth-open point.
```

## Conclusion

The absence of affine residue points at `3` and `7` is entirely a chart
artifact.  Both primes have certified open p-adic branches above infinity
boundary classes, so neither prime can obstruct a rational Orbit-`01`
halving point.
