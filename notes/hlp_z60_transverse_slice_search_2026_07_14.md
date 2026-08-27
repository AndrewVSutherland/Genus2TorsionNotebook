# Arithmetic search on a transverse line through the split `[60]` seed

## Outcome

The exact marked tangent calculation gives a genuine direction away from the
`(2,2)`-split Humbert surface, but the sharpest small coefficient line through
that direction contains no second rational order-`60` specialization of
projective height at most

```text
81,505,794.
```

This is a rigorous result on one transverse line, not a global nonexistence
statement for geometrically simple cyclic `[60]` Jacobians.

## The selected line

Start with the exact HLP control

```text
F0 = -46250000*x^6 + 1761500625*x^4
     -22332312000*x^2 + 94277468160
```

and take

```text
F_t = F0 + t*x.
```

The deterministic marked calculation in
`code/hlp_z60_marked_tangent.py` proves that `dF=x` has a unique rational
first-order lift through the order-5, order-3, and order-4 identities.  The
normal to the unique Humbert-4 branch at the seed is

```text
N = (0,81125,0,904800,0,9916416,0),
```

so `N(x)=81125 != 0`.  Thus this is not a coordinate direction and the line
is not locally contained in the split surface.

As an independent geometric control, the specialization `t=1` has an
absolutely simple reduction at `p=11`.  Its Frobenius polynomial is

```text
T^4 - 3*T^3 + 7*T^2 - 33*T + 121,
```

and the minimal polynomials of the first twelve powers of a Frobenius root
all have degree four.  Hence the Jacobian at `t=1` is geometrically simple
over `Q`.  This sample does not have rational order `60`; it only proves that
the selected line really reaches the simple locus.  Reproduce this check with
`code/hlp_z60_transverse_sample_simple.m`.

## Rigorous finite masks

Let `t=n/d` with `gcd(n,d)=1`.  At a prime `p` not dividing `60`, a rational
point of order `60` injects into the Jacobian of every good reduction.  In
particular the exponent of `J_t(F_p)` must be divisible by `60`.

The bounded script

```text
code/z60_split_seed_finite_masks.m
```

computes the complete affine `t mod p` mask by constructing the finite
Jacobian and checking its invariant factors.  Singular fibers and any
computation failure are retained in a separate `bad` set, so they are never
used to exclude a rational parameter.  The complete output through `p=181`
is

```text
data/z60_split_seed_tx_masks_p181.txt.
```

Ten primes have the especially strong result

```text
allowed affine t = {0},     bad = empty:

7, 11, 13, 19, 29, 73, 79, 107, 131, 149.
```

Their product is

```text
M = 6643194523826861.
```

The original seed passes every mask; for example the finite invariant factors
at `p=7,11,13` are respectively `[120]`, `[120]`, and `[4,60]`.

## Height certificate without enumeration

For each singleton prime `p` above:

- if `p` does not divide `d`, the affine reduction is defined and the mask
  forces `n/d = 0 mod p`, hence `p | n`;
- if `p | d`, the parameter reduces to the projective point at infinity and
  is conservatively retained.

Since `gcd(n,d)=1`, every singleton prime divides exactly one of `n,d`.
For nonzero `n` this gives

```text
M | n*d.
```

But

```text
floor(sqrt(M)) = 81505794
```

and `81505794^2 < M`.  Therefore no nonzero primitive `n/d` with

```text
max(|n|,d) <= 81505794
```

can pass these necessary order-`60` conditions.  The only specialization in
that box not excluded by these necessary conditions is `t=0`, whose exact torsion is `[60]` but whose Jacobian is split.

No numerator/denominator enumeration was required.  The full mask run took
about five seconds and about 32 MB in the recorded run, with a hard 2 GB
Magma memory limit.

## Interpretation

The deformation succeeds geometrically and fails arithmetically on this
particular line:

- the three marked primary layers have no infinitesimal obstruction;
- the line leaves the split locus and contains geometrically simple curves;
- nevertheless its finite marked pullback has no second rational order-`60`
  point in a very large height range.

The tangent lift is a rational formal branch, not a rational-function
parametrization of the order-5/order-4 pullbacks.  The most obvious exact
relation, `B5=q0`, was also tested: its smooth local germ is the coordinate
saturation of an even-sextic subgerm and is therefore contained in Humbert 4
to all formal orders.  See `notes/hlp_z60_linked_germ_2026_07_14.md`.
Consequently a useful next move must let `B5` and `q0` separate and analyze a
different marked auxiliary curve; merely extending the coefficient height on
this line is not competitive.

For comparison, `code/z60_split_seed_direction_scout.m` ranks the nearby
lines `F0+t*(a+x)` for small integer `a`.  They are all locally thin at small
primes; `a=0` gave the cleanest rigorous certificate above.
