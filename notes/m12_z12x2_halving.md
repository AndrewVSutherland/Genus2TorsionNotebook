# Halving the extra 2-torsion in the `M(12)` family

This tests whether the one-parameter family

```text
a = (1-r)/4
```

with known subgroup `Z/12Z x Z/2Z` can be promoted to `Z/12Z x Z/4Z` by halving the independent rational 2-torsion class.

## Setup

Move the rational Weierstrass point `x=2` to infinity.  The two finite rational roots in the resulting odd quintic are

```text
beta_div = (1-r)/(4r),
beta_ind = (2-r)/(4(r-1)).
```

Here `beta_div` is the class `6D`, where `D` is the order-12 class, so it is already divisible by `2`.  The class at `beta_ind` is the new independent rational 2-torsion class.

For `f(X)=(X-beta)g(X)`, a half of `[beta-infinity]` is equivalent to

```text
g(X) - (X-beta)(mX+n)^2 = c (X^2 + A X + B)^2,
```

where `c` is the leading coefficient of `g`.  This is the corrected tangent-polynomial identity; the tangent polynomial is `H(X)=(X-beta)(mX+n)`.

## Exact search

The script

```text
code/m12_z12x2_halving_search.m
```

uses Magma's exact `IsDivisibleBy` on integral models isomorphic to the odd quintics.  A run to height `40` gave

```text
checked 1955
good 1955
known_ok 1955
hits 0
```

The `known_ok` count is the sanity check that the `beta_div=6D` class is divisible by `2` in every tested case.

## Finite-field sieve

The script

```text
code/m12_z12x2_halving_finite_field_sieve.m
```

finds a good-reduction obstruction modulo `7`.

Modulo `7`, the good residues are

```text
r = 3, 4.
```

For both residues the known class `6D` is divisible by `2`, but the independent 2-torsion class is not.  The other residues are boundary or bad:

```text
r = 0,1,2  boundary/denominator
r = 5,6    discriminant zero
```

Thus there are no nonboundary good mod-`7` points on this one-parameter family with `Z/12Z x Z/4Z`.  Any rational example on this line would have to reduce to the boundary or bad locus modulo `7`.

## Conclusion

The direct halving route on the clean one-parameter `Z/12Z x Z/2Z` family does not currently produce a strictly simple `Z/12Z x Z/4Z` example.  Away from the mod-`7` boundary it is obstructed.

The remaining realistic continuation is a local boundary analysis at `p=7`, or else leaving this one-parameter line and imposing the halving equations on the full two-dimensional extra-Weierstrass surface inside `M(12)`.

## Boundary follow-up

The script

```text
code/m12_z12x2_halving_boundary_sieve.m
```

does the first boundary pass.

For residues where the odd-quintic tangent equations still make sense modulo `7`, the independent halving equations have the following mod-`7` behavior:

```text
r = 2,3,4,5  no tangent solutions
r = 6        two tangent solutions
```

The residues `r=0,1` are coordinate boundary values for this odd model.  The branch at `r=6` does lift modulo higher powers of `7`, so the mod-`7` boundary is not killed by a first-order Hensel check.

A combined residue sieve at

```text
7, 11, 19, 23, 31, 43, 47, 59, 67, 71
```

is much stronger.  Here a residue is allowed if it is bad/boundary or if the independent 2-torsion class is divisible over the finite field.  Up to height `300`, the only rational parameter surviving all these local residue constraints is

```text
r = -1.
```

But at `r=-1` the odd quintic is singular:

```text
f5 = (X + 3/8)(X + 1/2)^2(X^2 - 2X - 1).
```

So the boundary pass still finds no nonsingular `Z/12Z x Z/4Z` specialization on this one-parameter line.
