# Z/60 literature, database, and Nicholls order-30 seed audit

Audit date: 2026-07-10.

## Bottom line

No geometrically simple genus-2 Jacobian over `Q` with a rational cyclic
subgroup of order `60` was found in the literature or public databases audited
here.  This is a scoped search result, not a nonexistence theorem.

The published `Z/60` construction of Howe--Leprévost--Poonen is explicitly
split: every Jacobian in it is `(2,2)`-isogenous to a product of elliptic
curves.  The simple order-30 curve in Nicholls, Example 3.1.8, is also a dead
end for doubling the order: its unique rational 2-class is not divisible by 2
already over `F_5`.  The workspace simultaneous contact-5/contact-6 family is
a different family; its complete `Q_11` obstruction does not apply to the
Nicholls curve.

The best construction priorities remain:

1. impose rational `3`-torsion on a simple order-20 surface;
2. impose rational `5`-torsion on a full order-12 surface, keeping an order-4
   class present from the outset;
3. as a secondary route, impose an order-4 class on an order-15 family.

These routes avoid trying to halve the sole order-2 class of an order-30
specialization, which is locally fragile and is now ruled out for both the
workspace contact-(5,6) family and the Nicholls seed.

## What the published `Z/60` result actually gives

Howe--Leprévost--Poonen (HLP) prove that a positive-rank elliptic curve
parametrizes genus-2 Jacobians **containing** a rational subgroup isomorphic to
`Z/60Z`.  Their equation is

```text
(2*t^2-2*t+1)*(6*t^2-6*t+1)*y^2
    = (2*u-1)*(4*u^2-2*u-1).                    (HLP, equation (3))
```

Putting `t=1/3` gives Cremona curve `900A1`, of rank `1`; all but finitely many
of its rational points give the claimed Jacobians.  The construction glues an
elliptic curve with a rational 12-torsion point to one with a rational
10-torsion point along their 2-torsion.

This does not solve the present problem:

- HLP's theorem is subgroup containment, not a claim that the full rational
  torsion subgroup is exactly `Z/60Z`.
- Their genus-2 construction is built so that the Jacobian is `(2,2)`-isogenous
  to the product of the two elliptic curves.  It is therefore geometrically
  split.

Primary source:

- Everett W. Howe, Franck Leprévost, Bjorn Poonen, *Large torsion subgroups of
  split Jacobians of curves of genus two or three*, Forum Math. 12 (2000),
  315--364: <https://math.mit.edu/~poonen/papers/large.pdf>

HLP Table 1 lists `Z/60Z` and Section 3.5 gives equation (3), the specialization
`t=1/3`, and the `900A1` rank calculation.

## State of the record in later sources

Nicholls' 2018 thesis separates known genus-2 torsion orders into geometrically
simple, geometrically split, and unverified columns.  Its Table 3.1 places order
`60` only in the geometrically split column, citing HLP.  Example 3.1.8 gives a
geometrically simple order-30 seed, audited below.

Primary source:

- Chris Nicholls, *Descent methods and torsion on Jacobians of higher genus*,
  Oxford DPhil thesis (2018):
  <https://ora.ox.ac.uk/objects/uuid%3A04cef70a-2ab9-44c2-8bbe-ca2ac33bfe41/files/m231349cc1337ef2c59646e9d0d93cc29>

Elkies' public page also clearly distinguishes the split HLP constructions from
absolutely simple examples.  It records Leprévost's unpublished one-parameter
order-30 family, but no simple order-60 example:

- <https://people.math.harvard.edu/~elkies/g2_tors.html>

Two modern endomorphism-type results sharpen the strategy, but must be kept in
their exact scopes:

- Laga--Schembri--Shnidman--Voight prove that if an abelian surface over `Q`
  has geometric endomorphism ring a maximal order in a non-split quaternion
  algebra, then its rational torsion is killed by `12` and has order at most
  `18`.  Thus this maximal-order potential-QM locus cannot contain the desired
  `Z/60`.
  <https://arxiv.org/abs/2308.15193>
- Alessandrì--Coppola, Conjecture 4.5, predicts the possible torsion orders for
  simple non-CM `GL_2`-type abelian surfaces to be

  ```text
  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
  28,31,37,44,56.
  ```

  This excludes `60`, but it is explicitly a conjectural list and the authors
  say they do not know whether it is exhaustive.
  <https://arxiv.org/abs/2602.21047>

Consequently, a generic `End(J_bar)=Z` / `USp(4)` search is strategically more
plausible than forcing extra endomorphisms.  This is a heuristic prioritization,
not a theorem ruling out all RM, CM, or other special loci.

## LMFDB audit

As of the audit date, the LMFDB genus-2 database contains `66,158` curves in
`65,534` isogeny classes, all with absolute discriminant at most `1,000,000`.
Its torsion-order statistics list orders only through `39`, and the API query
for torsion order `60` returns an empty data array.

- Statistics: <https://www.lmfdb.org/Genus2Curve/Q/stats>
- Exact API query:
  <https://www.lmfdb.org/api/g2c_curves/?torsion_order=60&_format=json>

This proves only that the bounded LMFDB collection has no recorded order-60
curve; it says nothing about curves beyond that collection.

## Independent screen of Sutherland's 5-smooth-discriminant list

### Input

Andrew Sutherland supplies a second list of `487,493` global minimal genus-2
models over `Q` whose absolute discriminants are 5-smooth:

- Dataset description: <https://math.mit.edu/~drew/genus2curves.html>
- Raw data:
  <https://math.mit.edu/~drew/gce_genus2_hyperelliptic_5smooth.txt>

The file downloaded on 2026-07-10 had:

```text
lines       = 487493
bytes       = 34216695
sha256      = bb26b54405cb266296e21a5d280da129863a278744664c31745070659f37803d
```

Each line gives a model

```text
y^2 + h(x)*y = f(x)
```

and the current raw file also carries a Sato--Tate label as a fourth
colon-delimited field.

### Exact local test

At an odd good prime `p`, complete the square:

```text
Y^2 = h(x)^2 + 4*f(x).
```

The screen counted the curve over `F_p` and `F_(p^2)` directly.  If these
counts are `N1` and `N2`, then

```text
#J(F_p) = (N1^2 + N2)/2 - p.
```

Because every discriminant in the input is 5-smooth, every prime
`p >= 7` used below is good.  Also `p` does not divide `60`.  Hence a rational
`Z/60Z` would inject into `J(F_p)`, and the necessary condition

```text
60 | #J(F_p)
```

must hold at every tested prime.  The computation used exact finite-field
arithmetic; there were no bad or skipped records.

### Survivor counts

```text
prime just tested     survivors
-------------------   ---------
input                  487493
7                       22764
11                       1736
13                        205
17                         50
19                         25
23                         18
29                         18
31                         16
37                         16
41                         16
43                         16
```

All `25` survivors after prime `19` carry the dataset label `G_{3,3}`; the
later `16` survivors do as well.  In the genus-2 Sato--Tate classification,
`G_{3,3}` is the `SU(2) x SU(2)` decomposable type, so none is geometrically
simple.  Relevant primary references are:

- Sutherland's label/distribution page:
  <https://math.mit.edu/~drew/g2SatoTateDistributions.html>
- Fité--Kedlaya--Rotger--Sutherland, *Sato-Tate distributions and Galois
  endomorphism modules in genus 2*: <https://arxiv.org/abs/1110.6638>

The local test deliberately has false positives: divisibility of group orders
does not assert the presence of a compatible rational torsion point.  As a
spot check, the first survivor is LMFDB curve `4500.a.108000.1`.  LMFDB records
its actual rational torsion as `[6]`, its geometric endomorphism algebra as
`Q x Q`, and `is_simple_geom=false`:

- <https://www.lmfdb.org/api/g2c_curves/?label=4500.a.108000.1&_format=json>

Thus this screen rules out a geometrically simple `Z/60` **within this public
487,493-curve collection**, using the supplied Sato--Tate classifications.  It
does not establish a global nonexistence result.

## Direct audit of Nicholls Example 3.1.8

Nicholls' curve is

```text
C: y^2 = f30(x),

f30 = x^6 - (16/3)*x^5 + (70/9)*x^4
      + (131/27)*x^2 + (16/27)*x + 64/81.
```

The square change `Y=9y` gives the integral model

```text
Y^2 = 81*x^6 - 432*x^5 + 630*x^4 + 393*x^2 + 48*x + 64.
```

### Exact torsion and rational 2-torsion

Over `Q`, the sextic factors as

```text
(x^2 + 1/3)
*
(x^4 - (16/3)*x^3 + (67/9)*x^2 + (16/9)*x + 64/27),
```

with both factors irreducible.  The two Galois orbits of Weierstrass points
have sizes `2` and `4`, so there is exactly one nonzero rational 2-class:

```text
T2 = [x^2 + 1/3, 0].
```

Magma computes the exact rational torsion subgroup as

```text
J(Q)_tors = Z/30Z.
```

For any generator `P30`, the unique element of order 2 satisfies

```text
15*P30 = T2.
```

### Halving obstruction

The curve has good reduction at `p=5`, and Magma gives

```text
J(F_5) = Z/30Z.
```

Its order-2 point is not in `2*J(F_5)`; equivalently, this finite group has no
element of order `4`.  Therefore `T2` cannot be twice a rational point over
`Q`.  Magma's independent exact call `IsDivisibleBy(T2,2)` also returns
`false`.  Consequently the Nicholls order-30 point cannot be lifted to order
`60` by halving.

### Geometric simplicity

At the good prime `13`, the Frobenius polynomial is

```text
Phi_13(T) = T^4 - 4*T^3 + 6*T^2 - 52*T + 169.
```

It is irreducible, its Galois group has order `8` and transitive-group
description `D(4)`, and for a root `pi` one obtains

```text
deg minpoly(pi^n) = 4,  n=2,...,12.
```

This is the D4/root-power certificate used throughout the workspace and
independently confirms that the Jacobian is geometrically simple.

### It is not the workspace contact-(5,6) family

Every smooth member of the simultaneous contact-5/contact-6 family has an odd
quintic model and therefore a `Q`-rational Weierstrass point at infinity.  The
Nicholls sextic has no linear factor.  Its two rational points at infinity are
unramified and hence are not Weierstrass points.  A `Q`-isomorphism preserves
Weierstrass points, so this curve is not a specialization of that family up to
`Q`-isomorphism.

The family-wide `Q_11` obstruction in
`notes/contact30_to_60_halving_obstruction.md` is therefore separate.  For the
Nicholls seed, the simpler good-prime obstruction at `p=5` is decisive.

### Reproduction

Run:

```text
magma code/nicholls_order30_to60_halving.m
```

The script asserts every torsion, local-divisibility, factorization, and
geometric-simplicity claim above.

## Recommendation

Do not spend further time halving this Nicholls seed or the workspace
contact-(5,6) family.  The literature and database evidence do not suggest a
hidden known simple example.  Continue the order-20-plus-3 lane first, with the
full order-12-plus-5 surface as the strongest independent backup.  For every
candidate, demand exact `TorsionSubgroup`, a class-specific local check (not
only `60 | #J(F_p)`), and a D4/root-power simplicity certificate.

## Update: an explicit exact split specialization

The HLP point `(t,u,y)=(1/3,-1,9)` has since been converted into the explicit
genus-2 curve

```text
y^2 = -185*(125*x^2-1728)*(2000*x^4-48525*x^2+294912).
```

Magma certifies that its full rational torsion subgroup is exactly `[60]`.
This strengthens the positive control from subgroup containment to an exact
torsion computation for one specialization.  It does not alter the main
conclusion of this audit: the sextic is even, hence has an elliptic involution,
and its Jacobian is geometrically split.  Details and a reproducible verifier
are in `notes/hlp_z60_explicit_control_2026_07_14.md` and
`code/hlp_z60_explicit_verify.m`.
