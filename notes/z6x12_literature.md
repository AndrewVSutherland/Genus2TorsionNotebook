# Literature and construction audit for `[6,12]`

Audit date: 2026-07-10.

## Conclusion

No geometrically simple genus-2 Jacobian over `Q` with rational torsion
`Z/6Z x Z/12Z` was found in the literature or public databases audited here.
This is a scoped search result, not a nonexistence theorem.

Howe--Leprévost--Poonen (HLP) give a positive-rank family whose Jacobians
contain `Z/6Z x Z/12Z`, but every member is `(2,2)`-isogenous to a product of
elliptic curves.  Thus the known published construction is entirely on the
decomposable/Humbert locus.

The most direct simple construction is now very concrete: start with the
workspace's geometrically simple exact `[6,6]` curve and halve one of its
three nonzero rational 2-classes.  Equivalently, inspect its rational Richelot
neighbors for gained 4-torsion.  Either operation preserves the full rational
`[3,3]` subgroup and changes the 2-primary part from `[2,2]` to `[2,4]`, giving
exactly the desired invariant pattern `[6,12]`.

## Published split family

HLP Section 3.5 takes two elliptic curves with rational 12-torsion.  Their
2-torsion modules can be glued when

```text
(2*t^2-2*t+1)*(6*t^2-6*t+1)*y^2
  = (2*u^2-2*u+1)*(6*u^2-6*u+1).                (1)
```

Taking `(u,y)=(t,1)` as the zero section, HLP show that `(u,y)=(t,-1)` is an
infinite-order section.  Away from the usual degeneracy/isomorphism locus,
Proposition 4 glues `E_12^t` and `E_12^u` along 2-torsion and produces a
genus-2 Jacobian containing

```text
(Z/12Z x Z/12Z)/<(6P,6Q)> = Z/6Z x Z/12Z.
```

Every such Jacobian is `(2,2)`-isogenous to `E_12^t x E_12^u`, so it is
geometrically split.  HLP prove subgroup containment; they do not assert that
the full rational torsion is exactly `[6,12]`.

Primary source:

- Everett W. Howe, Franck Leprévost, Bjorn Poonen, *Large torsion subgroups of
  split Jacobians of curves of genus two or three*, Forum Math. 12 (2000),
  315--364: <https://math.mit.edu/~poonen/papers/large.pdf>

The relevant items are Tables 6--7, Proposition 4, and Section 3.5.

## A fully explicit HLP split seed

The paper does not print a specialized genus-2 equation for this entry, but
one follows cleanly from its formulas.  Set

```text
t = u = 1/3,    y = -1
```

in (1).  After the rational square scaling from HLP Table 6, both elliptic
factors are

```text
E: Y^2 = X*(X^2 - 183*X + 9216),
P = (24,360),    Order(P)=12,    6P=(0,0).
```

The two nonzero nonrational 2-points have `X`-coordinates

```text
(183 +/- 15*sqrt(-15))/2.
```

The gluing fixes `(0,0)` and swaps these two conjugate points.  Applying HLP
Proposition 4 and removing a rational square factor gives

```text
C_HLP:
y^2 = 183*(x^2+1)*(32*x^2+61*x+32)*(32*x^2-61*x+32).
```

This equation was derived from Proposition 4; it is not printed verbatim in
HLP.  It has two explicit degree-2 quotient maps to `E`:

```text
phi_1(x,y) = (
  3072/61*(x^2+1),
  3072/3721*y
),

phi_2(x,y) = (
  3072/61*(1/x^2+1),
  3072/3721*y/x^3
).
```

Direct substitution verifies both maps.  Pulling back `P=(24,360)` gives
rational degree-2 divisor classes supported by

```text
phi_1^{-1}(P):
  x^2 + 67/128 = 0,    y = 55815/128,

phi_2^{-1}(P):
  x^2 + 128/67 = 0,    y = -(55815/67)*x.
```

Subtracting the corresponding degree-2 fibers over the elliptic origins gives
two rational classes of order 12.  Their sixth multiples agree, so the group
they generate has invariants `[6,12]`.  This is a compact split seed with
marked torsion data for deformation away from the Humbert locus.

The independent exact Magma check

```text
code/contact6_m612_hlp_seed_verify.m
```

gives

```text
J(C_HLP)(Q)_tors = [6,12]
```

of order `72`.  Its good-prime local polynomials are squares of elliptic
quadratics, as expected from the geometric splitting.

## Public database audit

### LMFDB

The current LMFDB genus-2 collection has `66,158` curves of absolute
discriminant at most `1,000,000`.  Its torsion-order statistics contain no
order `72`, and the exact API query returns an empty data array:

- <https://www.lmfdb.org/Genus2Curve/Q/stats>
- <https://www.lmfdb.org/api/g2c_curves/?torsion_order=72&_format=json>

This only excludes the bounded LMFDB collection.

### Sutherland's 487,493-curve 5-smooth list

The larger public list is described at
<https://math.mit.edu/~drew/genus2curves.html>; raw data are at
<https://math.mit.edu/~drew/gce_genus2_hyperelliptic_5smooth.txt>.

For each model `y^2+h*y=f`, the screen completed the square and counted the
curve over `F_p` and `F_(p^2)`.  If the counts are `N1,N2`, then

```text
#J(F_p) = (N1^2+N2)/2-p.
```

All primes `p >= 7` are good because the input discriminants are 5-smooth.
Since `p` does not divide `72`, a rational subgroup of order `72` injects on
reduction, so `72 | #J(F_p)` is necessary.  The exact survivor counts were:

```text
input     487493
p=7       19627
p=11       2018
p=13        321
p=17        119
p=19         52
p=23         46
p=29         43
p=31         36
p=37         36
p=41         33
p=43         33
```

The 33 final records have supplied Sato--Tate labels

```text
G_{3,3}       28
E_1            2
J(E_1)         2
N(G_{1,3})     1
```

All four are geometrically decomposable types in the genus-2 Sato--Tate
classification.  Thus the necessary local screen leaves no geometrically
simple candidate in this public collection.  Divisibility is only a necessary
condition, so these 33 are expected to include many curves whose actual
torsion is much smaller.

Classification references:

- <https://math.mit.edu/~drew/g2SatoTateDistributions.html>
- Fité--Kedlaya--Rotger--Sutherland,
  <https://arxiv.org/abs/1110.6638>

## Theoretical constraints

There is no elementary symplectic obstruction to `[6,12]`.  Its primary
parts are

```text
3-primary: (Z/3Z)^2,
2-primary: Z/2Z x Z/4Z.
```

The rational 3-subgroup can be maximal isotropic: the Weil pairing takes
values in `mu_3`, and `Q` contains no nontrivial cube roots of unity.  The
2-primary pattern is also compatible with a principally polarized abelian
surface.  The target is therefore structurally plausible.

There are useful endomorphism-type restrictions:

- Laga--Schembri--Shnidman--Voight prove that an abelian surface over `Q`
  whose geometric endomorphism ring is a maximal order in a non-split
  quaternion algebra has rational torsion killed by 12 and of order at most
  18.  This rules out `[6,12]` on that maximal-order PQM locus.
  <https://arxiv.org/abs/2308.15193>
- Alessandrì--Coppola's Conjecture 4.5 does not include torsion order `72`
  among the predicted orders for simple non-CM `GL_2`-type abelian surfaces.
  This is conjectural evidence only; the authors explicitly do not claim an
  exhaustive theorem.  <https://arxiv.org/abs/2602.21047>

These results favor a generic `End(J_bar)=Z` construction rather than one
that deliberately forces extra endomorphisms.

## Workspace assets and recommended attacks

### 1. Halve a rational 2-class on the simple `[6,6]` seed

The best existing seed is

```text
C_66:
y^2 = 11389248*x^5 - 18252000*x^4 + 42399396*x^3
      - 10288044*x^2 + 29659500*x,

J(Q)_tors = [6,6].
```

It has contact parameters

```text
a=133/39, b=-7/13, L=29/16, U=-9/4, v=5/2,
```

The stronger verifier `code/contact6_m36_simple_hit_strong_verify.m` gives a
D4 geometric-simplicity certificate at `p=37`.
Its squarefree polynomial is
proportional to

```text
x*(48*x^2+8*x+39)*(39*x^2-69*x+125),
```

where both quadratics are irreducible.  Hence

```text
J[2](Q) = (Z/2Z)^2
```

with three explicit nonzero classes furnished by the three factors.

Because the odd-order part of multiplication by 2 is invertible, an order-6
class is divisible by 2 exactly when its order-2 part is.  Therefore it is
enough to test the three rational 2-classes `T` for

```text
T in 2*J(Q).
```

If any test succeeds, the existing `[3,3]` survives and the torsion contains
`[6,12]`.  The efficient workflow is class-specific local divisibility in
`J(F_p)` first, followed only for survivors by the compact quadratic-factor
halving cover and exact `IsDivisibleBy`.

Source: `notes/contact6_m36.md`, especially the `[1,2,2]` core-cover section.

### 2. Test rational Richelot neighbors

The factorization `x*q1*q2` supplies a rational maximal isotropic subgroup of
`J[2]`.  A degree-4 Richelot isogeny is prime to 3, so it preserves the full
rational `[3,3]`.  If a Richelot codomain has 2-primary torsion `[2,4]`, its
full torsion has the desired invariants `[6,12]`.

Prepared scripts are:

```text
code/contact6_m612_richelot_seed_test.m
code/contact6_m612_richelot_known_hits.m
```

The exact seed sweep now gives one rational Richelot neighbor of the simple
`[6,6]` seed.  Its torsion remains `[6,6]`; the complete rational 2-power
isogeny traversal likewise contains only one Jacobian, again with `[6,6]`.
The three smaller core seeds map to elliptic products.  A finite diagnostic
finds no good-open `[6,12]` source or distinguished dual at `p=5` or `p=7`,
although target-compatible residues are abundant from `p=11` onward.  Thus
this particular lane is forced onto simultaneous small-prime boundary disks.

### 3. Add independent 3-torsion to the simple `[2,12]` family

The workspace has a one-parameter family

```text
a=(1-r)/4
```

with torsion containing `Z/12Z x Z/2Z` and simple specializations.  In the
`M(2,12)` coordinates it is the line `z=+/-r`.  Adding an order-3 class
independent from the built-in `4D` produces

```text
Z/12Z x Z/2Z x Z/3Z = Z/6Z x Z/12Z.
```

The general `M(2,12)+extra3` search has a good-open obstruction modulo 5; its
only known exact `[3,12]` point is split and is not on this `[2,12]` line.
Thus the right next calculation is a one-variable, boundary-aware extra-3
search on `z=+/-r`, not another two-variable blind scan.

Sources:

```text
notes/m12_simple_route.md
notes/m212_three_torsion.md
code/m12_z12x2_family.m
code/m212_extra3_residue_search.m
```

### 4. Deform the explicit HLP seed

The displayed `C_HLP` has three rational quadratic factors and two marked
order-12 pullback classes.  It is therefore a convenient exact base point for
linearizing the marked-torsion equations and seeking a transverse direction
off the decomposable Humbert divisor.  Any proposed deformation should track
the two divisor classes explicitly and certify that the resulting Frobenius
quartic is absolutely irreducible; merely perturbing the coefficients will
not preserve torsion.

## Priority

1. Run class-specific finite halving tests on the three rational 2-classes of
   `C_66`.
2. In parallel, compute the distinguished Richelot codomain torsion.
3. If both fail locally, specialize the existing extra-3 residue machinery to
   the simple `[2,12]` line `z=+/-r`, with a full `p=5` boundary analysis.
4. Use the explicit HLP curve and marked pullbacks for a deformation attack
   only after these cheaper exact tests.
