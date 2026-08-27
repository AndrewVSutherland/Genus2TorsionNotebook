# Why HLP 2-gluing cannot produce `[3,21]`

## Exact scope

Let `E` and `F` be elliptic curves over `Q`, let

```text
psi : E[2] -> F[2]
```

be a Galois-equivariant isomorphism, and let `G` be its graph.  Assume that
`psi` is not induced by an isomorphism of the elliptic curves, so that the HLP
quotient

```text
J = (E x F)/G
```

is the principally polarized Jacobian of a genus-2 curve.  Thus the quotient
map is defined over `Q`, has degree `4`, and has kernel contained in
`(E x F)[2]`.

Under precisely these hypotheses, `J(Q)` cannot contain

```text
Z/3Z x Z/21Z.
```

This also applies to any abelian surface over `Q` connected to a product of
two elliptic curves over `Q` by a `2`-power-degree isogeny.  It is not a
statement about a product decomposition that exists only over `Qbar`, nor
about isogenies whose degree is divisible by `3` or `7`.

## Proof

Let

```text
pi : E x F -> J
```

be the HLP `(2,2)`-isogeny.  For every odd integer `n`, the restriction

```text
pi : (E x F)[n] -> J[n]
```

is injective because its kernel is `2`-primary.  The two finite group schemes
have the same order `n^4`, so this is a Galois-equivariant isomorphism.  Taking
Galois invariants, and then taking the union over odd `n`, gives

```text
J(Q)[odd] = E(Q)[odd] x F(Q)[odd].
```

Mazur's theorem says that the rational torsion of an elliptic curve over `Q`
is either cyclic of order `1 <= N <= 10` or `N=12`, or is
`Z/2Z x Z/2NZ` for `1 <= N <= 4`.  Consequently the odd part of the rational
torsion of one elliptic curve is one of

```text
0, Z/3Z, Z/5Z, Z/7Z, Z/9Z.
```

Now

```text
Z/3Z x Z/21Z = (Z/3Z)^2 x Z/7Z.
```

Two independent rational `3`-directions on `J` force both elliptic factors to
have nonzero rational `3`-torsion.  A rational `7`-torsion point on `J` forces
at least one of the two factors also to have a rational `7`-torsion point.
That factor would therefore have a rational point of order `21`, contradicting
Mazur's theorem.  Hence `[3,21]` is impossible in this HLP chart.

Equivalently, the only way to distribute `(Z/3)^2 x Z/7` between two cyclic
elliptic odd-torsion groups would put both a `3` and a `7` on the same factor.

## Why `[63]` is different

There is no analogous obstruction to cyclic `[63]`, since

```text
Z/63Z = Z/9Z x Z/7Z
```

and elliptic curves over `Q` may separately have rational points of orders `9`
and `7`.  HLP carry this out explicitly.  Their Section 3.6, equation (4), is

```text
y^2 = 897*x^6 - 197570*x^4 + 79136353*x^2 - 146398496.
```

Its Jacobian is `(2,2)`-isogenous over `Q` to the product of an elliptic curve
with `9`-torsion and one with `7`-torsion.  HLP prove that the full rational
torsion is cyclic `[63]`, using good reduction at `5`, where the Jacobian has
exactly `63` points.  The local exact verifier is

```text
code/z63_hlp_exact_verify.m
```

A Richelot or any other `2`-power isogeny cannot turn this cyclic `3`-primary
part `Z/9` into `(Z/3)^2`: all such isogenies preserve the complete odd
Galois module.  Changing `[63]` into `[3,21]` would require a construction with
degree divisible by `3`, such as a genuinely different `3`-gluing problem.

## Realistic route to `[3,21]`

The Lepr\'evost `[21]` family is the natural route because it is generically
geometrically simple; its `t=1` specialization has geometric endomorphism ring
`Z`.  The split obstruction above therefore does not apply.

Write `T=7*D21`, the built-in rational class of order `3`.  We seek a rational
class `E` of order `3` with `E` not in `<T>`.  At every good prime `p != 3`,
reduction is injective on rational `3`-torsion, so a necessary condition is

```text
dim_F3 J_t(F_p)[3] >= 2.
```

This gives a very cheap residue sieve on the one parameter `t`; it is
implemented in

```text
code/z21_extra3_residue_search.m.
```

After the finite sieve, exact `TorsionSubgroup` computations are needed only
for the rare survivors.  A hit must be checked for exact invariants `[3,21]`
and, if a nonsplit example is desired, by a geometric-simplicity certificate.

There is also a direct algebraic version on the normalized contact-7 model
`y^2=g_t(z)`.  Generically, a pair `{E,-E}` of nonzero `3`-torsion classes is
encoded by

```text
q = z^2 + U*z + V,
H = z^3 + A*z^2 + B*z + C,
H^2 - q^3 = L^2*g_t.
```

Comparing the six coefficients gives the `3`-division cover of the `t`-line.
There are `(3^4-1)/2 = 40` nonzero pairs `{E,-E}` geometrically.  One is the
known pair `{T,-T}`; the residual degree-`39` cover is exactly where an
independent rational `3`-class must lie (apart from lower-degree boundary
representatives).  If the height sieve finds persistent residue branches,
factoring and analyzing this residual cover is the next exact step.

## Prior local status

There was no previous direct `[3,21]` search or example in the repository.
The older `contact7_root_plus3` computation forced a rational Weierstrass root,
hence rational `2`-torsion, and searched only for the order-`42` necessary
condition.  It did not test `3`-rank two and is not a `[3,21]` attempt.

## Primary reference

Everett W. Howe, Franck Lepr\'evost, and Bjorn Poonen, *Large torsion
subgroups of split Jacobians of curves of genus two or three*, Forum Math. 12
(2000), 315--364:

<https://math.mit.edu/~poonen/papers/large.pdf>

Propositions 3 and 4 give the `(2,2)` quotient and its degree-`4`, `2`-torsion
kernel; Section 2 states the needed form of Mazur's theorem; Section 3.6 gives
the exact cyclic `[63]` curve.
