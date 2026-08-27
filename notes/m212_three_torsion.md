# Algebraic `3`-torsion on `M(2,12)`

The `M(2,12)` chart already has a marked order-`12` class.  Thus it
automatically has rational `3`-torsion: if `D` is the marked order-`12` class,
then `4D` has order `3`.

I made this explicit in

```text
code/m212_construct_3torsion.m
```

Start from the `M(12)` model

```text
y^2 + (x-r)(T+1)y = a x^2 T(T+1),
T = a x^2 - x + r.
```

On the `M(2,12)` chart, write

```text
1 - 4a(r+1) = z^2,
a = (1-z^2)/(4(r+1)).
```

Choose the rational root

```text
w = 2(r+1)/(1+z)
```

of `T+1` and move it to infinity by `x = w + 1/X`.  The image of the marked
point `P=(0,0)` gives a divisor class `D` of order `12` on the odd quintic
model.

Magma computes the following generic Mumford representative for `4D`:

```text
4D = [u(X), v(X)]
```

where

```text
u(X) = X^2 + z X + (z^2 - 1)/(4(r+1))
```

and

```text
v(X) =
  (-4 z^3 r^2 - 3 z^3 r + 4 z^2 r^2 + 7 z^2 r + 2 z^2
   - 5 z r - 4 z + r + 2)/(4(z+1)(r+1)) * X
 +(-z^3 r + 2 z^2 r + z^2 - z r - 2 z + 1)/(4(r+1)).
```

The script verifies symbolically:

```text
3*(4D) = 0,
4D != 0.
```

Thus this gives an explicit algebraic rational `3`-torsion point on the
generic `M(2,12)` chart.  It is not an additional independent `3`-torsion
condition; it is the `3`-primary part already contained in the marked
order-`12` class.

The chosen chart excludes the denominators `z+1=0` and `r+1=0`.  The conjugate
root of `T+1` gives the corresponding formula on the complementary chart by
replacing `z` with `-z`.


## Independent `3`-torsion

The construction above is only the built-in `3`-part of the marked order-`12`
class.  To look for an independent rational `3`-torsion point, I added

```text
code/m212_extra3_search.m
```

The finite diagnostic uses the necessary condition that, at every good prime
`p != 3`, the finite group `J(F_p)[3]` must have rank at least `2`.  This is
stronger than merely requiring `3 | #J(F_p)`, since the built-in class already
accounts for one rational `3`-direction.

The first diagnostic gives an immediate open-chart obstruction at `p=5`:

```text
p 5 total 25 good 6 bad 19 rank0 0 rank1 6 rank2 0 rank3plus 0
```

Thus any rational point with an extra independent `3`-torsion class must reduce
to the bad/boundary locus modulo `5`.

A height-`30` rational search with rank filters through `p=43` found such a
specialization:

```text
z = -5/3,   r = -3/5,   a = -10/9.
```

The conjugate value `z=5/3` gives the same curve.  On the odd integral model

```text
Y^2 =
  5668704 X^5 - 22143375 X^4 + 36098622 X^3
  - 30305259 X^2 + 12990780 X - 2259900,
```

Magma computes

```text
TorsionSubgroup(J) = [3,12].
```

The verifier

```text
code/m212_extra3_verify_hit.m
```

prints the built-in order-`3` class

```text
4D = (X^2 - 5/3 X + 10/9, 486 X, 2)
```

and an independent order-`3` class

```text
E = (X^2 - 17/12 X + 19/36, -243/4 X - 81/4, 2).
```

It verifies

```text
3E = 0,
E != 0,
E != 4D,
E != -4D.
```

So this is an explicit `M(2,12)` specialization with two independent rational
`3`-torsion directions, namely torsion `[3,12]`.


## Geometric simplicity of the extra-`3` hit

The extra-`3` specialization above is not geometrically simple.  I added

```text
code/m212_extra3_geom_simple_check.py
```

which runs Sage/Lombardo's endomorphism tests and compares Frobenius
polynomials with elliptic curves in the Cremona database.  Sage reports

```text
geometric_endomorphism_algebra_is_field(B=100) = False
geometric_endomorphism_ring_is_ZZ(B=100)       = False
```

Moreover, the genus-2 Frobenius polynomial matches the product of the
Frobenius polynomials of elliptic curves in the isogeny classes

```text
90c:  y^2 + x y + y = x^3 - x^2 + 13 x - 61,
510g: y^2 + x y     = x^3 + 25 x - 375.
```

The verifier checks this equality at all good primes tested up to `293`:

```text
P_p(J,T) = P_p(90c,T) P_p(510g,T).
```

So the extra-`3` hit is a split/nonsimple example.  More explicitly,

```text
code/m212_extra3_split_certificate.m
```

uses Magma's `Degree2Subcovers` to produce two degree-`2` maps from the
hyperelliptic curve to elliptic curves.  After minimization, the elliptic
factors are

```text
90c3:  y^2 + x y + y = x^3 - x^2 - 122 x + 1721,
510g1: y^2 + x y     = x^3 + 25 x - 375.
```

Thus the Jacobian is split over `Q`, with elliptic factors in Cremona classes
`90c` and `510g`.  This is therefore a positive example for the split target,
not a simple one.


## Target correction: simple examples only

The split certificate above is useful for understanding the construction, but it
is not a target example if we require geometric simplicity.  Under the corrected
criterion, the specialization

```text
z = -5/3, r = -3/5, a = -10/9
```

should be discarded: its Jacobian is split over `Q` with elliptic factors in
Cremona classes `90c` and `510g`.

I added a faster residue-table search

```text
code/m212_extra3_residue_search.m
```

which precomputes, for each good prime, the `(z,r) mod p` residues where
`J(F_p)[3]` has rank at least `2`.  This replaces the slower repeated finite
Jacobian group computations in `code/m212_extra3_search.m`.

A partial height-`30`, prime-bound-`43` run reproduced the split extra-`3` hit
and exact-checked subsequent residue survivors back down to torsion `[12]`.
The output is in

```text
data/m212_extra3_residue_h30_p43.txt
```

The only exact `[3,12]` hit in the partial run was the discarded split curve:

```text
EXTRA3_HIT z=-5/3 r=-3/5 a=-10/9 torsion [3,12]  split
```

The later exact tests printed before interruption all had torsion `[12]`, for
example:

```text
z=-29/4 r=23/14     torsion [12]
z=-23/5 r=-17/6     torsion [12]
z=-17/14 r=19/11    torsion [12]
z=-26/15 r=13/18    torsion [12]
z=-11/17 r=-12/25   torsion [12]
z=-11/19 r=26/29    torsion [12]
z=-1/21 r=-1/30     torsion [12]
z=-29/22 r=23/5     torsion [12]
z=-30/23 r=-7/29    torsion [12]
z=-6/23 r=10/11     torsion [12]
```

So, with the simple-only target, no acceptable `M(2,12)+extra3` example has
been found yet.  The next search should keep the residue-table sieve but add a
geometric-simple screen for every exact `[3,12]` hit.
