# The genus-2 discriminant quotient at `s=59/49`

This note analyzes the quotient

```text
C_D: z^2 = D_s(q),    s=59/49,
```

where `D_s(q)` is the discriminant of

```text
U^3+(1+q)U^2+sU-(1+q)(q-s).
```

The unconditional arithmetic calculation is in

```text
code/elkies22210_orbit12_genus2_discriminant.m.
```

The explicit fake-2-cover collection is in

```text
code/elkies22210_orbit12_genus2_two_cover_collection.m.
```

The latter script uses GRH class-group bounds, and every result depending
on that choice is labeled below.

## 1. Integral and globally minimal models

At `s=59/49`,

```text
D(q) = 4q^5 - (775/49)q^4 - (290/49)q^3
       + (84509/2401)q^2 + (11728/2401)q
       - 2752704/117649.
```

Put

```text
X=14q,    Y=1372z.
```

Then `C_D` has the primitive integral model

```text
Y^2 = 14X^5 - 775X^4 - 4060X^3 + 338036X^2
      + 656768X - 44043264.                         (1)
```

Magma's `MinimalWeierstrassModel`, followed by `ReducedModel`, gives the
globally minimal reduced model

```text
y^2+(x^2+x)y
  = 7x^5+86x^4-2228x^3-29613x^2+7680x+4096.       (2)
```

The direct map from the original quotient to (2) is

```text
x = 7q-8,
y = (343z-49q^2+105q-56)/2.
```

The minimal discriminant factors as

```text
-2^12 * 7^4 * 59^3 * 127 * 8191^3 * 2078789.
```

## 2. Known rational points

The five currently known points on (2) are

```text
infinity,
(0,64), (0,-64),
(3/7,2050/49), (3/7,-2080/49).
```

The first affine pair is the split seed

```text
q=8/7,    z=+-128/343.
```

The second is the first-radicand boundary

```text
q=s=59/49,    z=+-590/2401.
```

An exact Magma point search on (1) with `Bound=10^6` returned exactly
these five points.  This is a bounded search, not a proof that the list
is complete.

## 3. Jacobian arithmetic (unconditional)

Unconditional quintic-field class-group computations and 2-descent give

```text
2 <= rank J_D(Q) <= 3,
J_D(Q)_tors = 0,
Sel_2(J_D/Q) = (Z/2Z)^3.
```

The two displayed affine `q`-values give divisor classes

```text
D_seed     = [(X-16, 512)],
D_boundary = [(X-118/7, 2360/7)]
```

on the Jacobian of (1).  Their height regulator is

```text
21.30979127586030539760,
```

so they are independent.  Magma's genus-2 saturation routine leaves
their subgroup unchanged at every prime

```text
2,3,5,7,11,13,17,19,23,29,31.
```

This is a certified `p`-saturation statement for those primes, not a
certificate that the subgroup is globally saturated.

`HasSquareSha(J_D)` returns true.  Consequently, **conditional on
finiteness of Sha**, the alternating Cassels--Tate pairing and the
3-dimensional 2-Selmer group force

```text
rank J_D(Q) = 3.
```

This conditional parity is evidence, not part of the unconditional rank
interval.  Magma cannot certify the local L-series data at `p=2` for
this model (`v_2` of the minimal discriminant is 12), so no root-number
claim is made.

At `p=11`, the Frobenius polynomial is

```text
121T^4+11T^3-14T^2+T+1,
```

which is irreducible over `Q`.  Thus `J_D` is `Q`-simple; there is no
rational elliptic quotient that would turn this into elliptic Chabauty.

Since the proven rank lower bound equals the genus, classical Chabauty
on `C_D` is unavailable.

## 4. Fake 2-Selmer covering collection (GRH class-group bounds)

With

```text
A=Q[theta]/(F(theta)/14),
```

Magma's complete local test returns five fake-2-Selmer classes.  The
class-group enumeration used to assert that this is the entire set was
run with `SetClassGroupBounds("GRH")`.

Three distinct classes contain the known points.  In Magma's descent-map
convention they are represented by

```text
14,              infinity,
16-theta,        q=8/7,
118/7-theta,     q=59/49.
```

Two everywhere-locally-soluble classes remain unexplained.  Exact
representatives for them, and their defining quadrics, are printed by
the covering script.

For any representative `delta`, write

```text
delta*(u0+u1 theta+...+u4 theta^4)^2
     = c0+c1 theta+c2 theta^2+c3 theta^3+c4 theta^4.
```

The corresponding fake-cover quotient is

```text
c2=c3=c4=0
```

in `P^4`.  The script verifies that every one is a nonsingular complete
intersection of three quadrics, hence has degree 8 and genus 5.  Its map
to (1) is

```text
X=-c0/c1,
Y=+-sqrt(14 Norm(delta))*Norm(u0+...+u4 theta^4).
```

The full two-cover before quotienting by the lifted hyperelliptic
involution has genus 17.

`PointSearch` to height 100 finds no rational point on either unexplained
genus-5 cover.  It finds the expected point `[1:0:0:0:0]` on each of the
three known covers.  Again, these searches are diagnostic only.  There
is no local elimination: all five classes passed the two-cover local
conditions.

## 5. Consequence for the genus-7 attack

This quotient does not by itself finish the genus-7 curve:

1. its Jacobian rank is already at least its genus;
2. it is `Q`-simple, so there is no lower elliptic quotient;
3. two fake-cover classes are locally soluble but have no known rational
   points; and
4. even the three known genus-5 covers could contain further points.

The useful output is therefore the explicit five-cover collection.  A
rigorous completion needs one of:

- a full finite-index Mordell--Weil basis followed by a Mordell--Weil
  sieve on the five genus-5 covers;
- a second descent on those covers; or
- quadratic Chabauty using additional endomorphisms/correspondences, if
  such extra Neron--Severi classes can be proved.

The immediate arithmetic bottleneck is to decide whether the rank is 2
or 3 and, in the likely rank-3 case, to find a third Mordell--Weil
generator.
