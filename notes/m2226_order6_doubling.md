# M(2,2,2,6): order-6 classes and theta doubling

This uses the cleaner form of Theorem 9.4.  Write the even-degree model as

```text
y^2 = x (x + 2s^2 - sn)
        (x + 2s^2 + sm - 2sn - mn)
        (x + 2s^2 + sm - sn - mn)
        (2x - mn)
        (2x + 4s^2 - 4sn - mn).
```

Send the Weierstrass point `x=0` to infinity by

```text
X = 1/x,    Y = y/x^3.
```

Then

```text
Y^2 = L1 L2 L3 L4 L5,    P = (0,2),
```

where infinity is the old `x=0` Weierstrass point and

```text
L1 = 1 + (2s^2 - sn) X
L2 = 1 + (2s^2 + sm - 2sn - mn) X
L3 = 1 + (2s^2 + sm - sn - mn) X
L4 = 2 - mn X
L5 = 2 + (4s^2 - 4sn - mn) X.
```

The class `g = P - infinity` has order 6.

## Other order-6 classes

Let `A_i = L_i(0)`, so `(A_1,A_2,A_3,A_4,A_5) = (1,1,1,2,2)`.  Let `W_i` be the finite Weierstrass point cut out by `L_i=0`.

The six classes represented by `P +` a Weierstrass point are the obvious ones and are not useful for the theta-doubling test.  The other ten coset representatives are

```text
g + (W_i + W_j - 2 infinity),    1 <= i < j <= 5.
```

For these, the quadratic through `P,W_i,W_j` is

```text
ell_ij = 2 L_i L_j/(A_i A_j).
```

Put

```text
raw_U_ij = ( product_{k notin {i,j}} L_k
             - 4 L_i L_j/(A_i^2 A_j^2) ) / X.
```

Then the reduced Mumford representative is

```text
[ monic(raw_U_ij),  -ell_ij mod monic(raw_U_ij) ].
```

The condition that this class is represented by `2Q` for a point `Q` on the curve is

```text
Disc_X(raw_U_ij) = 0.
```

The pair `[3,5]` is exceptional: `raw_U_35` is already a scalar multiple of `X^2`, and Magma checks that the corresponding class has order 3, not order 6.  This is the expected `g + 3g = 4g` type exception in the coset.

## Rational-point result

After removing boundary factors from the discriminants, the remaining nine curves in `P^2_{s,m,n}` are genus-one curves.  Magma maps each one to the rank-zero elliptic curve

```text
y^2 = x^3 + 4x,
```

with torsion `Z/4Z`.  Pulling back all torsion points, and adding the rational singular/base points of the plane models, gives only points on the bad locus of the original family.  The bad locus factors as

```text
n, m, s, s-n, s-n/2, s+m/2,
sm-sn-mn,
sm+sn-mn/2,
s^2-sn/2+mn/4,
s^2+sm/2-sn-mn/4.
```

So this route gives no nondegenerate rational point: none of the other order-6 classes in the `M(2,2,2,6)` family is `2Q` for a rational point `Q`, away from the singular/degenerate boundary.

## Scripts

- `code/m2226_order6_doubling.m` prints the quintic model, the compact Mumford formulas, and the discriminant factors, then verifies the formulas on two rational specializations using Jacobian arithmetic.
- `code/m2226_search_doubling_curves.m` does a bounded primitive-triple search with nonsingularity and order-6 checks.
- `code/m2226_certify_doubling_curves.m` certifies the residual genus-one curves by mapping them to the rank-zero elliptic curve above and checking that all rational points are boundary points.
