# M(12) route with simple Jacobians

This follows the proposed five-step route, but with the extra requirement that the Jacobian be simple.

## 1. Completed-square model

Start from Theorem 9.4's `M(12)` model

```text
y^2 + (x-r)(T+1)y = a x^2 T(T+1),
T = a x^2 - x + r,
P = (0,0).
```

Completing the square with `Y = 2y + (x-r)(T+1)` gives

```text
Y^2 = W(x),
W = ((x-r)(T+1))^2 + 4a x^2 T(T+1).
```

Magma factors this as

```text
W = (T+1) Q4,
```

where

```text
Q4 = (4a^2+a)x^4 + (-2ar-4a-1)x^3
     + (ar^2+4ar+3r+1)x^2 + (-3r^2-2r)x + r^3+r^2.
```

## 2. Rational Weierstrass point

The quadratic factor `T+1` has discriminant

```text
1 - 4a(r+1).
```

So a rational Weierstrass point is obtained by writing

```text
1 - 4a(r+1) = z^2,
a = (1-z^2)/(4(r+1)).
```

For each rational root `w` of `T+1`, the script moves `w` to infinity via

```text
X = 1/(x-w),   Y_old = Y_new/X^3,
```

and verifies that the image of `P` gives a Jacobian point of order `12`.

## 3. Simple-Jacobian certificate

Magma's `Jacobian(C)` object here does not support `IsSimple` directly.  The script instead uses a sufficient certificate: find a good prime `p` such that the genus-2 `LPolynomial(C mod p)` is irreducible of degree 4 over `Q`.  If the Jacobian were `Q`-isogenous to a product of elliptic curves, the Frobenius polynomial would factor as a product of two quadratics at every good prime.  Thus this certifies `Q`-simplicity.

## 4. Results

`code/m12_simple_search.m` finds simple examples immediately.

A representative example is

```text
a = 13/2,   r = -25,
T = (13/2)x^2 - x - 25.
```

The model is

```text
y^2 + (13/2*x^3 + 323/2*x^2 - 49*x - 600)y
 = 2197/8*x^6 - 169/2*x^5 - 8255/4*x^4 + 637/2*x^3 + 3900*x^2.
```

The completed-square polynomial factors as

```text
W = (x - 2)(x + 24/13)(x + 50/27)
    (x^3 - 2/13*x^2 + 251/13*x - 600/13).
```

So the curve has three rational Weierstrass points.  Taking `w=2` and moving it to infinity gives the odd quintic

```text
Y^2 = -2600 X^5 + 517349 X^4 + 368149 X^3
      + 413001/4 X^2 + 30901/2 X + 4563/4.
```

Equivalently, after scaling `Y`,

```text
Y^2 = -10400 X^5 + 2069396 X^4 + 1472596 X^3
      + 413001 X^2 + 61802 X + 4563.
```

Magma verifies:

```text
Order(P - infinity) = 12.
```

There is an independent rational 2-torsion point from the extra rational root.  In the odd model above, the rational roots are

```text
X = -13/50,   X = -27/104.
```

The first is `6(P-infinity)`, while the second is independent.  Thus the Jacobian contains at least

```text
Z/12Z x Z/2Z.
```

The simplicity certificate for this example is at `p=11`:

```text
L_p(T) = 121 T^4 + 11 T^3 - 14 T^2 + T + 1,
```

which is irreducible.


## 4a. One-parameter `Z/12 x Z/2` family

The extra rational root condition has a simple component.  After imposing

```text
1 - 4a(r+1) = z^2,
```

the equation `Q4(u)=0` factors on the line `z=±r`; one factor is

```text
r u + 2r - 2u = 0.
```

Thus

```text
a = (1-r)/4,
u = -2r/(r-2)
```

gives a one-parameter family with three rational Weierstrass points.  In this family

```text
T+1 = (x-2)(x + (2r+2)/(r-1))
```

and

```text
W = (x-2)(x + 2r/(r-2))(x + (2r+2)/(r-1))
    (x^3 + 4/(r-1)x^2 + (-r^2-5r-2)/(r-1)x
       + (2r^2+2r)/(r-1)).
```

The script `code/m12_z12x2_family.m` verifies the family and the sample `r=-25`.

## 5. Full split attempt

Requiring all six Weierstrass points to be rational did not produce hits in the parametrized search.  More importantly, `code/m12_full_split_finite_field_sieve.m` shows a good-reduction obstruction modulo `3`:

```text
p 3 good 1 full_split 0
```

So any rational full-split point in this `M(12)` family would have to reduce to the boundary modulo `3`.  This is the same pattern seen in the split `M(2,2,2,6)` halving attempts.

Conclusion: this route does produce strictly simple Jacobians with larger torsion, concretely at least `Z/12 x Z/2`.  The stronger full-split `(2,2,2,12)` target still appears obstructed away from boundary.
