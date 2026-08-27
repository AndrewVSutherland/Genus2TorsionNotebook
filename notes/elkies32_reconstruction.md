# Reconstructing the Elkies `[32]` Family

Elkies' genus-2 torsion page lists a one-parameter family with a rational
point of order `32`, but the displayed page only gives a single member:

```text
y^2 = (15*x - 1)*(1056*x^4 + 156183*x^3 + 26297*x^2 + 649*x - 121).
```

The goal here was to reconstruct the hidden family from the printed member,
without doing a blind coefficient search.

## Normalized contact model

The printed curve has rational points `(0,11)` and `(-1,1440)`.  After scaling
`y` by `11`, it fits the model

```text
y^2 = h(x)^2 - a^2*x^5*(x+1),
h(x) = a*x^3 + b*x^2 + c*x + 1.
```

For this model,

```text
div(y-h) = 5*(0,1) + (-1,h(-1)) - 6*infinity.
```

Thus, if

```text
D = (0,1) - infinity,
E = (-1,h(-1)) - infinity,
```

then `E = -5D` in the Jacobian.

For Elkies' printed member the normalized parameters are

```text
a = -240,
b = -1323/11,
c = -112/11,
p = a - b + c - 1 = -1440/11,
r = 1/15.
```

Here `p=-h(-1)` and `r` is the finite Weierstrass root.

## Order-32 condition

If `W=(r,0)` and `T=W-infinity`, then `T` is a rational 2-torsion class.  To
force `T=16D`, use the relation

```text
4D + 4E + T = 0.
```

Since `E=-5D`, this relation is equivalent to `T=16D`.

A function with divisor

```text
4*(0,1) + 4*(-1,h(-1)) + W - 9*infinity
```

must lie in `L(9*infinity)`, with basis

```text
1, x, x^2, x^3, x^4, y, x*y, x^2*y.
```

So the order-32 condition is a concrete `9 x 8` rank drop, plus `f(r)=0`.

## Reduction

The common non-boundary part of the first rank condition becomes linear after
setting

```text
p = a - b + c - 1.
```

It gives

```text
c = (-a^4 + 4*a^3*p - 8*a^2*p^2 + 8*a*p^3
     + 8*a^2*p - 16*p^3)/(4*a^2*p),
b = a + c - 1 - p.
```

Set `z=a/p`.  Eliminating `p` from `f(r)=0` and one remaining rank minor gives
the component containing Elkies' printed point:

```text
0 = 3*z^4*r^4 + 9*z^4*r^3 - 16*z^3*r^4
  + 10*z^4*r^2 - 56*z^3*r^3 + 32*z^2*r^4
  + 5*z^4*r - 72*z^3*r^2 + 144*z^2*r^3 - 64*z*r^4
  + z^4 - 40*z^3*r + 208*z^2*r^2 - 224*z*r^3 + 80*r^4
  - 8*z^3 + 120*z^2*r - 288*z*r^2 + 160*r^3
  + 24*z^2 - 160*z*r + 160*r^2 - 32*z + 80*r + 16.
```

The printed example is the rational point

```text
(z,r) = (11/6, 1/15).
```

The scale `p` is the common root of two quadratic equations in `p`; the script
records the exact rational formula.  At the printed point it gives
`p=-1440/11`, recovering the displayed curve.

## Geometry of the component

Magma verifies that the projective plane closure of this component is
irreducible of genus `0`.  Its rational singular points are

```text
(2 : 0 : 1)  multiplicity 3,
(0 : 1 : 0)  multiplicity 4,
(1 : 0 : 0)  multiplicity 4.
```

Magma also gives a parametrization from the printed point, but the default
coordinate polynomials have degree `8` and very large coefficients.  So the
current reconstruction is best viewed as an implicit genus-0 family, not yet a
clean human-sized parametrization.

## Simple specialization and inherited `[16]`

The printed member has a quick local `Q`-simplicity certificate.  At `p=7`,
Magma gives

```text
L_7(T) = 49*T^4 + 7*T^3 + 6*T^2 + T + 1,
```

which is irreducible over `Q`.  Thus the reconstructed component is not
contained in the decomposable locus.

The same family also supplies an inherited infinite `[16]` row: if the marked
class `D` has order `32` generically, then `2D` has order `16` generically.
This is not a new moduli construction, but it was missing from the inventory.

## Files

- `code/elkies32_reconstruct_conditions.py` derives the rank-drop equations and
  the component containing the printed point.
- `data/elkies32_reconstruct_conditions.txt` contains the exact equations,
  factorization, and the formula for `p`.
- `code/elkies32_component_genus.m` verifies the geometry of the plane
  component in Magma.
- `data/elkies32_component_genus.txt` contains the genus and singularity output.
- `code/elkies32_simple_certificate.m` checks the printed member for an
  irreducible Frobenius certificate.
- `data/elkies32_simple_certificate.txt` records the `p=7` certificate.

## Conclusion

Yes: we effectively have reconstructed the Elkies `[32]` family, in implicit
form.  Rational points `(z,r)` on the genus-0 component, away from boundary and
singular/discriminant loci, give curves

```text
y^2 = (a*x^3 + b*x^2 + c*x + 1)^2 - a^2*x^5*(x+1)
```

with `D=(0,1)-infinity` satisfying `16D=(r,0)-infinity`; generically this gives
exact order `32`.

The best next step is not another coefficient search.  It is to find a smaller
parametrization of the genus-0 plane component, probably by resolving one of
the rational singular branches at infinity rather than using Magma's bulky
default parametrization.
