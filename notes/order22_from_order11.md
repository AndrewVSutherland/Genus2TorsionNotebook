# Order 22 from order 11 infinity-torsion families

This is a constructive follow-up to the internet sweep.  Instead of looking for
extra `2`-torsion on the `[23]` family, start from the one-parameter genus-2
families where the difference of the two points at infinity has order `11`.
If such an even sextic also has a finite rational branch point, then it gives a
cyclic point of order `22`.

## Divisor argument

Let

```text
C: y^2 = f(x)
```

be an even sextic with leading coefficient a square, so the two points at
infinity `infinity_+` and `infinity_-` are rational.  Suppose

```text
D_inf = infinity_+ - infinity_-
```

has order `11`.  If `W=(r,0)` is a rational finite Weierstrass point, then

```text
div(x-r) = 2W - infinity_+ - infinity_-.
```

Thus for `E = W - infinity_+`,

```text
2E = infinity_- - infinity_+ = -D_inf.
```

Since `D_inf` has order `11`, the class `E` has order `22`.

## Flynn family

Flynn's order-11 family is

```text
F_t(x) = x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3
         + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2.
```

For a candidate rational root `r`, the equation `F_t(r)=0` is quadratic in
`t`:

```text
(r-1)^2*t^2
+ 2*r*(r+1)*(r^2-r+1)*t
+ r^2*(r^2+r+1)^2 = 0.
```

Its discriminant is

```text
16*r^5.
```

Therefore set `r=s^2`.  For `s != 0,+/-1`, the two branches are

```text
t_eps(s) = (-s^2*(s^2+1)*(s^4-s^2+1) + 2*eps*s^5)/(s^2-1)^2,
eps = +/-1.
```

Then `x=s^2` is a rational branch point, so the Jacobian contains a rational
point of order `22`.

## Daowsud-Schmidt family

The Daowsud-Schmidt order-11 family is

```text
G_u(x) = x^6 - 4*x^5 + 8*(1+u)*x^4 - (10+32*u)*x^3
         + 8*(1+6*u+2*u^2)*x^2
         - 4*(1+6*u+16*u^2)*x + 64*u^2 + 1.
```

For a candidate rational root `r`, `G_u(r)=0` is quadratic in `u` with

```text
a = 16*(r-2)^2,
b = 8*(r-1)*r*(r^2-3*r+3),
c = (r-1)^2*(r^2-r+1)^2,
disc = 256*(r-1)^5.
```

Set `r=1+s^2`.  For `s != 0,+/-1`, the two branches are

```text
u_eps(s) = (-s^2*(s^2+1)*(s^4-s^2+1) + 2*eps*s^5)/(4*(s^2-1)^2),
eps = +/-1.
```

Then `x=1+s^2` is a rational branch point, again giving a rational class of
order `22`.

The two order-11 source families are known to be closely related; I am not
claiming the two resulting order-22 parameterizations are inequivalent moduli
families.  The important point for this search is that they give explicit
infinite `[22]` simple-Jacobian families with small formulas.

## Verification samples

For `s=2`, Magma gives exact torsion `[22]` on all four branch choices:

```text
Flynn eps=-1: t=-36
  y^2 = x^6 + 2*x^5 - 69*x^4 + 2*x^3 + 1297*x^2 - 2664*x + 1296
  TorsionSubgroup(J)(Q) = [22]
  irreducible Lp at p=5: 25*x^4 + 10*x^3 + 6*x^2 + 2*x + 1

Flynn eps=+1: t=-196/9
  primitive model:
  y^2 = 81*x^6 + 162*x^5 - 3285*x^4 + 162*x^3
        + 38497*x^2 - 80360*x + 38416
  TorsionSubgroup(J)(Q) = [22]
  irreducible Lp at p=5: 25*x^4 - 5*x^3 + 2*x^2 - x + 1

Daowsud-Schmidt eps=-1: u=-9
  y^2 = x^6 - 4*x^5 - 64*x^4 + 278*x^3 + 872*x^2 - 4972*x + 5185
  TorsionSubgroup(J)(Q) = [22]
  irreducible Lp at p=5: 25*x^4 + 10*x^3 + 6*x^2 + 2*x + 1

Daowsud-Schmidt eps=+1: u=-49/9
  primitive model:
  y^2 = 81*x^6 - 324*x^5 - 2880*x^4 + 13302*x^3
        + 17896*x^2 - 143404*x + 153745
  TorsionSubgroup(J)(Q) = [22]
  irreducible Lp at p=5: 25*x^4 - 5*x^3 + 2*x^2 - x + 1
```

Sage/Lombardo checks give

```text
geometric_endomorphism_algebra_is_field(B=100) = True
geometric_endomorphism_ring_is_ZZ(B=100) = True
```

for all four samples.  Thus these families are not contained in the decomposable
locus; the generic member is geometrically simple away from a proper exceptional
set.

## Files

```text
code/order22_from_order11_check.m
code/order22_from_order11_geom_check.sage
data/order22_from_order11_check.txt
data/order22_from_order11_geom_check.txt
```
