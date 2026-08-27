# Simultaneous Contact-5/Contact-6 Order-30 Family

This constructs an infinite family of genus-2 curves with a rational point of
order `30` on the Jacobian.

## Contact Setup

Seek an odd quintic `f` with two contact presentations

```text
f = h6^2 - (x-1)^6 = h5^2 - K*x^5,
```

where

```text
h6 = x^3 + A*x^2 + B*x + C,
h5 = e*x^2 + d*x + c.
```

Then `y-h6` gives a rational divisor class `D6=(1,h6(1))-infinity` of order
dividing `6`, and `y-h5` gives `D5=(0,c)-infinity` of order dividing `5`.
If the specialization is smooth and nonboundary, `D5+D6` has order `30`
generically.

The coefficient equations are

```text
C^2 - c^2 = 1,
2*B*C - 2*d*c + 6 = 0,
B^2 + 2*A*C - d^2 - 2*e*c - 15 = 0,
2*A*B - 2*e*d + 2*C + 20 = 0,
A^2 - e^2 + 2*B - 15 = 0,
K = -2*A - 6.
```

Parametrize `C^2-c^2=1` by

```text
C = (u^2+1)/(2u),     c = (u^2-1)/(2u).
```

Use `s=A+e` and `q=A-e`.  Eliminating `q` gives a genus-zero plane curve in
`(u,s)`.  After the substitution `u=t^3`, one rational branch is controlled by
the conic

```text
Y^2 = 5*t^2 - 6*t + 5.
```

Using the rational point `(t,Y)=(5,10)`, set

```text
t = (5*R^2 - 20*R + 19)/(R^2 - 5),
Y = -2*(5*R^2 - 22*R + 25)/(R^2 - 5),
u = t^3,
s = t^5 + t^4 + (5/2)*t^3 + (1/2)*t
    +/- t*(t-1/2)*(t+1)*Y.
```

Then recover

```text
q = (15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2 + 231*u^3
     + 2*u^2*s - 15*u*s^2 + 90*u^2 - 20*u*s + 15*u - 2*s)
    /(u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s - u*s^3 + u^2),

A = (s+q)/2,      e = (s-q)/2,
B = (15-s*q)/2,  d = (B*C+3)/c.
```

Finally

```text
C_R: y^2 = (x^3 + A*x^2 + B*x + C)^2 - (x-1)^6.
```

The excluded parameters are the finite set where a denominator vanishes,
`deg(f)<5`, or `disc(f)=0`.

## Verification

Scripts:

```text
code/contact5_contact6_order30_symbolic.m
code/contact5_contact6_order30_sumdiff.m
code/contact5_contact6_order30_curve_sage.py
code/contact5_contact6_order30_fiber_search.py
code/contact5_contact6_order30_verify_fibers.m
code/contact5_contact6_order30_family.m
```

Main outputs:

```text
data/contact5_contact6_order30_symbolic.txt
data/contact5_contact6_order30_sumdiff.txt
data/contact5_contact6_order30_fibers_h500.txt
data/contact5_contact6_order30_verify_h500.txt
data/contact5_contact6_order30_family_samples.txt
```

The Sage probe confirms that the eliminated core has genus `0`.  The height
`500` fiber search found the rational fibers

```text
u=125,   s=5415
u=125,   s=2715
u=1/125, s=831/3125
u=1/125, s=-69/3125
```

plus one singular parametrization-center point.

For these smooth fibers, Magma verifies

```text
Order(D5) = 5,
Order(D6) = 6,
Order(D5+D6) = 30,
TorsionSubgroup(J)(Q) = [30].
```

The two distinct curves from `R=5/2` are geometrically simple by irreducible
Frobenius certificates:

```text
branch -1: p=19, Lp = 361*T^4 - 2*T^2 + 1,
branch +1: p=11, Lp = 121*T^4 - 2*T^2 + 1.
```

Other tested parameters (`R=0,4,6,-1,1/2`) also gave exact torsion `[30]` and
simple certificates.

## Conclusion

This gives a concrete infinite rational family with generic rational
`30`-torsion.  Because the family has specializations with irreducible
Frobenius quartics, it is not contained in the decomposable locus; hence the
generic member is geometrically simple, with only a finite exceptional set of
non-simple specializations.
