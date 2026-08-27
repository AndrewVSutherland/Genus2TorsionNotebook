# How the `[2,2,20]` example was found

This is a short reconstruction of the search path.  The `[2,2,20]` example is
fully documented in the current workspace. 

## The geometrically simple `[2,2,20]` example

We started from the quintic-contact 5-torsion construction

```text
h(x) = 1 + t*x + ((t^2 - 1)/2)*x^2,
f(x) = h(x)^2 - ((t + 1)^4/4)*x^5.
```

The contact identity gives a rational 5-torsion class.  In this one-parameter
subfamily there is also an explicit divisor class

```text
H = [x^2 + 2*x/(t+1), (t+2)*x + 1]
```

with

```text
2H = [x-1, 0].
```

Thus a smooth specialization has a rational point of order `20`.

To get extra rational 2-torsion, we factored the residual quartic

```text
f(x)/(x-1).
```

The key simplification was to put

```text
u = t + 1,        y = u*x.
```

Then the residual quartic is equivalent to

```text
u*(y^4 + 4*y^3 + 8*y^2 + 8*y + 4) - 4*y*(y+1)^2.
```

The first extra-2 locus, where the residual quartic has a rational linear
factor, is parametrized by

```text
t = -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4).
```

Searching this parametrized locus found the special value

```text
z = -1/7,
t = -8233/7225,
b = 7790832/52200625.
```

At this point the residual quartic has factor type `1+1+2`:

```text
f/(x-1) = (x - 7225/1296)*(x - 7225/7056)
          *(x^2 - 917575/4032*x + 52200625/28224).
```

Magma computes

```text
J(Q)_tors = [2,2,20].
```

A reduced model for the curve is

```text
y^2 + (x^2 + x)*y =
    -391671*x^6 + 1894851*x^5 + 6846924*x^4
    -15133525*x^3 + 3904068*x^2 + 2625336*x + 254016.
```

This example is geometrically simple.  Sage/Lombardo's endomorphism test gives

```text
geometric_endomorphism_algebra_is_field(B=100) = True
geometric_endomorphism_ring_is_ZZ(B=100)       = True.
```

Equivalently, the good-prime certificate at `p=71` has Frobenius polynomial

```text
X^4 + 2*X^3 + 14*X^2 + 142*X + 5041,
```

whose 12th-power transform is irreducible.


