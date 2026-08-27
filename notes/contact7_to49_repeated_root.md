# Repeated-root `2P-K` route toward cyclic 49

Set `U=(x-p)^2` in the contact-7 norm equation.  On an odd-degree
genus-2 model, this is the Mumford class

```
Q = 2P-K = 2(P-infinity).
```

The required relation `7Q=+/-D7` is encoded by

```
A^2-B^2*f = -(x-1)*(x-p)^14,
deg(A) <= 7,  B monic of degree 5.
```

## Safe elimination and dimension

There are 16 variables

```
alpha,beta,p,a0,...,a7,b0,...,b4
```

and 15 coefficient equations.  A globally triangular order avoids
the expression swell seen in the generic division cover:

1. Degrees 14 through 10 solve `b4,...,b0`; every pivot is the
   constant `-2`.
2. After those substitutions, degrees 9, 8, 7 solve `a2,a1,a0` on
   the chart `a7 != 0`, with pivots `2*a7`.
3. Seven residual equations remain in the eight variables
   `alpha,beta,p,a3,...,a7`.

Thus the expected open repeated-root locus is one-dimensional.  This
does not prove irreducibility.  The chart `a7=0` can contain isolated
or exceptional pieces and was not eliminated by division.  The
specialization `p=1` has no open solution: there `P` is one of the
marked points and `7(2P-K)=0`, not `+/-D7`.  Weierstrass points and
solutions with `B(p)=0` are likewise outside the cyclic-49 open chart.

## Local obstruction

At every good prime a solution must satisfy the marked-point equation

```
14*[P-infinity] = +/-D7  in J(F_p).
```

At `p=3` none of the five smooth marked contact fibers has such a
point.  Consequently any rational repeated-root solution must meet a
bad or nonintegral 3-adic boundary.  At `p=5` the condition is also
very restrictive: the only points occur on `(alpha,beta)=(2,1)` and
have `(x(P),y(P))=(4,+/-1)`.

## Rational search

For `p != 0,1`, rational points were parameterized without a square
test by

```
c = (p-1)^7,
Y = (r+c/r)/2,
H = (c/r-r)/2,
beta = (H-1+(7/2)*p-alpha*p^2)/p^3,
y(P) = Y/p.
```

The separate `p=0` chart uses
`alpha=y(P)^2/2+35/8`.  Exact local masks through 43 screened all
rational parameters of height at most 8: 635,970 generic cases and
7,569 cases on the `p=0` chart.  Nineteen candidates survived the
finite masks, but exact rational checks showed that every one was a
singular curve.  No cyclic-49 example was found.  Peak memory was
22.5 MB.

This is a bounded negative search, not a proof that the expected curve
has no rational open points.  A realistic continuation would impose
the forced 3-adic bad-reduction charts before attempting lifting,
rather than enlarge the naive height box.

Implementation: `code/contact7_to49_repeated_root.m`.
Numerical summary: `data/contact7_to49_repeated_root_summary.txt`.
