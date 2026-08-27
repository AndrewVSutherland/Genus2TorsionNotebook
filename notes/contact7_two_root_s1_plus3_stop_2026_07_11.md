# Contact-7 two-root `s=1` branch plus 3: stop certificate

Date: 2026-07-11.

## Outcome

The removable `s=1` branch is a genuine one-parameter family, but its
3-torsion support cover fails the agreed low-genus stop rule.  A full-degree
modular resolvent certifies that the degree-40 support cover is connected,
and six simple nodal fibers alone imply genus at least `15`.  Independently,
a complete conservative search through parameter height `10000` found no
candidate.

## Exact family

Put

```text
a = 35/8,
b = -(8t^4+24t^3+48t^2+45t+15)/(8(t+1)^3),
h = 1-(7/2)x+a*x^2+b*x^3,
f = (h^2+(x-1)^7)/x^2.
```

Then

```text
f = x*(x+t^2-1)*g_t(x),
h(1) = -t^2*(8t^2+9t+3)/(8(t+1)^3).
```

Thus the roots `0` and `1-t^2`, together with the marked contact-7 class,
give visible `Z/2 x Z/14` away from the displayed degeneracies.

Up to a nonzero rational constant, the discriminant is

```text
(t-1)^2 * t^14 * (8t^2+9t+3)^7
* (16t^3+29t^2+20t+5)^2
* (4t^3+16t^2+12t+3)^2 * K(t) / (t+1)^36,
```

where

```text
K(t)=864t^6+1760t^5+920t^4-249t^3-249t^2+45t+45
```

is irreducible over `Q` and occurs with exponent one.  Its six roots are
therefore six simple nodal fibers.

## Connected support cover and genus stop

For the support of a nonzero 3-torsion class modulo sign, use

```text
w*(1+r1*x+r2*x^2+r3*x^3)^2
    = f + k*(x^2+U*x+V)^3.
```

After saturating the trivial component, the calculation over `F_5(t)` gives

```text
triangular basis degrees in M   39, 39, 40
M-resolvent degree              40
factor degrees                  [40]
L^2-M irreducible               true
```

Prime `5` is a good reduction for the generic equations and the degree does
not drop.  Thus the characteristic-zero degree-40 cover

```text
(J[3]-{0})/{+1,-1} -> P^1_t
```

is connected; the signed degree-80 cover is connected as well.  The exact
transcript is in `data/contact7_two_root_plus3_geometry_s1_f5.txt`.  At a
simple nodal fiber, Picard--Lefschetz monodromy is a
symplectic transvection.  On the 40 nonzero sign classes its cycle shape is

```text
1^13 3^9,
```

with ramification contribution `18`.  The six roots of `K` alone give

```text
2g-2 >= -2*40 + 6*18 = 28,
g >= 15.
```

Additional boundary or infinity ramification can only increase this bound.
On the 80 nonzero vectors a transvection contributes `36`, so the six nodes
also give genus at least `29` for the connected signed cover.  In particular,
neither cover has a genus `0`, `1`, or `2` component.

## Complete height-10000 sieve

The dedicated one-parameter sieve tested every reduced `t=n/d` with

```text
|n| <= 10000,  1 <= d <= 10000,
t != 0,+1,-1.
```

At each prime it discarded a residue only when the model had good reduction
and the prime-to-`p` part of `84` failed to divide `#J(F_p)`.  Denominator,
boundary, and bad-reduction disks were retained.  Using 52 primes from `5`
through `251` gave

```text
rational parameters  121,589,943
exact boundaries                3
tested              121,589,940
survivors                      0
```

In particular, at `p=5` the only good affine residues are `t=2,3`; both
have `#J(F_5)=28`, so any global hit is forced into a `5`-adic boundary
disk.  The full multi-prime sieve nevertheless treats those disks
conservatively.

## Files

```text
code/contact7_two_root_s1_plus3_monodromy.m
code/contact7_two_root_s1_plus3_modular_resolvent.m
code/contact7_two_root_s1_plus3_sieve.cpp
data/contact7_two_root_s1_plus3_h10000.txt
data/contact7_two_root_s1_plus3_h10000_survivors.m
data/contact7_two_root_plus3_geometry_s1_f5.txt
```

Decision: stop this branch.  It is arithmetically negative in a large
one-parameter box and geometrically far above the low-genus threshold.  This
does not prove that the full two-parameter surface has no order-84 point; the
full-surface result is the separate, bounded height-100 search.
