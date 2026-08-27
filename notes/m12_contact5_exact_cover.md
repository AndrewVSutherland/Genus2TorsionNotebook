# Exact contact-5 cover of the split M(12) surface

## Goal

The split-Weierstrass `M(12)` surface carries a rational class `D12` of exact
order 12.  A point-contact class `D5=P-infinity` of exact order 5 on the same
odd quintic gives

```text
D60 = D12+D5
```

of exact order 60.  Since a generic member of `M(12)` has torsion `[12]`, this
route can produce exact cyclic `[60]`; unlike the extra-2 Richelot route, it
does not force a second rational 2-direction.

## Full split-Weierstrass chart

Start from

```text
y^2+(x-r)*(T+1)*y = a*x^2*T*(T+1),
T = a*x^2-x+r.
```

The first implementation used

```text
a = (1-z^2)/(4*(r+1)),
w = 2*(r+1)/(1-z),
```

so `w` is a rational root of `T+1`.  Moving `w` to infinity gives an odd
quintic and a directly verified order-12 divisor.  The finite point-count
funnel in `code/m12_full_surface_plus5_order60_search.m` is locally open at
`p=7`, unlike the old one-parameter line `a=(1-r)/4`.

Its exhaustive height-15 run through primes at most 31 checked 82,369
`(r,z)` pairs.  There were 638 finite-mask survivors, 68 smooth survivors,
34 distinct curves after the `z`-sign duplication, and 34 exact torsion
computations.  Every exact group was `[12]`; there was no order-60 hit.

## Exact contact-5 covariants

For an odd quintic `f` and a proposed contact abscissa `u`, put

```text
Ai = f^(i)(u)/i!,
Q2 = 4*A0*A2-A1^2.
```

Writing a quadratic contact polynomial locally as

```text
h = c+d*(x-u)+e*(x-u)^2,
```

the first three Taylor equations determine

```text
c^2=A0,
d=A1/(2*c),
e=(A2-d^2)/(2*c).
```

The remaining two equations are exactly

```text
E3 = 8*A0^2*A3-A1*Q2 = 0,
E4 = 64*A0^3*A4-Q2^2 = 0,
c^2=A0 != 0.
```

These are necessary and sufficient: reconstruction then verifies

```text
f-h^2 = LeadingCoefficient(f)*(x-u)^5.
```

Thus `E3=E4=0` plus the square condition defines the exact one-dimensional
contact-5 cover in `(r,z,u)`.  The implementation is
`code/m12_contact5_exact_cover_search.m`.  It includes a planted positive
control at `t=2` in the standard contact-5 family; the script must recover
`u=0` and the exact contact identity before it begins any search.  Repeated
rational roots of the gcd are retained, since multiplicity means ramification,
not invalidity.

## Exact finite masks and height-50 result

The class-specific finite mask asks for a finite `u` satisfying `E3=E4=0`
with `f(u)` a nonzero square.  It is much sharper than the necessary condition
`5 | #J(F_p)`, but it remains nonempty on good fibers:

```text
p     allowed good base pairs / good base pairs
7          1 / 18
11         2 / 70
13         1 / 110
17         3 / 198
23         5 / 404
29         8 / 680
```

Consequently there is no good-open local obstruction at these primes.  The
corrected height-50 run used 3,095 reduced rational values for each of `r,z`,
hence 9,579,025 pairs.  Its counters were

```text
rigorous finite-mask survivors    19,128
smooth exact gcd/square tests     12,942
exact rational cover points            0
order-60 hits                          0
exact cyclic [60] hits                 0
```

This is evidence of a thin global cover, not impossibility.

## Smaller normalized root chart

For global geometry, `(r,z)` is unnecessarily large.  Choose the rational
root `w` directly and set

```text
r = w-1-a*w^2,
z = 2*a*w-1.
```

Equivalently, put `b=a*w`, move `w` to infinity, then scale the odd coordinate
by `xi=w*X` and the y-coordinate by `eta=w^2*Y`.  The order-12 point now has
fixed abscissa `xi=-1`.  With

```text
L = b+(2*b-1)*xi,
A = xi+w*(1+b*xi),
```

the odd quintic is the compact polynomial

```text
F = L*(L*A^2+4*b*(1+xi)^2*(w*L-xi^2)).
```

The open chart has `b*w*(1-b) != 0`; degree five additionally needs
`2*b-1 != 0`.  Its leading coefficient is `-4*b*(2*b-1)`, independent of
`w`, and the visible order-12 point is

```text
(xi,eta)=(-1,(1-b)*(w*(1-b)-1)).
```

This `(b,w)` chart is birational to the original split surface and is the
preferred chart for factoring and saturating the exact equations.  A bounded
height search here is also genuinely new: small `(b,w)` can map to `(r,z)` of
height far above 50.

## Next step

Build `E3,E4` from the compact `F`, clear denominators (there are none in this
normalization), factor out the loci

```text
b*w*(1-b)*(2*b-1),
u+1,                    // contact point collides with D12
L(u),                   // contact point is Weierstrass
discriminant(F),
```

and normalize the remaining common component(s).  Search rational points on
that normalized cover, rather than enlarging the `(r,z)` height box.
