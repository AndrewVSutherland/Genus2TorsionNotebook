# Local feasibility of the contact-6 route to `[6,12]`

## Scope

This is a bounded go/no-go test for extending the packaged geometrically
simple `[6,6]` example in `notes/contact6_m36_66_example.md`.  The marked
order-6 contact class `D` is halved while the independent cubic-contact class
and the `[1,2,2]` factorisation used by the successful `[6,6]` core are kept.
No rational height search is performed.

The reproducible files are

```text
code/contact6_m36_612_local_feasibility.sage
data/contact6_m36_612_local_feasibility_p5_p11.txt
```

Run them with

```bash
sage code/contact6_m36_612_local_feasibility.sage --primes 5,11 \
  > data/contact6_m36_612_local_feasibility_p5_p11.txt
```

## Intersection and screens

The variables are

```text
(a,b,L,U,v,q,r).
```

The script imposes the three cubic-contact equations `F1=F2=F3=0`.  For the
halving identity it uses `A=b+3`, `C=a+b+2` and eliminates the two variables
which are linear on `A*C != 0`:

```text
p = (A^2 + 2*(A+C-2) - r^2)/(4*A),
s = A*q^2/C - (A+C-2).
```

The remaining equations are the existing `K2=K3=0`, equivalently
`H1=H2=H3=H4=0` on this chart.  Every recovered solution is checked again in
the full polynomial identity

```text
ell^2 - f = -2*A*t*u^2.
```

The pointwise saturation removes

```text
A*C*L*v*(v^3-1)*(U^2-4*v^2)*q = 0,
disc(f) = 0, gcd(x^2+U*x+v^2,f) != 1, gcd(u,f) != 1.
```

It also requires the good-prime factor degrees to be compatible with a
rational factorisation `x*(quadratic)*(quadratic)`: all irreducible factors
must have degree at most two.

The finite-group screen is not based only on `72 | #J(F_p)`.  The script
enumerates all reduced Mumford representatives, recovers the full abelian
invariants from multiplication kernels, and applies the valuation-rank
criterion for embedding `Z/6 x Z/12`.  It then constructs `D`, its half `H`,
the cubic-contact class `E`, and rational 2-torsion classes, and verifies that
some `E+T` together with `H` generates 72 distinct classes.

## Results at 5 and modulo 25

The complete affine enumeration over `F_5` has 476 solutions of the five
intersection equations.  There are

```text
fully open core points:             0
full subgroup-containment passes:   0
explicit Z/6 x Z/12 passes:         0
smooth nonsplit core points:        0
```

Every residue point lies on at least one explicit degeneracy or bad-curve
locus.  This strengthens the earlier `M(2,12)` observation that the broader
`[3,12]` problem is forced to the `p=5` boundary.

Lifts modulo 25 are counted exactly from

```text
F(x_0+5*t)/5 = F(x_0)/5 + Jac(F)(x_0)*t mod 5.
```

There are 42 rank-five intersection points, each with `5^2=25` lifts modulo
25, as well as singular intersection points with larger tangent spaces.
None is a lift of a smooth nonsplit core stratum: all 476 special fibres have
already failed the nonboundary/core screens.  Thus the affine open route has
no `Z_5` point and no useful mod-25 seed.

## Independent prime

At `p=11` the full intersection has 4,972 affine points.  Exactly 56 are
fully open `[1,2,2]`-compatible points.  All 56 pass both the full ambient
subgroup test and the explicit generator test, and all have

```text
J(F_11) invariants = [12,12].
```

They are not nonsplit seeds.  Every one has

```text
P_11(T) = (T^2+11)^2,
P_11^[12](W) = (W-11^6)^4.
```

Hence these reductions are geometrically supersingular split.  Each point is
rank five on the intersection and has exactly `11^2=121` lifts modulo 121.
The algebraic intersection is therefore genuinely nonempty away from 5, but
the visible open stratum at this comparison prime is entirely split.

## Diagnostic

**No-go for a broad rational search on the current affine open cover.**  The
required open locus is absent modulo 5, so any rational specialization in
this construction must enter a 5-adic boundary chart.  The comparison prime
does not supply a nonsplit component to target: its 56 open target points are
all split supersingular.

This is not a global nonexistence theorem for simple `[6,12]`.  A rational
point with negative 5-adic valuations can reduce to the projective boundary,
and after rescaling its generic curve need not equal the singular affine
special fibre.  Continuing this route is justified only after blowing up or
normalising a specific `p=5` boundary component and exhibiting a smooth
generic fibre with independent torsion.  The present data do not justify a
height-box search.
