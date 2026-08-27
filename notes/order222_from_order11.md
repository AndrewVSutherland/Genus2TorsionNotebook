# Attempted `[2,22]` from order-11 infinity torsion

This note records the natural attempt to upgrade the explicit `[22]` family in
`notes/order22_from_order11.md` to `[2,22]`.

## Starting point

Use Flynn's order-11 family

```text
F_t(x) = x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3
         + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2.
```

For

```text
t_eps(s) = -s^2*(s^3-eps)^2/(s^2-1)^2,  eps = +/-1,
```

the point `x=s^2` is a rational branch point, and the divisor
`(s^2,0)-infinity_+` has order `22` on the Jacobian, away from the usual
degenerate values `s=0,+/-1`.

After substitution the sextic factors over `Q(s)` as

```text
F_t_eps(x) = (x-s^2) * q_eps(s,x) * c_eps(s,x)
```

up to harmless square/scalar denominator factors.  Generically this is a
`1+2+3` factorization, giving only one rational 2-torsion class.  An
independent extra rational 2-torsion class, hence torsion containing `[2,22]`,
would require one of the following:

- the residual quadratic `q_eps` splits over `Q`;
- the residual cubic `c_eps` has a rational root.

These are exactly the extra rational even subsets of branch points available
from a `1+2+3` factorization.

## Residual quadratic

For `eps=-1`,

```text
q_-(s,x) =
  -(s-1)^2*x^2
  +(s-1)*(s^3-s^2+s+1)*x
  -(s^2-s+1)^2
```

with discriminant

```text
Delta_- = (s-1)^5*(s^3+s^2-s+3).
```

Away from the boundary `s=1`, this is a square exactly when

```text
Y^2 = (s-1)*(s^3+s^2-s+3)
    = s^4 - 2*s^2 + 4*s - 3.
```

Magma converts this genus-one curve to

```text
y^2 + 8*y = x^3 + 4*x^2 + 16*x,
rank_bounds 0 0, torsion [3].
```

Its rational points are just the two points at infinity and the boundary point
`(s,Y)=(1,0)`.

For `eps=+1`,

```text
q_+(s,x) =
  -(s+1)^2*x^2
  +(s+1)*(s^3+s^2+s-1)*x
  -(s^2+s+1)^2
```

with discriminant

```text
Delta_+ = (s+1)^5*(s^3-s^2-s-3).
```

Away from `s=-1`, this is a square exactly when

```text
Y^2 = (s+1)*(s^3-s^2-s-3)
    = s^4 - 2*s^2 - 4*s - 3.
```

Magma converts this to

```text
y^2 - 8*y = x^3 + 4*x^2 + 16*x,
rank_bounds 0 0, torsion [3].
```

Again the only rational points are the two points at infinity and the boundary
point `(s,Y)=(-1,0)`.

Conclusion: the quadratic residual never splits at a nonboundary rational
parameter.

## Residual cubic

The `eps=-1` residual cubic is

```text
c_-(s,x) =
  s^6*x - s^6 - 2*s^5*x + 2*s^4*x^2 + 2*s^5 + 3*s^4*x
  - 4*s^3*x^2 + s^2*x^3 - 3*s^4 - 2*s^3*x + 5*s^2*x^2
  - 2*s*x^3 + 2*s^3 + 2*s^2*x - 4*s*x^2 + x^3
  - s^2 - 2*s*x + x^2 + x.
```

The other branch is `c_+(s,x)=c_-(-s,x)`.

Sage computes that the projective closure has geometric genus `3` for both
signs.  Thus this condition is not a rational or elliptic positive-dimensional
source of parameters.  A bounded search for rational points with `s` of height
at most `500` found no nonboundary instances where the cubic has a rational
root.

## Current conclusion on this chart

This natural construction does not produce an infinite `[2,22]` family.  More
strongly, inside the explicit `[22]` family obtained from order-11
infinity-torsion:

- the only low-genus extra-2 component is the residual-quadratic splitting
  curve, and it has rank `0` with only boundary rational points;
- the remaining extra-2 condition is a genus-3 curve, with no rational hits up
  to height `500`.

So `[2,22]` needs a different order-11 component or a genuinely different
construction, not this particular one-parameter `[22]` branch-point upgrade.
This is a statement about this chart, **not** evidence that a generic
geometrically simple realization is unlikely.  In particular, increasing the
height bound on the residual cubic is not the best next use of compute time.

## Strategic reset (2026-07-24)

The project table points to `19044.h.2` as a known `[2,22]` example, but that
example lies on an RM locus.  The right goal is therefore not to manufacture
the group for the first time; it is to escape the extra-endomorphism locus
while retaining the level structure.

There is an important dimension heuristic behind the reset.  Choosing a
nonzero rational 2-torsion point and a point of order 11 is finite level data
on the three-dimensional moduli space of principally polarized abelian
surfaces.  It should not, over an algebraic closure, force a Humbert surface.
The Flynn calculation became zero-dimensional only because it first restricted
to a very special one-parameter order-11 family and then imposed another
factorization condition.  It should not be treated as the ambient geometry.

### Lane A: broad three-root search (run first)

The new script `code/order222_three_root_sieve.m` searches the normalized chart

```text
y^2 = x(x-1)(x-r)(x^3+a*x^2+b*x+c),
```

with irreducible residual cubic.  The `1+1+1+3` branch-orbit pattern supplies
exactly two independent rational 2-torsion classes.  (The tempting
`1+1+4` chart supplies only one and therefore cannot contain `[2,22]`.)
Before asking Magma for an exact
rational torsion subgroup, the script computes a running gcd of `#J(F_p)` at
good primes and rejects a curve as soon as the gcd is not divisible by `44`.
Every exact `[2,22]` hit is printed together with an irreducible Frobenius
polynomial when one is found.

Run successive boxes, retaining the complete survivor logs:

```text
magma -b B:=6  exact_cap:=25  code/order222_three_root_sieve.m
magma -b B:=12 exact_cap:=100 code/order222_three_root_sieve.m
magma -b B:=24 exact_cap:=500 code/order222_three_root_sieve.m
```

For the larger boxes, shard on `r` before parallelizing.  Do not merely raise
`exact_cap`: first add more reduction primes to collapse false positives.  A
hit still needs the usual geometric-endomorphism-ring check; an irreducible
quartic Frobenius polynomial proves geometric irreducibility but does not by
itself prove `End(J_bar)=Z`.

### Lane B: deform the RM seed in level coordinates (highest-value lane)

Export from the database record `19044.h.2` all of the following, rather than
only its curve equation:

1. a sextic model and its two independent rational 2-torsion classes;
2. a Mumford representative `(u,v)` for a point of order 11 (or for the
   order-22 generator);
3. the RM discriminant and the certificate used to identify RM.

Put the curve in the same three-root chart as Lane A.  Introduce variables for
`r,a,b,c,u_1,u_0,v_1,v_0`, impose `u | (f-v^2)`, and impose `[11](u,v)=0`
using denominator-cleared Cantor arithmetic.  Then:

1. compute the Jacobian matrix at the RM seed modulo several good primes;
2. find tangent vectors not tangent to the seed's Humbert equation;
3. Hensel-lift those transverse residue directions to precision `p^k`;
4. perform simultaneous rational reconstruction with increasing height;
5. verify every reconstructed point exactly, never from its p-adic
   approximation.

This is the most direct way to turn the existing RM point from a dead end into
an anchor.  A transverse tangent proves the level scheme is not locally trapped
in that Humbert surface.  Failure of transversality at one seed means only that
the chosen point or chart is singular; repeat at a quadratic twist or an
isomorphic three-root normalization before drawing a global conclusion.

### Lane C: mine the database neighborhood

Search the database export beyond the single displayed record.  Filter first
by rational torsion divisible by `[2,22]`, then test every isomorphism class for
geometric endomorphism ring.  Also inspect curves with conductor near `19044`
and the same `1+1+1+3` branch-orbit pattern: small twists and neighboring models
can preserve useful mod-`p` level signatures even when rational torsion changes.
This lane is cheap and may solve the problem immediately; more importantly it
can provide several seeds for Lane B instead of overfitting one RM point.

### Stop/go rules for Fable 5

- **Stop** enlarging the old height-500 residual-cubic search until its
  genus-3 rational points are attacked structurally (2-cover descent plus
  Chabauty), because another rectangular scan adds little information.
- **Go** on Lane A until the reduction-gcd survivor rate is measured at two
  box sizes.  If it stays too high, add primes; if it is zero, enlarge the box.
- **Go** on Lane B as soon as the RM seed's torsion generator is exported.  The
  first deliverable is a modular tangent-rank/transversality report, not a huge
  elimination polynomial.
- **Accept** a realization only after exact `[2,22]`, smoothness, an absolutely
  irreducible/simple certificate, and an `End(J_bar)=Z` (or at minimum a
  certificate excluding the known RM field) are all checked independently.

## Files

```text
code/order222_from_order11_conditions.sage
code/order222_three_root_sieve.m
data/order222_from_order11_search_h500.txt
```
