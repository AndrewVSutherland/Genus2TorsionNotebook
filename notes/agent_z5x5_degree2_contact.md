# Agent notes: one contact-5 plus degree-2 contact/Mumford route

Files:

```text
code/agent_z5x5_degree2_contact_symbolic.m
code/agent_z5x5_degree2_contact_search.m
```

I pursued route (A): keep one rational 5-contact at `x=0`, and try to make a
second rational 5-torsion class from a degree-2 divisor.

## Normalized degree-2 contact subroute

Let `C: y^2=f(x)` be an odd quintic.  A rational contact-5 class at `x=0` is

```text
f = h^2 - K*x^5,     deg(h) <= 2.
```

For a degree-2 class `D=[U,V]`, with

```text
U = x^2 + s*x + t,
```

a restrictive but very concrete contact condition is that, after scaling `y`,
there is a monic quintic `H` such that

```text
f = H^2 - U^5 = h^2 - K*x^5.
```

Then `div(y-H)=5D-10*infinity`, so `5D=0`.  Canceling degrees `10` through `6`
in `H^2-U^5` forces

```text
H = x^5
  + (5/2)*s*x^4
  + ((15/8)*s^2 + (5/2)*t)*x^3
  + ((5/16)*s^3 + (15/4)*s*t)*x^2
  + (-(5/128)*s^4 + (15/16)*s^2*t + (15/8)*t^2)*x
  + m.
```

Write

```text
F = H^2 - U^5 = q5*x^5 + q4*x^4 + q3*x^3 + q2*x^2 + q1*x + q0.
```

On the open chart `q4 != 0`, the quartic tail is a square iff

```text
C1 = 8*q4^2*q1 - 4*q4*q3*q2 + q3^3 = 0,
C0 = 64*q4^3*q0 - (4*q4*q2 - q3^2)^2 = 0.
```

This is not the previous fixed quartic obstruction.  Here the equations live
in `(s,t,m)`, and Magma gives

```text
C1 = (s^2 - 4*t)
   * (s^5 - (40/3)*s^3*t + 80*s*t^2 - (256/3)*m)
   * (degree-14 factor),

gcd(C1,C0) = 1.
```

Eliminating `m` gives

```text
(s^2 - 4*t)^39
*(s^2 - (4/5)*t)^8
*(s^4 - 8*s^2*t + (16/5)*t^2)
*(degree-16 factor).
```

For `s != 0`, the last two dehomogenize to irreducible polynomials in
`Z=t/s^2`, so they give no rational ratio.  The rational resultant branches
are therefore:

```text
t = s^2/4       boundary: U has a double root,
t = 5*s^2/4     only common m is m=(45/32)*s^5.
```

The second branch is also a false lead for rational contact: it lands on the
`q4=0` boundary and the remaining tail is not a rational square.  The search
script checks this directly:

```text
s=1  t=5/4  m=45/32   model=quartic tail not square
s=-1 t=5/4  m=-45/32  model=quartic tail not square
s=2  t=5    m=45      model=quartic tail not square
s=-2 t=5    m=-45     model=quartic tail not square
```

The double-root branch is not a genuine degree-2 class:

```text
s=1  t=1/4  m=1/32   model=degree(f)!=5
s=1  t=1/4  m=3/64   model=quartic tail not square
```

## Smoke tests

Command:

```sh
magma -b prime_bound:=43 rational_height:=3 \
  full_field_search:=true \
  code/agent_z5x5_degree2_contact_search.m
```

Small affine rational smoke:

```text
integer s,t,m in [-3,3]:
tested=343 square_smooth=0 independent_hits=0
```

Finite-field smoke does produce honest independent pairs, so the construction
itself is not algebraically empty.  Examples:

```text
p=11, s=1, t=2, m=6
f = 10*x^5 + 5*x^4 + 4*x^3 + 9*x^2 + 9*x + 4
h = 4*x^2 + 6*x + 9
U = x^2 + x + 2
V = 10*x + 10
#J(F_11)=175
```

```text
p=31, s=1, t=3, m=3
f = 30*x^5 + 2*x^4 + 19*x^3 + 27*x^2 + 4*x + 14
h = 8*x^2 + 7*x + 18
U = x^2 + x + 3
V = 21*x + 10
#J(F_31)=1000
```

```text
p=41, s=1, t=5, m=16
f = 3*x^5 + 25*x^4 + 32*x^3 + 15*x^2 + 10*x + 1
h = 36*x^2 + 5*x + 1
U = x^2 + x + 5
V = 28*x + 36
#J(F_41)=1500
```

In each finite-field hit, Magma verifies `D0 != 0`, `D2 != 0`,
`5*D0=5*D2=0`, and no nontrivial relation `a*D0+b*D2=0` for
`a,b in F_5`.

## Verdict

The pure `y-H` degree-2 contact subroute is not viable over `Q`: the exact
elimination leaves only a double-root boundary branch and a `q4=0` branch whose
tail is not a rational square.  This is a different obstruction from the
irreducible quartic `X^5-(X-1)^5`; finite fields such as `F_11`, `F_31`, and
`F_41` have genuine independent pairs in the same equations.

The viable continuation of route (A) is the full degree-2 Mumford norm problem,
not the normalized `B=1` contact subcase.  The next computational object should
be the system

```text
f = h^2 - K*x^5,
A(x)^2 - B(x)^2*f = Lambda*U(x)^5,
U = x^2 + s*x + t,
deg(A) <= 5, deg(B) <= 2,
```

with the congruence `A + B*V == 0 mod U` used to recover the Mumford
representative `[U,V]`.  Equivalently, compute `[5][U,V]=0` by Cantor on the
one-contact family and then add modular/local filters.  This larger system is
where a second rational `5`-class can still hide; the `B=1` contact slice cannot
contain it.
# Degree-2 contact route to `[5,5]`

This is a bounded theory and finite-field study.  It replaces the second
rational point-contact condition in
`notes/agent_z5x5_contact5_contact5.md`; no rational height search or large
elimination was run.

Task-scoped files:

```text
code/agent_z5x5_degree2_contact_probe.sage
data/agent_z5x5_degree2_contact_probe_p11_p19.txt
```

## Correct degree-2 order-5 identity

Work over a field of characteristic different from `2,5`.  Keep the marked
point-contact family

```text
h0 = 1 + a*x + b*x^2,
f  = h0^2 - k*x^5.
```

For `k != 0` and squarefree `f`, the class
`D0=[(0,1)-infinity]` has exact order 5, since

```text
div(y-h0) = 5*((0,1)-infinity).
```

Let a second reduced class have Mumford polynomial

```text
q = x^2 + u*x + v.
```

The naive identity `H^2-f=q^5` is not the general order-5 condition.  A
function with pole at most 10 at infinity has the form `H-R*y`, because
`L(10*infinity)` also contains `x*y` and `x^2*y`.  Consequently the
correct normalized polynomial Pell/contact identity is

```text
R = r0 + r1*x + r2*x^2,
H = x^5 + s4*x^4 + s3*x^3 + s2*x^2 + s1*x + s0,
H^2 - f*R^2 = q^5.                                      (1)
```

The normalization making `H` monic loses no open degree-2 order-5 class:
the principal function for `5D` has pole order 10, whereas `R*y` has pole
order at most 9, so the coefficient of `x^5` in `H` is nonzero.

Conversely, suppose (1) holds and

```text
disc(q) != 0,  gcd(q,f)=1,  gcd(q,R)=1.
```

Then `w = H/R mod q` is well defined and satisfies `w^2=f mod q`.
The norm of `H-R*y` is exactly `q^5`, and the coprimality conditions keep
the conjugate factor nonzero over the support of `q`.  Therefore

```text
div(H-R*y) = 5*D,   D=(q,w),
```

up to replacing `D` by `-D` according to the sign convention.  Since a
reduced degree-2 class is nonzero, it has exact order 5.

## Coefficient equations and dimension

The coefficient of `x^10` in (1) cancels identically.  The ten equations
`E_i=[x^i](H^2-f*R^2-q^5)=0`, for `0 <= i <= 9`, are:

```text
E0 = -r0^2 + s0^2 - v^5
E1 = -2*a*r0^2 - 2*r0*r1 + 2*s0*s1 - 5*u*v^4
E2 = -a^2*r0^2 - 4*a*r0*r1 - 2*b*r0^2 - 2*r0*r2 - r1^2
     + 2*s0*s2 + s1^2 - 10*u^2*v^3 - 5*v^4
E3 = -2*a^2*r0*r1 - 2*a*b*r0^2 - 4*a*r0*r2 - 2*a*r1^2
     - 4*b*r0*r1 - 2*r1*r2 + 2*s0*s3 + 2*s1*s2
     - 10*u^3*v^2 - 20*u*v^3
E4 = -2*a^2*r0*r2 - a^2*r1^2 - 4*a*b*r0*r1 - 4*a*r1*r2
     - b^2*r0^2 - 4*b*r0*r2 - 2*b*r1^2 - r2^2
     + 2*s0*s4 + 2*s1*s3 + s2^2
     - 5*u^4*v - 30*u^2*v^2 - 10*v^3
E5 = -2*a^2*r1*r2 - 4*a*b*r0*r2 - 2*a*b*r1^2 - 2*a*r2^2
     - 2*b^2*r0*r1 - 4*b*r1*r2 + k*r0^2
     + 2*s0 + 2*s1*s4 + 2*s2*s3
     - u^5 - 20*u^3*v - 30*u*v^2
E6 = -a^2*r2^2 - 4*a*b*r1*r2 - 2*b^2*r0*r2 - b^2*r1^2
     - 2*b*r2^2 + 2*k*r0*r1 + 2*s1 + 2*s2*s4 + s3^2
     - 5*u^4 - 30*u^2*v - 10*v^2
E7 = -2*a*b*r2^2 - 2*b^2*r1*r2 + 2*k*r0*r2 + k*r1^2
     + 2*s2 + 2*s3*s4 - 10*u^3 - 20*u*v
E8 = -b^2*r2^2 + 2*k*r1*r2 + 2*s3 + s4^2 - 10*u^2 - 5*v
E9 = k*r2^2 + 2*s4 - 5*u
```

There are 13 variables

```text
a,b,k,u,v,r0,r1,r2,s0,s1,s2,s3,s4
```

and ten equations.  The expected dimension is therefore 3.  More
meaningfully, the ten variables after `a,b,k` should form a finite
level-structure cover of the point-contact base.  The finite tests below
verify this at two open points: the full Jacobian has rank 10 and its
10-by-10 block in the fiber variables also has rank 10.

## Open, boundary, and dependence loci

The open locus used here is

```text
char != 2,5,
k * disc(f) * disc(q) * Res(q,f) * Res(q,R) != 0.        (2)
```

The excluded factors respectively detect degree drop of `f`, a singular
curve, repeated support for `D`, support at a Weierstrass point, and failure
to recover `w=H/R mod q`.  Reducibility of `q` is not a boundary: split
`q` gives two rational support points, while irreducible `q` gives a
quadratic-conjugate pair whose Galois-invariant Mumford class descends.

Condition `disc(q) != 0` also proves independence from `D0` in this chart.
The only nonzero multiples of `D0` having reduced degree 2 are
`+/-2*D0`, whose unique Mumford polynomial is `x^2`; the classes
`+/-D0` have degree 1.  Thus a squarefree quadratic `q` cannot represent a
multiple of `D0`.

## Finite-field component tests

The script enumerates only the two finite fields `F_11` and `F_19`.  It
constructs the Jacobian classes with Sage, recovers (1) by a bounded scan of
the `p^3` possible quadratics `R`, and evaluates the equation Jacobian.

At `p=11` it found

```text
f = 7*x^5 + 1
q = x^2 + x + 7,              w = x + 3
H = x^5 + 2*x^4 + 4*x^3 + 4*x^2 + x + 2
R = 5*x^2 + 7
#J(F_11) = 125
full equation rank = 10, tangent dimension = 3
fiber rank = 10, fiber determinant = 10
```

At `p=19` it found

```text
(a,b,k) = (1,3,16)
f = 3*x^5 + 9*x^4 + 6*x^3 + 7*x^2 + 2*x + 1
q = x^2 + x + 3,              w = 4*x + 13
H = x^5 + 4*x^4 + 9*x^3 + 10*x^2 + 11*x + 18
R = x^2 + 9
#J(F_19) = 475
full equation rank = 10, tangent dimension = 3
fiber rank = 10, fiber determinant = 11
```

In both cases:

- (1) holds exactly over the residue field;
- all factors in (2) are nonzero;
- `q` is irreducible, so the support is a genuine conjugate quadratic pair;
- Sage verifies `5D0=5D=0` and no relation `D=c*D0`, `c in F_5`;
- the 12th-power transform of the Frobenius polynomial is irreducible.

The last item is the standard finite-field certificate that the reduction is
geometrically simple.  Hence these smooth expected-dimensional points are not
confined to an evident decomposable locus.

## What is proved and what is not

Proved:

- (1) is the correct general degree-2 Mumford/contact condition;
- its coefficient scheme has smooth open points of the expected dimension 3
  at two good primes;
- at both points the cover is locally finite over `a,b,k`;
- the two 5-torsion classes are independent and the reductions are
  geometrically simple.

Not proved:

- existence of a rational point on this cover;
- exact rational torsion `[5,5]` on any characteristic-zero specialization;
- generic geometric simplicity in characteristic zero.

The next bounded computation should Hensel-lift one or more of these
fiber-etale points through `p^2,p^3`, combine compatible lifts at two primes,
and try rational reconstruction in `a,b,k`.  Exact torsion and geometric
simplicity should be checked only after reconstruction.  A global Groebner
elimination or blind rational-height scan is not justified yet.
