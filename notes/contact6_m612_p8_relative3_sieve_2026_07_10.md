# A degree-12 support quotient, degree-24 signed cover, and parameter sieve on `P8` (2026-07-10)

Updated 2026-07-11 with the exact characteristic-zero cover and local disks.

## Outcome

The missing independent rational 3-direction on the rational endpoint-R3
`P8` family does not require the full degree-40 support quotient.
Symplectic pairing isolates a degree-12 relative **support** factor, and modular
elimination verifies the expected splitting

```text
40 = 1 + 12 + 27
```

after pullback to the `P8` parameter at `p=7,11,13`.

This factor is no longer only modular.  An exact characteristic-zero
reconstruction gives an irreducible polynomial `Phi12(v)` over `Q(e)`, exact
degree-`11` recovery maps for `M` and `U`, and the signed primitive polynomial

```text
Phi24(L) = Res_v(Phi12(v),L^2-M(v)).
```

Here `M` is nonsquare in the degree-`12` support field, `Phi24` is
irreducible of degree `24`, and `L` itself is primitive.  Coefficientwise
base change by the rational `P8` function `e=e(u)` gives the exact
characteristic-zero P8 cover; its good reduction at `7` is again irreducible
of degree `24`.

A rigorous finite-field parameter sieve through `p=71` was then run on every
reduced finite parameter

```text
tau=n/d,  |n|<=1000,  1<=d<=1000.
```

Among 1,216,767 parameters, the only mask survivors were

```text
tau=-3, 2,
```

and both have `e=0`, so they are boundary values rather than members of the
open `P8` family.  Thus there is no nonboundary extra-3 specialization in
this height box.

The height result is rigorous and bounded; it is not a proof that the exact
degree-24 cover has no rational point globally.  Moreover, subsequent local
analysis found smooth signed orthogonal branches in some forced disks at both
`p=7` and `p=17`, so neither prime supplies a global local obstruction.

## Why the relevant support quotient has degree 12

The marked order-6 class on the source supplies a nonzero rational
`P in J[3]`.  Nonzero 3-torsion classes modulo sign give

```text
(3^4-1)/2 = 40
```

cubic-contact supports.  One is the built-in support `{P,-P}`.

If a second class `Q` is rational, its Weil pairing with `P` is a rational
cube root of unity.  Since `mu_3(Q)={1}`, one must have

```text
e_3(P,Q)=1,
```

so `Q` lies in the three-dimensional symplectic hyperplane `P^perp`.
Removing the line spanned by `P` leaves 24 vectors, or 12 supports modulo
sign.  The nonorthogonal complement has 54 vectors, or 27 supports.  Hence
the relative support cover naturally decomposes into orbit sizes

```text
1 (built in), 12 (orthogonal and relevant), 27 (nonorthogonal).
```

The support forgets the sign of `Q`.  In the cubic-contact equations this is
the quotient `M=L^2`.  Recovering a signed Mumford class requires the square
cover `L^2=M`; the exact calculation below shows that this doubles the
orthogonal component from degree 12 to degree 24.

Richelot isogeny has degree four, prime to three, so it preserves the full
rational 3-primary module.  The same degree-24 signed condition may therefore
be tested on the simpler source curve instead of the dual.

## Exact and modular resolvent certificates

Over `Q(e)`, saturation by `M` gives the exact degree-`40` support quotient.
The universal length-`1` factor is the repeated support

```text
(M,U,v)=(1,-2,1),    U^2-4*v^2=0.
```

The exact degree-`12` factor `Phi12` is irreducible.  Its recovered maps
`M(v),U(v)` satisfy all three residual contact equations identically modulo
`Phi12`, and gcd checks remove `M*v*(U^2-4*v^2)`.  An exact square test gives
`M` nonsquare in `Q(e)[v]/(Phi12)`.  The resultant `Phi24(L)` above is
irreducible, and a linear gcd recovers `v` from `L`, hence also `U,N,R`.
Thus this is an exact signed cover, not a component inferred from modular
orbit sizes.  Full formulas and recovery maps are in
`notes/contact6_m612_p8_relative3_exact_2026_07_10.md` and
`code/contact6_m612_relative3_exact_reconstruct.m`.

For the `P8` normalization parameter `u`, substitute

```text
t = 4*(u^2+u-6)/(u^2+6),
e = -(25/3)*t^2/(t^4-25*t^2+1250/3),
a = 1/e,
b = 0
```

into the cubic-contact equations in `(M,U,v)`, where `M=L^2`.  Saturating
by `M` gives a zero-dimensional affine algebra of dimension 40 over
`F_p(u)`.  A lexicographic change of order produces a squarefree degree-40
resolvent in `v`.

At each of `p=7,11,13`, its factor degrees are exactly

```text
1, 12, 27,
```

all with multiplicity one.  At `p=7`, the linear factor is `v-1`, the
known repeated-support section, and the degree-12 factor is irreducible over
`F_7(u)`.  The degree-12 coefficient numerator/denominator degrees are at
most 32, so this is substantially smaller than carrying the raw degree-40
support quotient into subsequent boundary calculations.  Recovering `M`
from the lex basis and testing it in the degree-12 function field gives

```text
M is not a square in F_7(u)(support),
```

so adjoining `L` produces an irreducible degree-24 signed cover modulo 7.

The same good reduction certifies irreducibility of the coefficientwise
characteristic-zero P8 pullback.  No characteristic-zero primary
decomposition is needed for this conclusion.

## Exact genera

The exact fields before P8 base change have

```text
support genus = 5,
signed genus  = 10.
```

For the signed double cover, the divisor of `M` has exactly two odd
degree-`1` places, both over `e=0`; Riemann--Hurwitz therefore gives genus
`10`, not `9`.  This lower-genus signed curve is an exact
characteristic-zero target.

For the P8 base change, the exact relative different has degree `32` and its
pullback through `e=e(u)` has degree `168`.  Riemann--Hurwitz therefore gives

```text
support genus = 73,
signed genus  = 145.
```

Exact ramification of `e=e(u)` shows index `2` over `e=0`, making the
pulled-back divisor of `M` even; hence the P8 sign lift is unramified over its
support curve.  Exact irreducibility certifies connectedness, so its genus is
`2*73-1=145`.  Independent computations in characteristics `7` and `13`
reproduce the same pair `73/145`.

## Rigorous finite-field masks

For a good reduction prime `p != 3`, rational independent 3-torsion injects
into `J(F_p)[3]`.  Therefore a necessary condition is

```text
dim_F3 J(F_p)[3] >= 2.
```

For each projective residue `tau in P^1(F_p)`, the sieve constructs

```text
f_tau = x*(3*x^2+(a-3)*x+2)*(2*x^2-3*x+(a+3))
```

and computes the exact finite abelian group of its Jacobian.  Smooth fibers
with 3-rank below two are rejected.  Singular fibers, parameter poles, and
chart-boundary residues are retained, so the filter has no false negative
caused by the chosen model.

At the five primes

```text
p=7,17,19,23,41,
```

there is no good residue with 3-rank at least two.  Thus any global extra-3
point is forced into a bad-reduction residue disk at every one of these
primes.  For example, at `p=7` the only retained projective residues are

```text
tau = 1,2,4,6 mod 7;
```

the four good residues all have exactly one 3-direction.

Using all primes

```text
7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71
```

leaves only the two `e=0` boundary parameters through height 1000.

As a separate audit, using only primes through 43 at height 500 left 107
parameter values.  After removing `e=0` and duplicate `a`-fibers, exact
`TorsionSubgroup` computations were run on all 84 distinct curves.  Every
one had source torsion

```text
[2,2,6],
```

so none had the second rational 3-direction.  Adding primes 47 through 71
then removed all of those nonboundary survivors at the mask stage.

The projective parameter `tau=infinity` gives the same `e`-fiber as
`tau=12`; it is therefore covered by the finite scan and has no extra
3-direction.

## Forced local disks at `7` and `17`

The four retained residues at `p=7` split as follows:

```text
tau=2,-3: endpoint zero disks;
tau=1,-1: parameter-pole disks.
```

The primitive E3 layer, every uniformly scaled E3 depth, and the previously
enumerated endpoint fan through `v(e)=8` are empty in the zero disks, but
valuation signatures beyond that bounded fan are not exhausted.  Each
punctured pole disk has two smooth signed branches.  An explicit Weil-pairing
calculation proves that they are on the degree-`12` orthogonal support and
degree-`24` signed cover.  Hence `p=7` does not close the route.

At `p=17`, the ten retained disks have three types.  The two endpoint disks
have no uniformly scaled E3 branch, with other weight signatures unresolved;
the four `a=-2` common-root disks have no affine signed contact point, with
nonintegral blowups unresolved; and the four `DC=0` disks each have two
smooth signed orthogonal branches.  Thus `p=17` also does not close the
route.  See the two dated bad-disk notes for the precise scope of every
negative statement.

## Reproduction

From `torsion_jac`:

```bash
magma -b p:=7  code/contact6_m612_p8_relative3_modular.m
magma -b p:=11 code/contact6_m612_p8_relative3_modular.m
magma -b p:=13 code/contact6_m612_p8_relative3_modular.m
magma -b p:=7  do_genus:=true code/contact6_m612_p8_relative3_modular.m
magma -b p:=13 do_genus:=true code/contact6_m612_p8_relative3_modular.m

magma -b code/contact6_m612_relative3_exact_reconstruct.m
magma -b relative_only:=true do_signed:=false \
  code/contact6_m612_p8_relative3_exact_genus.m
magma -b p8_rh_only:=true \
  code/contact6_m612_p8_relative3_exact_genus.m

magma -b mode:=search height:=1000 prime_bound:=71 max_exact:=500 \
  code/contact6_m612_p8_extra3_residue_sieve.m

magma -b code/contact6_m612_p8_p7_bad_disks.m
magma -b code/contact6_m612_p8_p17_bad_disks.m
```

The exact signed genus and P8 base-ramification checks are in the continuation
scripts `code/contact6_m612_relative3_exact_m_divisor_continuation.m` and
`code/contact6_m612_relative3_exact_p8_base_ramification_continuation.m`.

The modular resolvent jobs take about 12--14 seconds each.  The height-1000
mask scan takes about 24 seconds on the current machine.

## Next exact step

Work with the exact degree-24 signed cover, not the full degree-40 quotient.
The first global target should be the exact genus-`10` signed curve over the
`e`-line: determine its rational points, Jacobian, and useful `S3` quotients
before attacking the much larger P8 pullback.  In parallel, finish the exact
global analysis of the genus-`145` P8 signed cover and the forced
bad-reduction disks at `p=19,23,41`; `p=7` and `p=17` already have positive
target branches and cannot obstruct the route alone.
The unresolved high-weight endpoint fans and nonintegral common-root blowups
must not be reported as empty.  Another undirected height increase is lower
value than these global and local component calculations.
