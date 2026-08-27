# Full contact-5/order-20 family plus rational 3-torsion

## Why this is a distinct order-60 lane

The family

\[
h_t=1+tx+\frac{t^2-1}{2}x^2,
\qquad
f_t=h_t^2-\frac{(t+1)^4}{4}x^5
\]

has a marked 5-torsion class and a marked order-4 class, hence a rational
point of order 20 on every smooth fiber.  Its discriminant is

\[
\frac{(t+1)^{11}(t+3)^4
(32t^3+152t^2+173t+37)}{256}.
\]

Generically

\[
f_t=(x-1)\,q_4(x)
\]

with irreducible quartic `q4`, so the generic rational 2-rank is one.  A
rational 3-class therefore gives an element of order 60 without forcing an
extra rational 2-direction.  Exact cyclic `[60]` remains possible.

The repository previously searched two proper subloci:

- the independent-extra-2 loci, targeting `[2,20]+3` and hence groups such
  as `[2,60]`;
- the order-40 cover, targeting `[40]+3` and order divisible by 120.

The older larger-torsion scan exact-tested only fibers where the residual
quartic was reducible.  Consequently the full one-parameter `t+3` problem on
the generic `[20]` locus had not been searched.

## Exact finite masks

At a prime `p != 2,3` of good reduction, reduction is injective on rational
3-torsion.  Hence a rational 3-class requires

\[
3\mid \#J_t(\mathbf F_p).
\]

Conversely, for a finite abelian group, divisibility by 3 is equivalent to
the existence of a rational element of order 3.  Thus this is the exact
finite-field existence mask, not merely a loose point-count heuristic.

The counts can be computed dependency-free from

\[
\#J(\mathbf F_p)=\frac{N_1^2+N_2}{2}-p,
\]

where `N1=#C(Fp)` and `N2=#C(Fp^2)`.  As a control, every smooth residue in
the computation has `20 | #J(Fp)`, exactly as the marked order-20 class
requires.

The first masks are:

| `p` | good `t` | allowed `t` with `3 | #J` | bad `t` |
|---:|---:|---|---|
| 7 | 5 | `{5}` | `{-3,-1}` |
| 11 | 9 | `{0,2}` | `{-3,-1}` |
| 13 | 11 | `{0,2,6}` | `{-3,-1}` |
| 17 | 14 | `{2,8,13,15}` | `{3,-3,-1}` |
| 19 | 16 | `{2,3,10,12,13,14,15}` | `{5,-3,-1}` |

Every tested prime through 101 has at least one smooth allowed residue.  The
full `t+3` locus is therefore locally alive; unlike several extra-2 and
order-40 subloci, it is not forced entirely onto the boundary at one prime.

Two rational-height results must be kept distinct:

1. Using masks only through `p=47`, the naive height-1000 box has 55
   survivors.  Two are `t=-1,-3`; the other 53 are smooth mask survivors.
   The first smooth survivor already appears at `t=-194/77` in height 200.
2. Using all 19 primes
   `7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79`, the full
   height-1000 box contains 1,216,767 reduced rationals and leaves exactly
   `t=-1,-3`.  Both are singular.  Thus there is no smooth rational fiber
   with a 3-part in this stated height box.

The 53 weaker-mask survivors are useful regression inputs, but every one is
killed by a later good prime.

## Cubic-contact quotient curve

### Low-degree `(s,z)` normalization

For computation, first put

\[
s=\frac{t-1}{t+1},\qquad
t=\frac{1+s}{1-s},\qquad
x=(1-s)z.
\]

Then

\[
g_s(z)=f_t((1-s)z)
=\bigl(1+(1+s)z+2sz^2\bigr)^2-4(1-s)z^5.
\]

This is the preferred chart: no `y`-scaling is needed, the constant
coefficient remains one, and

\[
\begin{aligned}
c_0&=1,&c_1&=2(s+1),&c_2&=s^2+6s+1,\\
c_3&=4s(s+1),&c_4&=4s^2,&c_5&=4(s-1).
\end{aligned}
\]

All coefficients now have degree at most two.  Moreover,

\[
\operatorname{disc}(g_s)=256(1-s)^2(2-s)^4
(8s^3-59s^2-18s+197).
\]

Here `s=1` is the `t=infinity` degree-drop boundary, `s=2` is the singular
fiber `t=-3`, and `s=infinity` is the singular fiber `t=-1`.  The last one
must be remembered separately from the affine `s`-chart.

The rational Weierstrass point `x=1` becomes `z=1/(1-s)`, so generic factor
type remains `[1,4]`.

In this chart the quotient identity is

\[
H^2-q^3=M g_s,
\]

and the constant square equation is simply

\[
E^2-V^3=M=L^2.
\]

This is substantially smaller than working directly in `t`.

### Recursive quotient equations

Write, in the normalized `(s,z)` chart,

\[
q=z^2+Uz+V,
\qquad
H=z^3+Az^2+Bz+E.
\]

A rational 3-class is represented by

\[
H^2-q^3=Mg_s,
\qquad M=L^2\ne0.
\]

Let `ci` be the coefficient of `z^i` in `g_s`.  Explicitly,

\[
\begin{aligned}
c_0&=1, & c_1&=2(s+1), & c_2&=s^2+6s+1,\\
c_3&=4s(s+1), & c_4&=4s^2, & c_5&=4(s-1).
\end{aligned}
\]

The coefficients in degrees 5, 4, and 3 determine

\[
\begin{aligned}
A&=\frac{Mc_5+3U}{2},\\
B&=\frac{Mc_4+3(U^2+V)-A^2}{2},\\
E&=\frac{Mc_3+U^3+6UV-2AB}{2}.
\end{aligned}
\]

The quotient contact curve in `(s,M,U,V)` is cut out by

\[
\begin{aligned}
G_2&=B^2+2AE-3(U^2V+V^2)-Mc_2=0,\\
G_1&=2BE-3UV^2-Mc_1=0,\\
G_0&=E^2-V^3-Mc_0=0.
\end{aligned}
\]

This is a curve before imposing that `M` is a rational square.  The open
part should be saturated by

\[
M(s-1)(s-2)(8s^3-59s^2-18s+197)
(U^2-4V)\operatorname{Res}_z(q,g_s).
\]

The square pullback is `M=L^2`.

### Repeated support must be retained

The divisor support polynomial need not be squarefree.  The locus

\[
U^2-4V=0,
\qquad q=(z-\rho)^2,
\qquad U=-2\rho,\qquad V=\rho^2
\]

represents a legitimate Mumford divisor `2P-2 infinity`.  Indeed, the
repository's validated order-30 cubic contact is of exactly this nonreduced
form.  It must therefore be audited before `disc(q)` is saturated.

For the exact-cyclic target, however, this entire branch is ruled out by its
factorization.  The contact identity becomes

\[
M g_s=\bigl(H-(z-\rho)^3\bigr)\bigl(H+(z-\rho)^3\bigr).
\]

Because `H` and `(z-rho)^3` are both monic cubics, the first factor has
degree at most two and the second has degree three.  Thus `g_s` has a
rational `2+3` factorization.  It also has the fixed rational Weierstrass
root `z=1/(1-s)`, so one of these factors splits again.  On a smooth fiber
there are at least three rational irreducible branch factors, hence rational
2-rank at least two.  A repeated-support hit can target `[2,60]`, but never
exact cyclic `[60]`.

Consequently the generic exact-cyclic computation may rigorously saturate by
`disc(q)` after recording this argument; solving the repeated scheme is only
a side computation for noncyclic groups.

The sub-slice `rho=0`, namely `q=z^2`, is itself impossible on the open
`M=L^2 != 0` locus.  Indeed, writing `E=eps*L`, the equations in degrees one
and two give

\[
B=\epsilon L(s+1),\qquad A=2\epsilon Ls.
\]

The recursive degree-three formula then gives `E=0`, contradicting
`E=eps*L` and `L != 0`.

For completeness, on the remaining repeated-support branch the workbench
imposes

\[
G_2=G_1=G_0=0,\qquad M=L^2
\]

in `(s,M,rho,L)` and saturates `M*rho`, the singular family factors, and
`Res(q,g_s)`.

## Exact elimination of the `V=0` branch

When `V=0`, the support polynomial is `q=z(z+U)`.  The constant equation
gives `E=eps*L`, with `eps=+1` or `-1`.  On the open `L*U != 0` locus, the
three remaining equations simplify successively to

\[
B=\epsilon L(s+1),\qquad
A=2\epsilon Ls,\qquad
U^3=2\epsilon L.
\]

The recursive formula for `B` then gives

\[
s=\frac3U-1,\qquad L=\frac{\epsilon U^3}{2}.
\]

Substitution in the recursive formula for `A` leaves

\[
\begin{aligned}
0={}&2U^5-3U^4-2U^2+6U-3\\
  ={}&(U-1)^3(2U^2+3U+3).
\end{aligned}
\]

The quadratic factor has discriminant `-15`.  Therefore the only rational
solution is

\[
U=1,\qquad s=2,\qquad L=\epsilon/2.
\]

It satisfies the contact identity with `M=1/4`, but `s=2` is the singular
fiber `t=-3`.  Thus the whole rational `V=0` branch contributes no smooth
order-60 candidate.

## Exact elimination of the `E=0` slice

On the square cover, `E=0` forces `V^3=-L^2`; this is another natural
low-degree slice, not the same as `V=0`.  Independent elimination in the
normalized `(s,z)` chart gives an `s`-resultant

\[
(s-1)^{24}P_{10}(s)P_{20}(s),
\]

where `P10` and `P20` are irreducible over `Q` of the indicated degrees.
They have no rational roots, and `s=1` is the degree-drop boundary.  Hence
the entire rational `E=0` square-cover slice is empty off boundary.

## The generic square cover `V != 0`

The constant equation is

\[
E^2-L^2=V^3.
\]

It is rational.  Put

\[
w=\frac{E+L}{V}.
\]

Then

\[
L=\frac{V(w^2-V)}{2w},\qquad
E=\frac{V(w^2+V)}{2w},\qquad
M=\frac{V^2(w^2-V)^2}{4w^2}.
\]

After substituting these formulas, impose the recursive equality for `E`
and the two equations `G2=G1=0`.  This gives three equations in
`(s,U,V,w)`, again a curve, with open boundary

\[
Vw(w^2-V)(U^2-4V)\operatorname{Res}(q,g_s)
(s-1)(s-2)(8s^3-59s^2-18s+197)\ne0.
\]

This formulation builds the rational-square condition into the geometry and
is preferable to a blind search in five contact variables.

## Fast projection to `(s,M)`

The low-degree equations admit a much faster elimination than a four-variable
primary decomposition.  The equations `G1` and `G2` are quadratic in `V`.
If

\[
G_1=a_2V^2+a_1V+a_0,
\qquad
G_2=b_2V^2+b_1V+b_0,
\]

then

\[
b_2G_1-a_2G_2=D_VV-N_V.
\]

Away from `D_V=0`, recover `V=N_V/D_V`, clear denominators, and eliminate
`U`.  The exact run is

[`code/contact5_order20_plus3_sx_fastelim.m`](../code/contact5_order20_plus3_sx_fastelim.m).

Put

\[
D=b_2a_1-a_2b_1,\qquad N=a_2b_0-b_2a_0,\qquad
C=a_1b_0-a_0b_1.
\]

For `R=N^2+DC`, direct polynomial identities give

\[
\begin{aligned}
D^2G_1(N/D)&=a_2R,\\
D^2G_2(N/D)&=b_2R.
\end{aligned}
\]

Since `D != 0` implies that `a2,b2` are not both zero, the two equations
`G1=G2=0` are equivalent on this chart to

\[
R=0,\qquad V=N/D.
\]

Consequently the exact generic plane projection is the single resultant

\[
\operatorname{Res}_U\bigl(R,D^3G_0(N/D)\bigr).
\]

The completed run first asserts both displayed polynomial identities.  This
exact resultant has total degree 288 and 12,143 terms.  As an independent
cross-check, it agrees exactly with the gcd of the two older pairwise
resultants: the reported extra degree is zero.  Thus there is no
``different `U` root'' or leading-coefficient artifact in the result below.

After removing the factors `M^56(s-1)^16`, the exact plane resultant factors
as follows:

- a degree-20 factor of multiplicity six, bidegree `(12,8)`, supported on
  the `D_V=N_V=0` denominator branch;
- one irreducible factor `P(s,M)` with bidegree `(56,40)`, total degree 96,
  and 2,107 terms.

The denominator branch is harmless for rational points.  Imposing
`D_V=N_V=0` together with both pairwise `V`-resultants against `G0` gives a
zero-dimensional scheme with no rational `(s,M,U)` points.  These equations
define a superset of the true denominator lifts, so emptiness of its rational
points is a rigorous exclusion.

Pulling back the other factor by `M=L^2` gives an irreducible polynomial of
bidegree `(56,80)`, total degree 136, still with 2,107 terms.  It is the sole
non-denominator component of the exact square-cover projection.  On the
generic open chart it recovers `U` and then `V=N/D`, so it is the genuine
generic cubic-contact curve, not merely a necessary sieve.  Exceptional
points of its plane model should still be lifted and checked in the original
equations individually.

## Reproducible geometry workbench

The preferred low-degree static script

[`code/contact5_order20_plus3_sx_geometry.m`](../code/contact5_order20_plus3_sx_geometry.m)

implements the `(s,z)` chart, including separate `vzero` and `repeat` modes,
generic quotient decomposition, `(s,M)` elimination, and the direct square
parametrization.  The earlier direct-`t` cross-check

[`code/contact5_order20_plus3_quotient_geometry.m`](../code/contact5_order20_plus3_quotient_geometry.m)

contains:

- symbolic coefficient assertions for the contact identity;
- a separate exact repeated-support mode for `q=(x-rho)^2`, including
  `rho=0`;
- the exact `V=0` factorization and singular control point;
- saturated decomposition of the quotient curve;
- lexicographic elimination to the `(t,M)` plane, followed by `M=L^2`;
- direct decomposition after the rational `V != 0` square parametrization.

The next computation should run the cheap `vzero` self-check, then decompose
the generic parametrized square cover.  The height-1000 finite exclusion is
strong evidence of sparsity, but it is not a global nonexistence proof; the
curve geometry is the decisive next step.

## Height-100000 sieve and exact survivor audit

The finite sieve was subsequently rewritten as the CRT-compressed C program

[`code/contact5_order20_plus3_height_sieve.c`](../code/contact5_order20_plus3_height_sieve.c).

It uses exact necessary masks through `p=113` in both the original parameter
`t` and the normalized parameter `s=(t-1)/(t+1)`.  At height `100000`, the
`t` chart reduces roughly 12.2 billion raw numerator/denominator pairs to
380,162,713 CRT candidates and 223,486,975 primitive candidates.  Besides the
two singular parameters, only four smooth parameters survive:

```text
-22485/10693, 75226/34441, 37925/71877, -69764/79661.
```

Fresh exact Magma reconstruction gives torsion `[20]` for all four.  A later
good prime kills the necessary 3-divisibility in every case, while D4
Frobenius certificates prove the four Jacobians geometrically simple.  The
`s` chart gives no additional smooth candidate.  Two earlier survivors of the
height-10000 sieve likewise have exact torsion `[20]`, not `[60]`.  Thus the
search is negative through height `100000` in two birational base charts, but
this remains a bounded-height statement.

## Fixed-Weierstrass normalization

The model

\[
 F_s(x)=((1-s)^2+(1-s^2)x+2sx^2)^2-4x^5
\]

has the fixed simple branch point `x=1`, with
`F_s'(1)=-4(s-2)^2`.  Shifting `z=x-1`, write

\[
 q=z^2+az+r^2,\qquad H=z^3+Az^2+Bz+r^3.
\]

The constant and linear contact equations, together with the rational-square
condition, have the rational parametrization

\[
 L=rj/2,\qquad A=(3a-r^2j^2)/2,
 \qquad B=(3ar^2-(s-2)^2j^2)/(2r).
\]

Only three sparse equations remain.  They have total degrees `9,8,10`, term
counts `10,15,17`, and degrees `2,3,2` in `a`.  Saturation by the smooth open
boundary gives a one-dimensional scheme with a 68-element Groebner basis.
The exact derivation and self-tests are in

[`code/contact5_order20_plus3_fixed_weierstrass.m`](../code/contact5_order20_plus3_fixed_weierstrass.m).

Eliminating `a` by exact quadratic recovery produces two equations in
`(s,r,j)`.  Their resultants have one large irreducible factor in each useful
projection:

- `P116(s,j)`, bidegree `(60,80)`, with 2,120 terms;
- `P44(s,r)`, bidegree `(22,40)`, with 671 terms and multiplicity two.

All equations are even in `j`.  On `P44`, putting `k=j^2` and dividing the
quartic recovery equation by the quadratic one leaves a linear remainder.
Consequently `k` is recovered rationally, and direct substitution modulo
`P44` verifies both original recovery equations.  After removing a common
factor, the exact sign cover is

\[
 P_{44}(s,r)=0,\qquad A_{43}(s,r)j^2+B_{43}(s,r)=0,
\]

where `A43` and `B43` have respectively 654 and 595 terms.  This is a much
smaller exact generic model than the earlier bidegree `(56,80)` projection.
The degree `40` over the base is consistent with the 40 pairs `\{T,-T\}` of
nonzero geometric 3-torsion points; the quadratic cover retains the sign
needed for an actual rational 3-torsion point.  This interpretation is a
structural consistency check, not yet a proof that the displayed plane model
has no exceptional components.

The remaining audit is finite and explicit: classify the low-degree resultant
factors, intersect the linear-recovery denominator with `P44`, and check every
exceptional lift in the original three equations.  Only after that audit may
`P44` plus its quadratic cover be treated as the global open contact curve.
The workbench is

[`code/contact5_order20_plus3_fixed_fastelim.m`](../code/contact5_order20_plus3_fixed_fastelim.m).
