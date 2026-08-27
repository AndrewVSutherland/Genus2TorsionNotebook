# Fast sign quotient for the `b2=0` norm cover

The direct six-variable saturation in
`code/m12_general5_b2zero_geometry.m` is avoidable.  On `b1 != 0`, the
simultaneous sign change

```text
(b0,b1) -> (-b0,-b1)
```

is a free involution.  It sends the order-5 class to its negative and fixes
the underlying `M(12)` curve.  The two open points in each of the existing
`p=7,11,13` logs are exactly one orbit of this involution.

## Root-of-`B` normalization

Put `c=b0/b1` and `z=x+c`, so `B=b1*z`.  Evaluating

\[
 A^2-B^2F=q^5
\]

at `z=0` gives a rational coordinate

\[
 d=\frac{A(-c)}{q(-c)^2},\qquad q(-c)=d^2.
\]

On `d != 0`, set `z=dZ` and write

\[
 q=d^2(Z^2+eZ+1),\qquad b_1=\lambda d^4,
 \qquad G(Z)=F(dZ-c)=\sum_{i=0}^5g_iZ^i.
\]

The sign involution is `lambda -> -lambda`; its quotient coordinate is
`tau=lambda^2`.  Let

\[
 p=\frac{5e}{2},\quad q_0=\frac52+\frac{15e^2}{8},\quad
 \delta=(2-e)^3
\]

and

\[
\begin{aligned}
 \ell_3&=pg_0-g_1,&\ell_4&=q_0g_0-g_2,\\
 \ell_5&=q_0g_0-g_3,&\ell_6&=pg_0-g_4,\\
 \ell_7&=g_0-g_5.
\end{aligned}
\]

After coefficients `0,1,2,8,9,10` are forced, the remaining norm
coefficients are

\[
\begin{aligned}
E_3&=\frac{5\delta}{8}+\tau\ell_3,\\
E_4&=\frac{5\delta(19e-6)}{64}+\tau\ell_4
       +\frac{\tau^2g_0^2}{4},\\
E_5&=\frac{\delta(32e^2-33e+58)}{32}+\tau\ell_5,\\
E_6&=\frac{5\delta(19e-6)}{64}+\tau\ell_6,\\
E_7&=\frac{5\delta}{8}+\tau\ell_7.
\end{aligned}
\]

On the squarefree-`q` open set, `delta*ell3 != 0`, and eliminating `tau`
gives only four equations in `(b,w,c,d,e)`:

\[
\begin{aligned}
C_1&=\ell_3-\ell_7,\\
C_2&=(32e^2-33e+58)\ell_3-20\ell_5,\\
C_3&=(19e-6)\ell_3-8\ell_6,\\
C_4&=4(19e-6)\ell_3^2-32\ell_4\ell_3+5\delta g_0^2.
\end{aligned}
\]

Recovery is

\[
 \tau=-\frac{5(2-e)^3}{8\ell_3}.
\]

## Four-variable quotient chart

If `g0 != 0`, then `C1` is linear in `e`.  Put

\[
 h=g_0+g_1-g_5,\qquad L_0=g_0-g_5.
\]

Then

\[
 e=\frac{2h}{5g_0},\qquad
 \tau=-\frac{(5g_0-h)^3}{25g_0^3L_0}.
\]

The quotient curve is cut out in `(b,w,c,d)` by

\[
\begin{aligned}
K_3={}&(38h-30g_0)L_0-40g_0(h-g_4),\\
K_2={}&(128h^2-330hg_0+1450g_0^2)L_0
 -1250g_0^3-150g_0h^2+500g_0^2g_3,\\
K_4={}&(95h-75g_0)L_0^2-10(25g_0^2+3h^2)L_0
 +100g_0g_2L_0+(5g_0-h)^3.
\end{aligned}
\]

The rational order-5 class itself exists only on the double cover

\[
 25g_0^3L_0\lambda^2+(5g_0-h)^3=0.
\]

These identities are implemented with symbolic reconstruction assertions in
`code/m12_general5_b2zero_rootquotient.m`.  A pure-Python finite-field
cross-check is `code/m12_general5_b2zero_rootquotient_local.py`.

## Local evidence

The reduced implementation reproduces the full Magma counts at the three
available primes and extends them cheaply:

```text
p     quotient points lifting to lambda     signed points
7                    1                           2
11                   1                           2
13                   1                           2
17                   2                           4
19                   5                          10
23                   5                          10
29                   9                          18
```

There are no `g0=0` open lifts for `p=7,11,13,17,19`.  The growth after
`p=13` is evidence for a genuine positive-dimensional cover rather than a
zero-dimensional scheme.  It does not distinguish a rational/low-degree
component from a higher-genus component: these counts include only quotient
points for which `tau` is a square, and the open set has many punctures.
Thus the three singleton lifting-orbit counts at very small primes alone do
not support a low-degree-component inference.

## Caveats and component strategy

The conditions `Disc(q) != 0` and `Res(q,F) != 0` in the original script are
genericity conditions, not consequences of a valid Mumford representation.
Repeated-support divisors and divisors meeting a Weierstrass point can be
valid.  They require separate valuation/divisor checks because the norm
identity can split its vanishing between the two hyperelliptic branches.
Thus the current saturation must not be advertised as the complete
`b2=0` locus until those boundaries (and `Res(B,q)=0`) are audited.

For the generic component, first saturate only by the cheap factors

```text
b*w*(b-1)*(2*b-1)*d*g0*L0*(h-5*g0)*(h+5*g0).
```

Decompose the resulting three-equation ideal, and only then test whether a
component is contained in `Disc(F)` or `Res(q,F)`.  Multiplying the large
degree-38 discriminant and degree-16 resultant into the first saturation is
needlessly expensive.  The `g0=0` chart is also cheap: it satisfies

```text
g0=0, g1=g5 != 0, g2=g4,
e=(6*g1+8*g4)/(19*g1), plus C2.
```

The root swap on the `M(12)` presentation is not an additional involution of
this slice: as a binary quadratic, a linear `B` has roots at infinity and at
`-c`; the relevant Mobius root swap sends infinity to a finite point and
generically produces a quadratic `B`.  It belongs to the full `b2,b1,b0`
surface, not to `b2=0`.

## Boundary reinterpretation: a point-order-10 lane

The `Disc(q)=0` and `Res(q,F)=0` boundaries should be attacked rather than
discarded.  They have a direct point-torsion interpretation, subject to the
local-valuation caveat below.

If `q=(x-r)^2` and the Mumford congruence is the tangent congruence at a
non-Weierstrass rational point `P`, then

\[
 [q,v]=2(P-\infty).
\]

Consequently an exact order-5 class `[q,v]` makes `P-infinity` have order 5
or 10.  The order-5 alternative is part of the existing point-contact-5
lane; the order-10 alternative is a new confluent point-contact lane.

If `q` is squarefree but shares one root with `F`, its common factor over
`Q` is generically linear and gives a rational Weierstrass point `W`.  The
other point in the Mumford support is a rational point `P`, and

\[
 D=[q,v]=(W-\infty)+(P-\infty).
\]

Writing `T=W-infinity`, if `D` has exact order 5 then

\[
 E=P-\infty=D-T,\qquad 5E=T.
\]

so `E` has exact order 10.  Thus the shared-Weierstrass boundary naturally
encodes a rational point-order-10 class whose order-2 part is `T`.

This is especially well matched to the compact `M(12)` family.  On a fiber
with rational 2-rank one, every nonzero rational 2-class is the same class.
If `R` is the marked order-12 class, then its order-2 element is `6R`; hence

\[
 5E=6R.
\]

The 2-primary subgroup generated by `R` and `E` is therefore the cyclic
order-4 subgroup generated by the 2-primary part of `R`, while their 3- and
5-primary parts are cyclic of orders 3 and 5.  It follows that

\[
 \langle R,E\rangle\simeq \mathbf Z/60.
\]

For a repeated-`q` point of exact order 5, combining it with `R` also gives
cyclic order 60 because 5 and 12 are coprime.  For a repeated-`q` point of
order 10, the same common-order-2 check `5E=6R` is required; rational
2-rank one makes it automatic.

These conclusions are not licensed by the raw norm identity alone on the
boundary.  When support points collide or lie above a branch point, the
vanishing of `A+B*y` and its conjugate can merge or split, so the generic
argument assigning all five zeros to one branch no longer applies.  Every
candidate must be rebuilt as a Mumford divisor, checked by exact Cantor
arithmetic to have order 5 or 10 as claimed, checked for the common order-2
relation, and finally checked to have exact torsion `[60]` and a geometric
simplicity certificate.  Fibers where the quartic factor of `F` splits can
have rational 2-rank greater than one, in which case order 12 plus order 10
may instead generate a group with an extra independent order-2 factor.
