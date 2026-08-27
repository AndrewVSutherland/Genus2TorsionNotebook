# Shared-Weierstrass `M(12)+5`: exact cyclic-`60` audit

## Outcome

The shared-Weierstrass boundary is locally alive and its divisor
interpretation is exact, but no rational cyclic-`60` example was found.
Exchanging the two rational roots in the split `M(12)` presentation also
identifies this lane birationally with the already studied point-contact-`5`
lane; it is not a genuinely new global cover.
The computation reduces the boundary to one irreducible high-degree curve
and proves the following bounded exclusions on its smooth open:

- no point with both compact parameters `b,z` of height at most `100`;
- more strongly, no point with `height(b) <= 5000`, with **no height bound on
  the other square-cover coordinate**.

The latter search checks `30,401,828` reduced rational `b` values and is empty
after the exact projection masks through `p=103`.

## Exact divisor identity

Write the compact odd model as

\[
 y^2=F(x),\quad L=b+(2b-1)x,
\]

\[
 F=L\bigl(LH_0^2+4b(1+x)^2(wL-x^2)\bigr),\qquad
 H_0=x+w(1+bx).
\]

The visible Weierstrass point is

\[
 W=\left(-\frac b{2b-1},0\right).
\]

Set `X=L` and multiply the right side by the square `(2b-1)^4`.  The
resulting square-equivalent model is

\[
 Y^2=\widehat F(X)=XG(X),\qquad \deg G=4,
\]

with `W=(0,0)`.  For a non-Weierstrass point `P=(R,H(R))`, put
`E=P-infinity` and `T=W-infinity`.  On the smooth open,

\[
 5E=T
\]

if and only if there is a cubic `H`, with `H(0)=0`, such that

\[
 \widehat F-H^2=\kappa X(X-R)^5.                 \tag{1}
\]

Indeed,

\[
 \operatorname{div}\left(\frac{Y-H}{X}\right)
   =5P-W-4\infty.
\]

This proves the desired relation directly and avoids applying the generic
Mumford norm interpretation on a common-support boundary.

There is also a useful valuation check on the old norm equation.  If
`q=(x-omega)(x-r)` shares the simple root `omega` with `F`, then
`A^2-B^2F=q^5` forces

\[
 \operatorname{ord}_{x-\omega}(A)\ge3,\qquad
 \operatorname{ord}_{x-\omega}(B)=2.
\]

Thus this boundary cannot occur in the previously studied linear-`B`
subcover; it lies naturally in the quadratic-`B` chart.

## Small sign quotient and rational square cover

Write

\[
 H=X(aX^2+cX+d),\quad s=a^2,\quad C=ac,\quad D=ad,
\]

and `G=sum(g_i X^i)`.  Cancelling `X` in (1) gives the five equations

\[
\begin{aligned}
 g_4-2C-5sR&=0,\\
 sg_3-C^2-2sD+10s^2R^2&=0,\\
 sg_2-2CD-10s^2R^3&=0,\\
 sg_1-D^2+5s^2R^4&=0,\\
 g_0-sR^5&=0,
\end{aligned}
\]

where

\[
 g_0=-4b^3(b-1)^2.
\]

Consequently `s` is a rational square precisely when `-bR` is a square.
Two useful birational square-cover coordinates are

\[
 R=-\frac{z^2}{b},\qquad
 a=\frac{2b^4(b-1)}{z^5},
\]

and, with `t=b/z`,

\[
 R=-\frac b{t^2},\qquad
 a=\frac{2(b-1)t^5}{b}.
\]

After solving the first two equations linearly for `C,D`, the last two are
of degrees `2` and `4` in `w`.  Exact reconstruction checks the full identity
(1).

## Geometry

In the `(b,w,z)` model the two primitive open equations have shapes

```text
N2: total degree 28, degree_w 2,  74 terms
N1: total degree 40, degree_w 4, 165 terms
```

Eliminating `w` factors as

```text
z^12 * b^62 * (b-1)^22 * (b-1/2)^8 * G(b,z),
```

where `G` is irreducible over `Q`, has `446` terms, and has bidegree

```text
(degree_b, degree_z) = (64,58).
```

It is even in `z`; the quotient in `(b,z^2)` remains irreducible and has
bidegree `(64,29)`.  There is no rational or low-degree component hidden by
the boundary factors.  Two endpoint coefficients are particularly useful:

\[
 \operatorname{lc}_z(G)=-\frac{64}{25}(b-\tfrac12)^6,
 \qquad G(b,0)=b^{42}(b-1)^{22}.                 \tag{2}
\]

Thus at every prime other than `2,5`, a point with good `b` cannot disappear
at either `z=0` or `z=infinity`.  This makes projection of the affine modular
masks to the `b`-line rigorous.

The birational `t=b/z` model gives the same irreducible component with
bidegree `(22,58)` in `(b,t)`; its endpoints are supported only on the same
known `b=0,1,1/2` boundaries.

## Root swap: equivalence with point-contact `5`

In the pre-normalized split model, let `w_0` be the selected rational root of
`T+1` and put `b=a*w_0`.  The other root is

\[
 w_1=\frac1a-w_0=\frac{w_0(1-b)}b.
\]

Choosing `w_1` instead gives the involution

\[
 (b,w)\longmapsto
 \left(1-b,\frac{w(1-b)}b\right).              \tag{3}
\]

It exchanges the old point at infinity `O` with the visible Weierstrass
point `W`.  If `E=P-O`, `T=W-O`, and `5E=T`, then in the swapped
presentation `O_new=W` and

\[
 E_new=P-O_new=E-T,\qquad 5E_new=0.
\]

Since `E` has exact order `10`, `E_new=E-5E=-4E` has exact order `5`.
Therefore the shared-Weierstrass point-order-`10` construction is precisely
the point-contact-`5` construction after root swap.  Formula (3) changes the
height of `b` by at most a factor of two, so it does not turn arbitrarily
high points into bounded ones, but it rules out treating this as an
independent new cover.

## Local exact checks

The direct cubic-Taylor covariants were enumerated over finite fields and
then checked in the Jacobian, rather than trusted merely as polynomial
identities:

```text
p    exact point classes E of order 10    order(D12+E)
7                   1                         60
11                  2                         60
13                  1                         60
17                  3                         60
19                  6                         60
```

Every sample satisfies

```text
Order(E)=10,  5*E=T,  Order(D12)=12,  Order(D12+E)=60.
```

Hence there is no small-prime obstruction to the construction.

## Bounded rational searches

The first projected search used exact `(b,z)` masks through `p=71`:

```text
height(b),height(z) <= 100
148,181,928 pairs tested
0 survivors
```

Using (2), the stronger one-coordinate search leaves `z` completely
unbounded.  Exact masks through `p=103` give

```text
height(b) <= 5000
30,401,828 reduced nonboundary b-values tested
0 survivors
```

The last survivor through `p=97`, `b=-2419/3437`, is removed by `p=103`.
This is a bounded exclusion in the chosen compact chart, not a proof that the
high-degree curve has no rational points.

The `b`-only sieve is rigorous rather than heuristic.  It passes residues
where the denominator of `b` is divisible by `p`, or where `b` reduces to
`0,1,1/2`.  Otherwise (2) forces `z` to be a p-adic unit.  The leading
`w^4` coefficient in the second reconstructed equation is a nonzero
open-unit multiple of `-b^4*(2b-1)^4/4`, so `w` cannot reduce to infinity.
The masks enumerate every affine `w` and retain singular reductions; hence
no good rational point is lost.

## Reproducibility

- `code/root_m12_point10_covariant_probe.m`: direct covariants and exact
  finite-field Jacobian checks;
- `code/m12_shared_weierstrass_order10_geometry.m`: five-equation sign
  quotient and full reconstruction;
- `code/m12_shared_weierstrass_order10_param.m`: rational `z` square cover,
  resultant, irreducibility, and endpoint factors;
- `code/m12_shared_weierstrass_order10_tparam.m`: lower-`b`-degree `t=b/z`
  model;
- `code/m12_shared_weierstrass_order10_sieve.py`: finite projection masks and
  the height-100 `(b,z)` search;
- `code/m12_shared_weierstrass_order10_b_sieve.py`: rigorous `b`-only search;
- `data/m12_shared_weierstrass_order10_local_2026_07_14.txt`;
- `data/m12_shared_weierstrass_order10_geometry_2026_07_14.txt`;
- `data/m12_shared_weierstrass_order10_h100_2026_07_14.txt`;
- `data/m12_shared_weierstrass_order10_b_h5000_2026_07_14.txt`.

## Assessment

The shared-Weierstrass construction is mathematically valid and locally
viable, but it is root-swap equivalent to the point-contact-`5` cover and the global locus is a single high-degree curve with no point in
substantial bounded searches.  It should not receive a larger blind height
search next.  A future continuation would need arithmetic of this curve (a
lower-genus quotient, descent, or Mordell--Weil/Chabauty information), not a
wider three-parameter enumeration.
