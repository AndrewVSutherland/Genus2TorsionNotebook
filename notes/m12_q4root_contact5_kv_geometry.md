# The Q4-root contact-5 cover of M(12)

## Compact rational surface

Write the completed-square M(12) model as

\[
T=ax^2-x+r,\qquad
W=(T+1)Q_4,\qquad
Q_4=(x-r)^2(T+1)+4ax^2T.
\]

Force a rational root `v` of `Q4`, and put

\[
d=v-r,\qquad \alpha=av^2,\qquad k=\alpha/d,
\qquad D=1-k,\qquad q=2k-1.
\]

Then

\[
\frac{Q_4(v)}{d^2}=q^2-Dd,
\]

so, away from the degenerate boundary,

\[
d=\frac{q^2}{D},\qquad
a=\frac{kq^2}{Dv^2},\qquad
r=v-\frac{q^2}{D}.
\]

This is a rational two-parameter cover of M(12).  Its affine open is

\[
kDqv\ne0.
\]

Equivalently one can use the branch-free coordinate `q=2k-1`; the three
excluded `k`-values become `q=-1,0,1`.

Move `v` to infinity with `X=v/(x-v)`.  Define

\[
\begin{aligned}
P&=4kD^2X^2+(2kq^2-vD)X+kq^2,\\
R&=-q^2DX^2+(2kq^2-vD)X+kq^2,\\
M&=q^2X+vD,\\
Q&=M^2P+4kq^2D(X+1)^2R.
\end{aligned}
\]

The quartic coefficient of `Q` cancels identically.  Hence

\[
F=PQ
\]

is an odd quintic, generically of factor type `[2,3]`.  More precisely,

\[
P=DX^2(T+1),\quad R=DX^2T,\quad M=DX(x-r),\quad
F=D^4X^6W(v+v/X).
\]

The last multiplier is a square, so this normalization preserves the
Jacobian and all contact-square conditions.  The marked order-12 point is at
`X=-1`, where

\[
F(-1)=\big((vD-q^2)(vD-q^2+D)\big)^2.
\]

The generic `[2,3]` factorization gives exactly one rational 2-direction, so
it is compatible with an exact cyclic `Z/60` target.  This is the main reason
to use the Q4-root cover rather than imposing a split `T+1`.

## Exact contact equations

For a possible contact abscissa `u`, put

\[
A_i=F^{(i)}(u)/i!,\qquad Q_2=4A_0A_2-A_1^2.
\]

A quadratic `h` with

\[
F-h^2=\lambda(X-u)^5
\]

exists precisely on the open `A0 != 0` locus when

\[
E_3=8A_0^2A_3-A_1Q_2=0,
\qquad
E_4=64A_0^3A_4-Q_2^2=0,
\]

and `A0` is a rational square.  The square condition is the final double
cover over the curve cut out by `E3=E4=0` in `(k,v,u)`.

Two slices are especially transparent:

- At `u=-1`, the prospective order-5 point collides with the marked
  order-12 point.  Smooth solutions should therefore be absent; the geometry
  script compares this slice with the discriminant boundary.
- At `u=0`,

  \[
  F(0)=(kq^2)^2D(v^2D+4kq^2),
  \]

  so the remaining square cover is the rational conic

  \[
  z^2=D(v^2D+4kq^2).
  \]

## Projective finite-mask correction

For a quintic with leading coefficient `c5`, the generic leading
coefficients in `u` are

\[
\operatorname{lc}(E_3)=5c_5^3,
\qquad
\operatorname{lc}(E_4)=95c_5^4.
\]

A rational `u` reducing to infinity modulo `p` requires both homogenized
equations to vanish there.  Therefore a residue is a projective boundary only
when both degrees drop, not when either one drops.  In particular, modulo 19
only `E4` drops generically; `E3` remains nonzero at infinity, and `p=19`
supplies a valid affine mask.  The corrected count is

\[
212\text{ good pairs},\qquad 2\text{ allowed pairs},\qquad
149\text{ bad/boundary pairs}.
\]

The two allowed triples are `(k,v,u)=(15,2,0)` and `(18,9,11)` modulo 19.
An independent finite-field implementation reproduces all Magma counts for
`p=7,11,13,17,19`.

## Audited searches

The direct compact `k`-height-20 run used 511 rational values in each
coordinate, hence 261,121 pairs.  With masks at `p=7,11,13,17`, it had

- 47,845 mask survivors;
- 47,753 smooth survivors;
- no exact common contact root and no order-60 candidate.

The stronger branch-free `q`-height-50 run used 3,095 rational values in each
coordinate and exact masks at

\[
7,11,13,17,19,23,29,31,37,41,43,47.
\]

Its audited counters are

- 9,579,025 parameter pairs;
- 12,672 raw mask survivors;
- 295 states after removing `q=-1,0,1` and `v=0`;
- 247 singular states;
- 8 degree-drop states;
- 40 smooth states, all with no common `E3/E4` root;
- no contact point.

Smooth contact residues do occur at every tested prime, so this is not a
local obstruction.  It is a meaningful low-height negative result on what is
expected to be a curve.

## Next computation

The next affine step is

[`code/m12_q4root_contact5_kv_geometry.m`](../code/m12_q4root_contact5_kv_geometry.m).

It factors the contact covariants, strips `kDqv`, compares with the quintic
discriminant, eliminates `u`, factors the projected curve in `(k,v)`, and
tests low-degree components for rational points before imposing
`F(u)=square`.  This should be run after the separate `a=-1/4` boundary
elimination.  A larger undirected height search is secondary: the corrected
masks make it feasible, but the resultant geometry can reveal whether the
remaining curve has rational components, positive-genus components, or only
boundary points.
