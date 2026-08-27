# The contact-6 `[6,6]` core as a route to `[6,12]`

## Executive conclusion

The geometrically simple contact-6 `[6,6]` construction is an exact source
for `[6,12]`.  Its generic primary torsion is

```text
[6,6]_(2) = [2,2],       [6,6]_(3) = [3,3].
```

Therefore one must make exactly one of the three nonzero rational `2`-classes
divisible by `2`, and then verify that this direction stops at order `4`.
This changes the `2`-primary part from `[2,2]` to exact `[2,4]`, while leaving
the already-constructed `[3,3]` untouched, and hence gives `[6,12]`.

There are two versions of this move:

1. halve one of the three classes on the contact-6 source;
2. pass through the distinguished rational Richelot isogeny and halve one of
   the three dual-kernel classes on the codomain.

The distinguished dual of the one known simple `[6,6]` seed remains `[6,6]`,
so there is no immediate hit.  The route is nevertheless locally alive from
`p=11` onward.  At both `p=5` and `p=7`, however, the good affine chart has no
finite `[6,12]` point on either the source or the distinguished dual.  A
targeted rational search should therefore combine the core-cover solver with
the exact halving tests and treat the `5`- and `7`-adic boundary charts
explicitly.

Current computational status (2026-07-10): no geometrically simple target
has yet been found.  We do have an exact split control with full torsion
`[6,12]`, and the simple search has now been reduced to two compact, exact
halving covers.  The `T0` core pullback has been exhausted through slice
height `8`; the `TB` cover has been scanned directly through parameter height
`30`.  The detailed negative results below are useful rather than generic
height failures: every candidate is forced onto the `5`- and `7`-adic
boundary wall.

The focused boundary audit later on 2026-07-10 sharpens this substantially.
The smooth affine `DB/DC` square-cover cones have no point on the exact
`R3/R2` halving chart.  The only smooth weighted `E9+R3` disks lie on a
genus-one component whose normalization is

```text
Y^2 = 2*Z^4-6*X^4;
```

it has no `Q_2` or `Q_3` point.  The second endpoint component is rational
and gives a genuine simple near miss with source torsion `[2,2,6]` and dual
torsion `[2,12]`.  Its explicit parameterization is now the preferred lane
for imposing the missing second rational `3`-direction.  See
`notes/contact6_m612_focused_audit_2026_07_10.md`.


## The contact-6 core and its `[3,3]`

Put

```text
h = 1 + a*x + b*x^2 + x^3,
f = h^2 - (x-1)^6.
```

The factorization that matters for `[6,12]` is

```text
f = x*B*C,
B = (b+3)*x^2 + (a-3)*x + 2,
C = 2*x^2 + (b-3)*x + (a+3).
```

On the open set where `f` is squarefree, `b != -3`, and both quadratics are
irreducible, the branch-factor type is `[1,2,2]` and

```text
J[2](Q) = [2,2].
```

The marked point

```text
P = (1,a+b+2),       D = P-infinity
```

has order `6` on the exact hits.

An independent order-`3` class is constructed from

```text
q = x^2 + U*x + v^2,
h3 = (1/L)*x^3 + (B3/(2*L))*x^2
                  + (Delta3/(8*L))*x + v^3/L,
M = L^2,
```

where, writing `c_i=[x^i]f`,

```text
B3     = c5*M + 3*U,
Delta3 = 4*c4*M + 12*(U^2+v^2) - B3^2,
```

and the three equations are

```text
B3*Delta3 + 16*v^3 - 8*c3*M - 8*U^3 - 48*U*v^2 = 0,
Delta3^2 + 64*B3*v^3 - 64*c2*M
           - 192*(U^2*v^2+v^4) = 0,
Delta3*v^3 - 4*c1*M - 12*U*v^4 = 0.
```

Thus the core cover is an expected surface: five variables
`(a,b,M,U,v)` and three equations, together with the square cover `M=L^2`.
On the generic branch the last equation is linear in `a`; fixing `(b,v)`
reduces the computation to a zero-dimensional solve in `(M,U)`.  This is the
solver in `code/contact6_m36_core_slice_search.m`.


## The known simple seed and all three rational 2-classes

The verified simple specialization is

```text
a = 133/39,       b = -7/13,
L = 29/16,        U = -9/4,       v = 5/2,
```

with integral model

```text
y^2 = 11389248*x^5 - 18252000*x^4 + 42399396*x^3
      - 10288044*x^2 + 29659500*x.
```

Magma gives exact torsion `[6,6]`.  The stronger verifier
`code/contact6_m36_simple_hit_strong_verify.m` gives the D4 certificate

```text
p=37: T^4 - 10*T^3 + 54*T^2 - 370*T + 1369,
Galois group D(4).
```

This is the project-standard geometric-simplicity certificate.

A height-`20` simple-only core run checked `260099` `(b,v)` slices and found
`1102` algebraic slice solutions.  Its six verified hits are all this same
curve, represented by three cubic-contact triples and the sign of `L`; there
is no second simple `[6,6]` curve through height `20`.  The exact summary is
in `data/contact6_m36_core_slice_h20_simple_summary.txt`.

For this specialization the two quadratic factors, normalized to be monic,
are

```text
Bmon = x^2 + (1/6)*x + 13/16,
Cmon = x^2 - (23/13)*x + 125/39.
```

The three nonzero rational `2`-classes in Mumford form are exactly

```text
T0 = [x,0],
TB = [x^2 + (1/6)*x + 13/16,0],
TC = [x^2 - (23/13)*x + 125/39,0],
T0 + TB = TC.
```

An exact Jacobian computation also gives

```text
3*D = TB.
```

This identifies the old marked-class halving cover with one of the three
precise `[6,12]` covers: `D` is divisible by `2` if and only if `TB` is
divisible by `2`.  Indeed, one implication follows by multiplying a half of
`D` by `3`; conversely, if `2K=TB`, then

```text
2*(K+2D) = TB+4D = D.
```


## Exact halving conditions on the source

For a group with `2`-primary part `[2,2]`, the following patterns are the only
ones possible:

```text
no nonzero T is halved       -> exponent 2,
exactly one nonzero T halved -> [2,2^e] for some e >= 2,
two are halved               -> all three are halved, at least [4,4].
```

The last assertion follows because the sum of halves of two distinct classes
is a half of the third.  Thus exact `[6,12]` asks for exactly one of `T0,TB,TC`
to be divisible by `2`, with `e=2` rather than a further lift to order `8` or
higher (followed by an exact full-torsion check).

### Halving `T0=[x,0]`

Write

```text
g = B*C = c4*x^4+c3*x^3+c2*x^2+c1*x+c0,
```

where

```text
c4 = 2*(b+3),
c3 = b^2+2*a-15,
c2 = 2*a*b+22,
c1 = a^2+2*b-15,
c0 = 2*(a+3).
```

The class `T0` is divisible by `2` if and only if there exist rational
`u,v,m,n` such that

```text
g - x*(m*x+n)^2 = c4*(x^2+u*x+v)^2.
```

Equivalently,

```text
c0                         = c4*v^2,
c1 - n^2                   = 2*c4*u*v,
c2 - 2*m*n                 = c4*(u^2+2*v),
c3 - m^2                   = 2*c4*u.
```

The first necessary square condition is particularly cheap:

```text
v^2 = (a+3)/(b+3).
```

It already rejects the known simple seed, for which the right side is
`125/48`.  This square condition is not sufficient.

When the equations hold, a half is represented by

```text
[x^2+u*x+v, (m*u-n)*x+m*v].
```

#### Reduced exact `T0` cover

The remaining tangent equations reduce to one quartic.  To avoid collision
with the cubic-contact variable `v`, write the first square root as `omega`
and put

```text
s = b+3,
a = s*omega^2-3,
b = s-3,
W = m^2.
```

Retain the coefficient names

```text
c3 = b^2+2*a-15,
c2 = 2*a*b+22,
c1 = a^2+2*b-15,
```

and define

```text
N  = c1-omega*c3+omega*W,
R8 = 8*s*(c2-4*s*omega)-(c3-W)^2,
F  = R8^2-256*s^2*W*N.
```

On the open chart `s*m != 0`, the exact criterion is

```text
F(s,omega,W)=0,       W=m^2 in Q.
```

The remaining tangent variables are then recovered rationally as

```text
u = (c3-W)/(4*s),
n = R8/(16*s*m).
```

Indeed, the second coefficient equation gives `n^2=N`, while the third gives
`2*m*n=R8/(8*s)`; eliminating `n` is exactly `F=0`.  Conversely the displayed
reconstruction proves sufficiency.  The branch `m=0` must be retained
separately as a lower-degree tangent boundary.

Magma finds that `F` is irreducible in `Q[s,omega,W]`, has degree `4` in `W`,
total degree `12`, and `69` terms.  After `W=m^2`, it remains irreducible of
degree `8` in `m`.  Thus the initial square layer is a rational
parameterization, but the full generic halving surface does not acquire an
obvious rational parameterization from this elimination.

The special slice `omega=1` factors into two quadratics in `W` and can be
parametrized completely.  It has `a=b`, however, and

```text
x^6*f(1/x)=f(x).
```

The involution `(x,y) -> (1/x,y/x^3)` makes its Jacobian geometrically split.
It is therefore a useful formula check but not a source of simple `[6,12]`.
The exact derivation and factor checks are in
`code/contact6_m612_t0_reduced.m`.

#### Pulling back the `[3,3]` core

The direct intersection is obtained simply by substituting

```text
a=s*omega^2-3,  b=s-3,  M=L^2
```

in the three cubic-contact equations and adjoining

```text
F(s,omega,m^2)=0.
```

This gives four equations in

```text
(s,omega,m,L,U,v3),
```

so it is again an expected surface.  A concrete zero-dimensional slicing
strategy is to fix `(omega,v3)` and solve for `(s,m,L,U)`.  This enforces the
halving algebraically before any exact torsion computation and is the cleanest
alternative to enumerating core points and applying `IsDivisibleBy(T0,2)`.

The optimized solver `code/contact6_m612_t0_core_slice_search.m` fixes
`(omega,v3)`, solves the zero-dimensional cubic-contact system first, and
then solves the independent quartic in `W`.  At height `8` it found

```text
7140 / 7140 slices zero-dimensional,
28 rational core points,
12 square-M lifts,
4 points on both required p=5 and p=7 boundaries,
4 rational W-roots, none a rational square,
0 verified candidates and 0 hits.
```

Moreover, the four boundary points are only sign variants of two singular
curves:

```text
(a,b)=(5,-5/2):    f=x*(x+2)^2*(x^2-11*x/4+4),
(a,b)=(-5/2,5):    f=x*(x+1/2)^2*(x^2-11*x/16+1/4).
```

Thus there is no smooth rational `T0`-halved core point through height `8`.
See `data/contact6_m612_t0_core_slice_h8_summary.txt` and
`code/contact6_m612_h8_core_audit.m`.

### Halving `TB=3D`

Set `t=x-1` and

```text
A = b+3,
C0 = a+b+2,
B0 = A+C0-2.
```

A half of the marked class `D`, equivalently a half of `TB`, is represented
by

```text
uD   = t^2+p*t+q,
ellD = r*t^2+s*t-C0,
ellD^2-f = -2*A*t*uD^2.
```

The exact coefficient equations are

```text
A*q^2 - C0*(B0+s) = 0,
-2*A*C0 + 4*A*p*q - B0^2 - 2*C0*r + s^2 = 0,
-2*A*B0 + 2*A*p^2 + 4*A*q - 2*C0 + 2*r*s = 0,
-A^2 + 4*A*p - 2*B0 + r^2 = 0.
```

These are already generated by `code/contact6_m36_halveD_symbolic.m`.  They
should be fiber-producted with the three cubic-contact equations above, not
replaced by the generic `M(2,12)` chart.

There is also a substantially smaller exact chart for the same `TB` cover.
Put

```text
b=2*s^2-3,       K=m^2,
A3=a-3+4*s^2*r-2*K.
```

Then the open halving locus is cut out by

```text
H1 = s*((a-3)*r^2+4*r-K*(a+3))-r*A3 = 0,
H2 = 8*s^2*(2*s^2*r^2+2*(a-3)*r+2-K*(2*s^2-6))
     -A3^2-32*s^3*r = 0.
```

The first equation is linear in `a`.  Eliminating `a` and stripping the open
factor `s^2*K` leaves one irreducible polynomial `H(s,r,K)` of total degree
`8`, with `31` terms and degree only `3` in `K`.  A half is recovered as

```text
G=x^2+A3/(4*s^2)*x+r/s,
Mtool=2*s^2/m,       Ntool=2*s^2*r/m.
```

This derivation is checked coefficient-by-coefficient in
`code/contact6_m612_tb_lowchart.m`; the search is
`code/contact6_m612_tb_cover_search.m`.

The height-`30` run checked `1,232,100` rational `(s,r)` pairs.  The residue
masks left `6,315` pairs; the cubic had `1,706` rational roots, of which
`146` were rational squares.  These yielded `81` smooth curves, but no curve
was locally compatible with `[6,12]` at every good prime through `31`.
The stronger fixed-curve computation then solved the full three-equation
cubic-contact ideal in `(M,U,v)` for every one of the `81` curves.  All `81`
saturated ideals were zero-dimensional and none had a rational point.  Thus
the complete height-`30` TB stream has no rational independent `3`-direction,
not merely no local survivor.  At height `10`, all seven smooth candidates
were independently checked in the Jacobian and the reconstructed class
really was a half of `TB`; their rejection is therefore an independent-3
obstruction, not a defect in the halving formulas.  See
`data/contact6_m612_tb_cover_h10_p31_summary.txt` and
`data/contact6_m612_tb_cover_h30_p31_summary.txt`, and the stronger exact
summary `data/contact6_m612_tb_core_h30_summary.txt`.

The discovery-oriented fiber product is also computationally tractable.
Fixing `(s,v)` and solving the saturated ideal in `(r,K,M,U)` leaves a
zero-dimensional system.  Five smoke slices, including `(1,-1)`, `(1,-2)`,
`(2,-1)`, `(1/2,-1)`, and `(1,1/2)`, all completed in `0.11--0.52` seconds;
none had a rational point.  The driver is
`code/contact6_m612_tb_core_sv_slice.m`.

The complete height-`4` driver
`code/contact6_m612_tb_core_sv_search.m` then checked all `484` nonzero
`(s,v)` slices.  Of these, `482` saturated ideals were zero-dimensional and
two (`s=+/-1, v=1`) were empty; there were no positive-dimensional or
unresolved exceptional slices.  The search found `11` rational fiber-product
points.  Seven had `M=L^2`, but none had `K=m^2`, so no point reached the
Jacobian verification stage and there was no hit.  Runtime was `395.47`
seconds.  An unfiltered height-`6` run is projected near `30` minutes and was
not started; the next extension should add projective finite-field masks.
See `data/contact6_m612_tb_core_sv_h4_summary.txt`.

#### The recovery boundary is a genuine new one-parameter family

The denominator used to recover `a` is

```text
aden=s*(r^2-K)-r.
```

The simultaneous boundary `aden=anum=0` has two components.  The first is
`r=1/s, K=0`; there `a=4s+3` and `B=2(sx+1)^2`, so the genus-2 curve is
singular.  The second is

```text
r=s*(2s+3),
K=(s+1)^2*(2s-1)*(2s+3).
```

On this component the remaining equation is quadratic in `a`, with

```text
Disc_a(H2)=1024*s^2*(s+1)^5*(s-1/2).
```

Requiring both `K` and this discriminant to be rational squares gives a
genus-zero fiber product.  One convenient two-step parametrization is

```text
s=(t^2-2*t+4)/(4*t),
z^2=t^2+2*t+4,
t=(3-u^2-2*u)/(2*u).
```

It yields two explicit rational `a(u)` branches.  Their genus-2
discriminants are not identically zero.  At `u=2`, both branches are smooth,
have exact torsion `[2,12]`, and are geometrically simple: the two D4
certificates occur at `p=29` and `p=37`, respectively.  Thus this boundary is
a genuine simple order-`12` family, not another split or singular artifact.
It was therefore a natural one-variable place to impose the independent
cubic-contact `3`-direction.
The exact formulas and factorization are in
`code/contact6_m612_tb_recovery_boundary.m` and
`data/contact6_m612_tb_recovery_boundary.txt`; the sample verification is
`code/contact6_m612_tb_recovery_sample_verify.m`.

That one-variable extra-`3` attack is now complete to a substantial bound.
After a simplification, both branches are even rational functions of `u`,
with boundary supported on

```text
u*(u^2-1)*(u^2-9)=0.
```

Exact finite target masks through `p=251`, followed by a conservative
projective CRT enumeration of every primitive `u=n/d` with
`|n|,d <= 100000`, leave only

```text
u=0, 1, 3
```

on each branch.  All three are parameter boundaries.  Every denominator
residue, every bad-reduction residue, and projective infinity were retained,
so this bounded exclusion does not lose a boundary chart.  The formulas are
in `code/contact6_m612_tb_recovery_family.m`, the sieve is
`code/contact6_m612_tb_recovery_height_sieve.c`, and the exact output is
`data/contact6_m612_tb_recovery_sieve_h100000_p251.txt`.

The full `(s,v)` fiber-product search has also been run through slice height
`4`: `484` slices, `482` zero-dimensional and `2` empty, with no unresolved
or positive-dimensional slice.  There were `11` rational points, but none
had square `K`; hence none reached the simultaneous half/core verification.
The summary is `data/contact6_m612_tb_core_sv_h4_summary.txt`.

A first finite-field projectivization now includes the square covers
`K=m^2` and `M=L^2`.  Fiberwise saturation gives very sharp priority masks:

```text
p=5:   3 / 25 fixed (s,v) fibers,
p=7:  10 / 49,
p=11: 25 / 121.
```

Together they retain only `7/2116` height-`6` slices and `311/15876`
height-`10` slices (with every base-denominator residue retained).  This is
not yet a rigorous exclusion mask: before saturation the projective closure
has a point in every fiber, and specialization need not commute with
fiberwise saturation.  Thus every rejected slice is currently classified as
**saturation-boundary-deferred**, not impossible.  The exact sets and the
seven height-`6` priority slices are in
`data/contact6_m612_tb_core_sv_projective_masks_p11.txt`.

The correct strengthening is prepared in
`code/contact6_m612_tb_core_generic_projective_masks.m`: saturate first over
`Q(s,v)`, take the relative projective closure in `(r,m,L,U)`, and only then
specialize its equations modulo small primes.  Coefficient-denominator and
base-infinity residues are retained automatically.  A direct six-variable
universal saturation was abandoned after its preliminary dimension step
used over one gigabyte without finishing.  The smaller generic-fiber version
was then launched successfully, with input equation shapes
`(6,14),(6,19),(12,96),(8,30)`, but its saturation still had not completed
after more than four minutes and about four gigabytes; it was interrupted
when another shared Magma job restarted.  Thus the only fully conservative
projective mask currently remains the raw all-fibers mask.

The seven height-`6` priority slices from `p=5,7,11` have nevertheless been
solved exactly over `Q`.  All seven saturated ideals were zero-dimensional,
and none had a rational point (`8.67` seconds total).  Exact good-affine
enumeration at four further primes gives masks of sizes

```text
p=13: 33/169,   p=17: 12/289,
p=19: 49/361,   p=23: 66/529.
```

These reduce the height-`6` priority count `7 -> 2 -> 0` at `p=13,17`.
At height `10`, the progressive count is

```text
311  --p=13--> 59 --p=17--> 4 --p=19--> 0.
```

So there is no height-`10` slice with good affine reduction in every listed
mask, and hence no open survivor to run through the rational solver.  This is
not a bounded nonexistence result: every failed slice is conservatively
deferred to at least one denominator, infinity, saturation, or bad-reduction
chart.  The exact enumeration and transcript are
`code/contact6_m612_tb_core_affine_masks.py` and
`data/contact6_m612_tb_core_affine_masks_h10_summary.txt`.

### Halving `TC` (and a uniform check for either quadratic class)

For `Q` equal to `B` or `C`, the class `[Q,0]` is divisible by `2` if and only
if there are

```text
uQ = x^2+u*x+v,
ellQ = l3*x^3+l2*x^2+l1*x+l0,
kappa != 0
```

such that

```text
ellQ^2 - f = kappa*Q*uQ^2.
```

Equating degrees `0..6` gives seven equations in the seven displayed
unknowns.  This is a generically finite `16`-sheeted halving cover.  It is the
most direct uniform implementation for `TC`; for `TB`, the four marked-class
equations above are much smaller.

Each of the three geometric halving covers has degree `16` away from the
discriminant.  Consequently their union remains a surface over the
two-dimensional `[6,6]` core; halving is a finite-level condition, not a
codimension-one equation.


## Distinguished Richelot dual

The factorization `f=x*B*C` defines a pointwise-rational maximal isotropic
subgroup

```text
K = {0,T0,TB,TC} in J[2].
```

Put

```text
DeltaR = (a+3)*(b+3)-4 = a*b+3*a+3*b+5.
```

For `DeltaR != 0`, the distinguished Richelot codomain has, up to a rational
square scaling and the standard twist sign, polynomial

```text
gR = DeltaR*R1*R2*R3,
R1 = (b^2-2*a-3)*x^2
     +(2*a*b+6*a+6*b+10)*x +(a^2-2*b-3),
R2 = 2*x^2-(a+3),
R3 = 2-(b+3)*x^2.
```

This follows from the bracket formula

```text
gR = DeltaR*[B,C]*[C,x]*[x,B].
```

The three factors give the pointwise-rational dual kernel and hence
`J'[2](Q)=[2,2]` generically.  The degree-`4` isogeny is an isomorphism on
`3`-torsion, so every rational `[3,3]` on the source survives on `J'`.
Therefore halving exactly one of the three dual factor classes gives the
candidate `2`-primary pattern `[2,2^e]` with `e>=2` on `J'`.  The target is
`[6,12]` only when the exact torsion check gives `e=2`, with no further lift
to order `8` or higher.

There is a useful new necessary condition for the two simple even factors.
Put

```text
A0=a+3,       S=b+3,       Delta=A0*S-4.
```

If the class of `R2=2*x^2-A0` is divisible by `2`, write the exact halving
identity in the form

```text
R2*L^2-Delta*R1*R3 = k*u^2.
```

In the quadratic algebra with `alpha^2=A0/2`, one has

```text
R3(alpha) = -Delta/2,
R1(alpha) = Delta*((S-6)/2+2*alpha).
```

Taking norms of the identity modulo `R2` therefore gives

```text
Norm(-Delta*R1(alpha)*R3(alpha))
  = (Delta^3/4)^2 * ((S-6)^2-8*A0).
```

The left side must be a rational square, so

```text
DC=(b-3)^2-8*(a+3)
```

is a rational square.  Symmetrically, halving `R3` forces

```text
DB=(a-3)^2-8*(b+3)
```

to be a rational square.  These are precisely the discriminants of the
original factors `C` and `B`.  This explains the equal and comparatively
large finite counts for `R2,R3`: the relevant source reduction acquires a
split quadratic factor.  It also gives a cheap characteristic-zero search:
retain all rational cubic-contact core points, including those rejected by
the old `[1,2,2]` filter, test `DB` and `DC` for squares, and exact-test their
Richelot neighbors.  The driver is
`code/contact6_m612_dual_split_core_search.m`.

The same norm argument for the mixed factor gives one more cheap condition.
Indeed

```text
Res(R1,R2)=Delta^2*DC,
Res(R1,R3)=Delta^2*DB,
```

so halving `R1` forces `DB*DC` to be a rational square.  Unlike the two
individual square conditions, this can hold while both `B` and `C` remain
irreducible (their discriminants have the same nonsquare class).  The same
driver therefore retains the union of the three covers
`DB=square`, `DC=square`, and `DB*DC=square` before doing any Richelot or
torsion computation.

The displayed norm proof for `R1` is on the open set
`b^2-2*a-3 != 0`, where `R1` is quadratic.  On the degree-drop locus
`b^2-2*a-3=0`, it becomes a linear factor.  The odd-degree linear-class
criterion gives the same necessary square class: writing its root as
`r=(3-b)/4`, halving forces `-3*(b-5)/(b+3)` to be a square, while
`DB*DC` differs from this by the square `(b+3)^2`.  Thus the prefilter also
remains valid on the smooth linear boundary `b != -3`.

When `B` or `C` splits, the search driver deliberately enumerates every
rational Richelot neighbor rather than only the distinguished `x|B|C`
quotient.  Any `[6,12]` hit would still be valid; attribution to a particular
dual class requires the separate exact `R1/R2/R3` divisibility test.

The exact run is complete through core slice height `10`, but a later audit
found that the original raw counters were dominated by the degree-drop fiber
`b=-3` and singular curves.  After skipping `b=-3`, classifying empty ideals
correctly, and verifying the cubic-contact class before interpreting the
square covers, the corrected results are:

```text
height 6:  formal lifts 7,  verified independent cores 2,
           verified cores on a dual square cover 0;
height 10: formal lifts 16, verified independent cores 6,
           verified cores on a dual square cover 0.
```

The six height-10 verified presentations are two presentations on each of
the three known `DeltaR=0` sources; none has `DB`, `DC`, or `DB*DC` square.
Thus no point reaches the Richelot-torsion stage and the negative open
conclusion is unchanged, but the earlier explanation in terms of dependent
order-3 directions was incorrect.  See
`notes/contact6_m612_core_audit_2026_07_10.md` and the reusable diagnostic
`code/contact6_m612_dual_split_core_audit.m`.

For the known simple seed, Magma finds only this one rational Richelot
codomain.  Its quadratic factors can be normalized as

```text
x^2 - (1493/604)*x - 3691/3624,
x^2 - 125/39,
x^2 - 13/16,
```

and its exact rational torsion is again `[6,6]`, not `[6,12]`.  It is
geometrically simple because it is isogenous to the simple source.

The three smaller core hits

```text
(a,b)=(-19/9,3/2), (-43/25,1/8), (-15/8,5/9)
```

all satisfy `DeltaR=0`; their degenerate Richelot codomains are products of
elliptic curves.  This explains uniformly why those three `[6,6]` examples
failed the simplicity screen.


## Finite local diagnostic

The script `code/contact6_m612_richelot_finite.m` enumerates the affine
`(a,b)` chart.  It requires the source finite group to contain `[6,6]`, builds
the distinguished nondegenerate dual, and tests source and dual compatibility
with `[6,12]`.

```text
p    source66   source612   dual612
5        1          0           0
7        3          0           0
11      14          3           6
13      33         15          16
17      23         18          11
19      69         31          31
23      52         26          24
29      81         53          43
31     246        138         155
```

Thus the construction is not generically locally obstructed: from `p=11`
onward there are many source and dual target residues.  But there are no
good-affine target residues at `p=5` or `p=7`.  Any rational solution in this
normalized chart must therefore have a parameter pole or discriminant
boundary at both primes.  This makes an unrestricted height enlargement a
poor search; the `5`- and `7`-adic boundary strata should be parameterized or
blown up and then combined by CRT.

The stronger class-by-class test is in
`code/contact6_m612_class_finite.m`.  Its mask records which of the three
distinguished classes is actually divisible by `2`; this removes finite
groups whose order-`4` direction comes only from extra splitting after
reduction.  The single-class counts are:

```text
       source single halves        dual single halves
p      T0       TB       TC         R1       R2       R3
13      2        4        4          0        5        5
17      7        1        1          0        3        3
19     11        4        4          2       11       11
23      3        4        4          1        5        5
29     16        7        7          0       11       11
31     27       22       22          4       48       48
```

At `p=5` and `p=7` every distinguished-class mask is zero.  The data suggest
testing the cheap source `T0` cover first.  On the dual, the two simple even
factors `R2` and `R3` are substantially more promising than the mixed factor
`R1`.

For integral affine parameters, the discriminant boundary that must absorb
the reductions at `5` and `7` is supported on

```text
b+3,
a+3,
DB = (a-3)^2-8*(b+3),
DC = (b-3)^2-8*(a+3),
RR = ((b+3)*(a+3)-4)^2
     -(b^2-2*a-3)*(a^2-2*b-3),
```

together with the parameter-denominator/infinity charts.  These are,
respectively, degree loss/collision, the two quadratic discriminants, and the
remaining cross-factor collision.  A rational target with affine-integral
parameters must lie on one of these components modulo each of `5` and `7`.


## Why the generic `M(2,12)` chart is secondary

The previously used halved-contact model

```text
m = (1-z^2)/(4*(r+1)),
T = m*x^2-x+r,
W = (x-r)^2*(T+1)^2+4*m*x^2*T*(T+1)
```

does carry a point of order `12`, but generically its rational torsion is
only `[12]`.  For example `(z,r)=(-8/3,5)` gives factor pattern `[1,1,4]` and
exact torsion `[12]`.  The two rational roots give only one nonzero rational
`2`-class in the even sextic model.

To reach `[6,12]` on this chart one must impose both:

1. an independent rational `3`-torsion direction;
2. an additional rational `2`-class, for example by splitting the residual
   quartic as `2+2` (or `1+3`).

The existing height-40 `[3,12]` search already shows that the first extra
condition is forced to the `p=5` boundary.  All `45` simple-certified exact
boundary survivors had torsion `[12]`; the three exact `[3,12]` points were
split/nonsimple and lay on `Rinf+Z0`.  Since `[6,12]` is a still smaller
subcover, this is computationally less attractive than the contact-6 core,
where `[3,3]` and the extra rational `2` are present from the outset and only
one halving condition remains.

The strongest one-variable component `z=+/-r`, equivalently
`a=(1-r)/4`, has now been closed much farther.  Exact cubic-contact masks
through `p=89` and a projective CRT scan of every primitive

```text
r=n/d,    |n|,d <= 100000,
```

leave only `r=0,1,-1`, all degenerate.  Bad affine reductions and every
projective-infinity residue were retained.  See
`notes/m612_m212_line.md` and
`data/m612_m212_line_height_h100000_summary.txt`.


## Recommended search

### Focused priority after the boundary audit

For the next run, the rational endpoint `P8` family supersedes a blind
two-dimensional enlargement of the core search.  Its exact height-`30`
intersection with the cubic-contact core tested `858` distinct nonzero
`e`-fibers and found no rational contact-open quotient point, with no solver
or exceptional-fiber failures.  The follow-up exact calculation has now
isolated the relevant component: the orthogonal support factor is irreducible
of degree `12`, its recovered `M` is nonsquare, and the signed cover has an
irreducible degree-`24` primitive polynomial in `L`.  A height-`1000` modular
sieve tested `1,216,767` reduced parameters and left only the boundary values
`u=-3,2` with `e=0`.

The forced local disks at `p=7` and `p=17` have also been analyzed.  Both
primes have smooth signed orthogonal branches in some disks, so neither is a
global local obstruction.  Negative statements in the remaining disks are
chart-bounded: higher endpoint weights and, at `17`, nonintegral common-root
blowups remain unresolved.  The next useful work is therefore:

1. attack rational points, the Jacobian, and `S3` quotients of the exact
   genus-`10` signed curve over the `e`-line before the larger P8 pullback;
2. develop global rational-point methods for the exact genus-`145` P8 signed
   cover after exploiting the lower-genus curve and quotients over `e`;
3. analyze the forced disks at `p=19,23,41`, rather than extending the
   undirected height scan;
4. treat the affine `s4=0` and remaining lower-rank weighted charts only as a
   secondary lane;
5. do not extend the one-parameter HLP `G_A` line, whose contact cover has
   genus `51`.

The exact parameterization, the simple `[2,12]` near miss, the obstructed
genus-one component, the exact relative-`3` cover, the bounded sieve, and the
local-disk status are collected in
`notes/contact6_m612_focused_audit_2026_07_10.md` and
`notes/contact6_m612_p8_relative3_status_2026_07_11.md`.  The exact pre-P8
support/signed genera are `5/10`, and the exact P8 support/signed genera are
`73/145`; good reductions at both `p=7` and `p=13` independently reproduce
the latter values.

The general checklist below remains useful when rational core points arise
from any lane.

The next search should modify the existing core-slice solver rather than scan
the generic `M(2,12)` chart:

1. Enumerate rational `(b,v)` slices and solve the zero-dimensional core
   equations for `(a,M,U)`, retaining `M=L^2`.
2. Before any simplicity or full-torsion computation, require bad/boundary
   reduction at both `5` and `7`, using `b+3`, `a+3`, `DB`, `DC`, `RR`, and
   the denominator-infinity charts.  At the other primes apply the stronger
   distinguished-class masks, not only the aggregate `[6,12]` group mask.
3. On each exact core point, test the cheap `T0` square condition first, then
   the full four tangent equations; test `TB` with the marked four-equation
   cover; test `TC` with the uniform quadratic-factor identity.
4. Before constructing a dual, test the cheap discriminant covers: `R2`
   can halve only if `DC` is a square, and `R3` only if `DB` is a square.
   Do this on every rational core point, including sources whose factor type
   is no longer `[1,2,2]`; those split-source points are exactly the intended
   inputs here.  Construct their Richelot neighbors and exact-test torsion.
   The direct exact tests for all three classes, including the mixed class
   `R1`, are now in `code/contact6_m612_dual_class_exact.m`; the universal
   open-chart equations and their size audit are in
   `code/contact6_m612_dual_halving_equations.m`.
5. Only after one class halves, compute `TorsionSubgroup`; require exact
   invariants `[6,12]`.  A source simplicity certificate also certifies every
   Richelot codomain.
6. In parallel, derive the `p=5` and `p=7` projective boundary charts of the
   core equations.  The finite table says this boundary calculation is not
   optional: every target point is forced there in the present affine chart.

Reproducibility files added for this audit are

```text
code/contact6_m612_class_audit.m
code/contact6_m612_class_finite.m
code/contact6_m612_richelot_known_hits.m
code/contact6_m612_richelot_seed_test.m
code/contact6_m612_richelot_finite.m
code/contact6_m612_dual_discriminant_covers.m
code/contact6_m612_dual_split_core_search.m
code/contact6_m612_t0_reduced.m
code/contact6_m612_t0_core_pullback.m
code/contact6_m612_tb_core_tools.m
code/contact6_m612_tb_core_verifier.m
code/contact6_m612_tb_core_search.m
code/contact6_m612_tb_core_sv_slice.m
code/contact6_m612_tb_core_sv_search.m
code/contact6_m612_tb_core_sv_projective_masks.m
code/contact6_m612_tb_core_generic_projective_masks.m
code/contact6_m612_tb_core_affine_masks.py
data/contact6_m36_core_slice_h20_simple_summary.txt
data/contact6_m612_t0_core_slice_h8_summary.txt
data/contact6_m612_tb_core_h30_summary.txt
data/contact6_m612_tb_core_sv_smoke.txt
data/contact6_m612_tb_core_sv_h4_summary.txt
data/contact6_m612_tb_core_sv_projective_masks_p11.txt
data/contact6_m612_tb_core_affine_masks_h10_summary.txt
data/contact6_m612_dual_split_core_h10_summary.txt
```
