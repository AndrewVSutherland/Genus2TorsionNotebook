# Point of order 60 attempts

The target is an actual rational Jacobian point of order `60`.  The two most
direct constructions tested were:

```text
M(12) + 5-torsion       -> lcm(12,5) = 60,
contact-5 [2,20] + 3   -> lcm(20,3) = 60.
```

## `M(12)` one-parameter line plus `5`

I added

```text
code/m12_plus5_order60_search.m
```

for the one-parameter family `a=(1-r)/4`, which has known torsion containing
`Z/12Z x Z/2Z`.  The script applies the necessary good-reduction condition

```text
5 | #J(F_p)
```

at all good primes `p != 5` up to `prime_bound`, allowing bad or boundary
residues through.  Exact torsion is computed only for residue survivors; if the
torsion exponent is divisible by `60`, the script prints a Mumford
representative of an order-`60` class.

The run

```text
magma -b height:=1000 prime_bound:=251 max_exact:=200 max_print:=50 \
    progress_interval:=200000 code/m12_plus5_order60_search.m
```

gave:

```text
checked           1216767
smooth            0
residue_survivors 0
exact_tests       0
hits              0
```

The first obstruction is already visible at `p=7`:

```text
p=7: good 3, allowed5 0, bad 4
```

Thus the tested one-parameter `M(12)` line gives no order-`60` candidate
through height `1000`; all rational parameters are killed by good-prime
`5`-divisibility residue conditions before exact torsion.

## Contact-5 `[2,20]` plus `3`

The direct order-`60` route is the existing live-boundary search

```text
code/contact5_extra2_plus3_live_boundary_search.m
```

because every smooth point on the extra-2 contact-5 loci has an order-`20`
class, and an independent rational `3`-torsion class would combine with it to
produce order `60`.

The previously recorded all-branch live-boundary run to height `5000` and
`prime_bound=251` found only the singular boundary survivor `z=-2, t=-3`.
I extended the two cheaper linear live branches to height `10000`.

For `linear_tminus3`:

```text
unique_t             4342196
exact_tests          0
plus3_hits           0
checked_crt_loop     6955529
primitive_crt_params 4342196
survivors            1
```

The unique survivor is again the singular point `z=-2, t=-3`, bad at every
filter prime, so it never reaches exact torsion.

For `linear_pole`:

```text
unique_t             4342216
exact_tests          0
plus3_hits           0
checked_crt_loop     6955525
primitive_crt_params 4342216
survivors            0
```

Thus neither linear live branch produces a smooth `[2,20]+3`, hence no
order-`60` point, through height `10000`.

## Current status

No point of order `60` was found.

The remaining plausible computational directions at that stage were:

```text
1. Extend the larger qq_tminus3 live branch beyond the recorded height 5000.
2. Build a two-parameter M(12)+5 sieve instead of staying on a one-parameter line.
3. Do component-wise boundary analysis for the M(12)+5 obstruction at p=7.
```

The `qq_tminus3` branch has since been refined locally at `p=11`.  Its four
live residue disks are `r=3,5,7,8`; the new `[2,2,20]` specialization
`r=-1/3` lies in the double-root disk `r=7`, but has no 3-part.  The
cubic-contact cover is ramified/singular at every one of its 32 boundary
points.  More decisively, the two contact points reducing to

```text
(m,U,V)=(2,-1,0), (-2,-1,0)
```

project onto every lift of all four `r` residues modulo `11^2` and `11^3`.
Thus a deeper `p=11` CRT sieve gives no additional parameter restriction, and
another blind height extension is low-value.  The next meaningful attack on
this branch is normalization of the saturated global cubic-contact
components.  Reproducible diagnostics are
`code/contact5_qq_plus3_p11_etale_diagnostic.m` and
`code/contact5_qq_plus3_p11_lift_projection.py`.

## Contact-(5,6) order-30 plus halving

A third direct construction was tested: halve the distinguished rational
2-class in the exact simultaneous contact-5/contact-6 order-30 family, then
combine its order-4 half with the visible order-15 class.

This route is now ruled out completely by a local obstruction at `p=11`,
including every bad/boundary residue disk.  The exact cover, projective finite
sieve, and boundary analysis are recorded in
`notes/contact30_to_60_halving_obstruction.md`.  In particular, this is not
merely another negative height search.

## Extra-2 Richelot lift

Forcing the cubic factor in the order-30 family to have a rational root gives
a distinguished Richelot neighbor whose finite fibers very frequently have
exponent 60.  This is locally alive, but the pointwise-rational dual Richelot
kernel forces rational 2-rank at least two.  Consequently it targets
`[2,60]`, not exact cyclic `[60]`.  Its compact root-cover geometry is recorded
in `notes/contact30_extra2_richelot_route.md`.

## Full M(12) exact contact-5 cover

The old one-parameter `M(12)+5` line has been replaced by the full
split-Weierstrass `M(12)` surface and, more importantly, by the exact
point-contact-5 equations.  The cover is good-open locally at every tested
prime, so it remains a viable route to exact cyclic `[60]`.  A
positive-controlled height-50 search checked 9,579,025 `(r,z)` pairs and
found no rational cover point.  The smaller normalized `(b,w)` chart and the
full derivation are recorded in `notes/m12_contact5_exact_cover.md`; component
normalization, not a larger `(r,z)` box, is the next step.

## Exact cyclic order-20 plus rational 3-torsion

The earlier extra-2 contact-5 searches were not the correct generic target:
they force rational 2-rank at least two and therefore naturally produce
`[2,60]`, not exact `[60]`.  A distinct one-parameter family has a marked
order-5 class and a marked order-4 class while generically having rational
2-rank one:

\[
 h_t=1+tx+\frac{t^2-1}{2}x^2,
 \qquad f_t=h_t^2-\frac{(t+1)^4}{4}x^5.
\]

The sum of the marked classes has order 20, so an independent rational
3-class would give exact cyclic order 60 without forcing an extra rational
2-direction.  This is therefore a genuine cyclic-60 lane.

Necessary good-reduction masks through `p=113`, implemented in
`code/contact5_order20_plus3_height_sieve.c`, exclude every smooth rational
parameter of height at most `100000` in both `t` and
`s=(t-1)/(t+1)`, after exact testing of the few finite-mask survivors.  Every
survivor has exact torsion `[20]`; independent Frobenius certificates also
show that these negative controls are geometrically simple.

The full cubic-contact equations have now been eliminated exactly.  The first
generic square-cover projection is one irreducible curve of bidegree
`(56,80)`.  A fixed-Weierstrass normalization reduces this to an irreducible
plane curve `P44(s,r)` of bidegree `(22,40)`, together with an explicitly
recovered quadratic sign cover

\[
 A(s,r)j^2+B(s,r)=0.
\]

Direct substitution modulo `P44` verifies the recovery identities.  Several
natural special slices (`V=0`, `E=0`, `U=0`, `A=0`, `B=0`, `M=1`, and `V=1`)
have also been eliminated exactly and contain no rational open candidate.
The low-degree resultant factors and the finite denominator locus of the
fixed-Weierstrass model still require explicit lift checks; this caveat is
important before calling the degree-44 model globally complete.  Full details
are in `notes/contact5_order20_plus3_full.md`.

## General Mumford 5-torsion over M(12)

The point-contact-5 cover is not the full order-5 problem.  A general Mumford
divisor of order 5 is represented by the exact norm identity

\[
 A^2-B^2F=q^5,
\]

with `A` monic of degree 5, `deg(B)<=2`, and `q` monic quadratic.  This gives
a two-dimensional cover of the compact two-parameter `M(12)` chart, rather
than the thin `P55` curve.

A new independently point-counted C sieve tests the necessary condition
`5 | #J(F_p)` through `p=79`, with an asserted order-12 positive control at
every good residue.  Its mask counts through `p=19` agree exactly with a
separate Magma implementation.  The run checks 2,393,618,932 nontrivial
rational `(b,w)` pairs of height at most 200.  Ten finite-mask survivors remain;
fresh exact reconstruction gives torsion `[12]` for every survivor, with D4
geometric-simplicity certificates.  Thus there is no cyclic-60 point in this
box, even allowing arbitrary rational 5-torsion support.

The first global symbolic subcover, `deg(B)<=1`, has five exact residual
equations after recursively eliminating `A`.  Its saturation and factorization
are now the preferred `M(12)+5` geometry computation.  See
`notes/m12_general5_order60.md`.

## External record and an isolated order-30 seed

The current literature/database audit found no geometrically simple cyclic-60
example.  The published HLP order-60 family is geometrically split.  A screen
of Sutherland's 487,493-curve 5-smooth-discriminant list leaves only curves
with decomposable Sato--Tate type after the necessary good-prime tests.

Nicholls' geometrically simple order-30 curve is not a member of the workspace
contact-(5,6) family, but it also cannot be lifted: at the good prime `5`, its
finite Jacobian is cyclic of order 30, so its unique order-2 class is not a
double.  The exact curve audit and carefully scoped primary-source record are
in `notes/z60_literature_and_order30_seed_audit.md`.

## Exact split `[60]` control obtained

The HLP specialization `(t,u,y)=(1/3,-1,9)` has now been reconstructed
explicitly and checked with Magma.  The curve

```text
y^2 = -185*(125*x^2-1728)*(2000*x^4-48525*x^2+294912)
```

has exact rational torsion subgroup `[60]`; an explicit Mumford generator has
order `60`.  The even sextic has the non-hyperelliptic involution `x -> -x`,
so its Jacobian is geometrically split.  Thus this is an exact positive
control, not a solution to the geometrically simple target.  See
`notes/hlp_z60_explicit_control_2026_07_14.md` and
`code/hlp_z60_explicit_verify.m`.

The subsequent split-quadratic and irreducible-quadratic general-`B` searches
found no simple rational example in their stated boxes.  Their scoped results
are in `notes/m12_general5_fullquad_d53_2026_07_14.md`,
`notes/hermite5_order60_modular_2026_07_14.md`, and
`notes/m12_general5_fullquad_irred_2026_07_14.md`.

## Transverse deformation of the exact split `[60]` control

The fixed cyclic generator has simultaneous exact order-5 norm, order-3 contact, and order-4 halving identities. Their normalized auxiliary ranks are `11,7,7`; the full incidence has rank `25` in `32` variables and projects etale onto all seven sextic coefficients. The unique Humbert-4 tangent normal is `(0,81125,0,904800,0,9916416,0)`, so `dF=x` is a genuine transverse marked deformation.

On the simplest line `F_t=F_0+t*x`, complete finite-Jacobian masks give `allowed={0}, bad=empty` at `7,11,13,19,29,73,79,107,131,149`. Their product is `6643194523826861`, which excludes every nonzero primitive `t=n/d` with `max(|n|,d)<=81505794` by a direct divisibility argument. The line does reach the simple locus (`t=1` has absolutely simple reduction at `11`), but that sample fails the order-60 masks.

See `notes/hlp_z60_marked_identities_and_tangent_2026_07_14.md` and `notes/hlp_z60_transverse_slice_search_2026_07_14.md`. This rules out one natural transverse pullback, not all simple cyclic-60 Jacobians.
The tempting exact link `B5=q0` was checked separately and is locally a dead end: its smooth germ is the coordinate saturation of an even-sextic subgerm, hence remains in Humbert 4 to all formal orders. Any further transverse construction must let these quadratics separate; see `notes/hlp_z60_linked_germ_2026_07_14.md`.
