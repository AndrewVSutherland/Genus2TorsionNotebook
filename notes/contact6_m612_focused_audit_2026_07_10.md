# Contact-6 `[6,12]`: focused audit and next work

Date: 2026-07-10.
Updated: 2026-07-11 with the exact relative-`3` cover, sieve, and local disks.

## Bottom line

No geometrically simple `[6,12]` curve was found.  The focused run did,
however, turn the previously broad `p=5` boundary problem into a short list
of exact conclusions.

The most useful surviving object is the rational `P8` component of the
endpoint `R3`-halving cover.  It gives a one-parameter family of genuine
halves and the explicit simple near miss

```text
(a,b)=(-409/200,0),
source torsion = [2,2,6],
dual torsion   = [2,12],
halves(R1,R2,R3) = [false,false,true].
```

It misses `[6,12]` only because it has one rational `3`-direction instead of
two.  Intersecting this family with the independent cubic-contact `[3,3]`
core was therefore the best next calculation.  That intersection is now an
explicit characteristic-zero cover: an irreducible degree-`12` support
factor and its connected irreducible degree-`24` signed lift.  A height-`1000`
residue sieve has no nonboundary survivor, while local analyses at `7` and
`17` find genuine target branches in some forced disks.  Thus the lane is
sharply defined but neither solved nor locally obstructed.


## Audited decision table

| lane | exact result | decision |
|:--|:--|:--|
| corrected affine core search | height `6`: 2 verified cores; height `10`: 6; no split-cover point or hit | height alone is not the issue; work on boundary charts |
| source `T0` at `p=5` | the unique escaping rank-3 cone dies after the second deflation; compensated pole charts stay on `U+2v=0` through mod `125` | deprioritize; lower-rank weighted poles remain formally unresolved |
| smooth affine `DB/DC` cones | 16 cones survive the necessary square covers, but zero survive the exact `R3/R2` halving equations on `s4 != 0` | closed on the valid square-quartic chart; `s4=0` remains a degree-drop boundary |
| weighted endpoint `E9+R3` | four smooth `Z_5` core disks survive and satisfy exact `R3` locally | local lane is real, so globalize its component |
| `E9` component `P16` | normalization is `Y^2=2*Z^4-6*X^4`; no `Q_2` or `Q_3` point, but a `Q_5` point | globally closed |
| rational component `P8` | explicit rational parametrization; genuine simple `[2,12]` dual near miss | highest-priority exact family |
| `P8` relative-3 cover | exact irreducible degree-`12` support and degree-`24` signed cover; height-`1000` sieve has no nonboundary survivor; target local branches at `7` and `17` | pursue global geometry and remaining forced primes, not another blind scan |
| HLP deformation line `G_A` | irreducible degree-40 contact cover of genus `51` | do not search this line further; only a better two-parameter slice is competitive |


## Exact weighted geometry

After saturation by `s4`, the endpoint `R3` equations have two prime curve
components.

The component carrying the `E9` residue disks is `P16`.  Resolving a
degree-two place and applying exact Riemann--Roch gives the minimized model

```text
Y^2 = 2*Z^4 - 6*X^4,
Jacobian: y^2 = x^3 + 3*x.
```

A primitive congruence calculation gives no point modulo `16` and no point
modulo `9`; Magma independently reports no `Q_2` or `Q_3` point.  This is a
global obstruction to every certified `E9` `Q_5` disk, not merely a failed
rational search.

The other component, `P8`, normalizes to

```text
y^2 = 100-6*t^2.
```

Writing

```text
t(u) = 4*(u^2+u-6)/(u^2+6),
y(u) = 2*(u^2-24*u-6)/(u^2+6)
```

gives exact rational functions `(e(u),mu(u),nu(u))` satisfying the full
`R3` covariants.  At `u=12`,

```text
(e,mu,nu)=(-200/409,-36320/167281,38136/167281).
```

The exact torsion audit at `a=1/e`, `b=0` gives the near miss displayed in
the bottom line.  Both source and dual are geometrically simple, certified
at `p=23` by the irreducible Frobenius polynomial

```text
529*T^4 + 22*T^2 + 1.
```


## Bounded `P8`-core intersection

The exact `P8` parameter was scanned through rational height `30`.  Duplicate
`e`-fibers were removed before solving the cubic-contact equations.  Every
fixed fiber was saturated by

```text
M*v*(U^2-4*v^2),
```

then rational quotient points were required to have `M=L^2`, and every
survivor would have been checked by exact Jacobian arithmetic.

The result was

```text
rational u parameters                 1111
distinct nonzero e fibers              858
duplicate e fibers                     251
e=0 boundary parameters                  2
parameter poles                           0
rational contact-open quotient points    0
square-M points                           0
verified independent [3,3] cores          0
saturation / dimension / Variety failures 0
hits                                      0
```

The projective parameter `u=infinity` was audited and duplicates the
`e=-200/409` fiber on the other conic branch.  At that fiber the unsaturated
core contains only the universal degenerate section

```text
(M,U,v)=(1,-2,1),    U^2-4*v^2=0,
```

so it is correctly removed from the contact-open core.


## Exact relative cover and subsequent sieve

The generic degree-`40` support quotient is now resolved as

```text
1 + 12 + 27.
```

The degree-`1` factor is the repeated marked support, the degree-`12` factor
is the orthogonal support relevant to a second rational `3`-class, and the
degree-`27` factor is nonorthogonal.  Exact characteristic-zero reconstruction
gives an irreducible `Phi12(v)`, degree-`11` recovery maps, and

```text
Phi24(L)=Res_v(Phi12(v),L^2-M(v)).
```

Here `M` is nonsquare and `Phi24` is irreducible of degree `24`, with `L`
primitive.  The exact pre-P8 support and signed fields have genera `5` and
`10`: `M` has exactly two odd degree-`1` places.  After P8 base change, the
exact relative different pulls back with degree `168`, giving support genus
`73`; the connected sign lift is unramified, giving signed genus `145`.
Good reductions `7` and `13` independently reproduce `73/145`.

The height-`1000` modular sieve tests `1,216,767` reduced parameters and
leaves only `u=-3,2`, both on `e=0`.  This is a rigorous bounded exclusion,
not a global no-point theorem.

At `p=7`, the two parameter-pole disks have smooth signed orthogonal branches;
at `p=17`, the four `DC=0` disks do as well.  Hence neither prime obstructs
the route.  Negative results in the other disks retain their bounded-chart
caveats: higher endpoint weights at `7` and `17`, and nonintegral common-root
blowups at `17`, are not exhausted.  The consolidated current status is in
`notes/contact6_m612_p8_relative3_status_2026_07_11.md`.


## Recommended next work

1. Use the exact genus-`10` signed curve over the `e`-line as the first global
   target: determine its rational points/Jacobian and exploit its `S3`
   quotients before the larger P8 pullback.
2. Develop global rational-point methods for the exact genus-`145` P8 signed
   cover after exploiting the lower-genus curve and quotients over `e`.
3. Analyze the forced bad-reduction disks at `p=19,23,41`.  The completed
   `p=7` and `p=17` analyses contain positive target branches and cannot close
   the route alone.  Do not replace this with another undirected height scan.
4. Revisit the affine `s4=0` and lower-rank weighted strata only as a
   secondary task.  They are unresolved degree-drop charts, whereas `P8`
   already supplies exact rational halving geometry.
5. Do not extend the one-parameter HLP `G_A` line.  Its genus-`51` contact
   cover makes a carefully chosen two-parameter surface the only sensible
   HLP continuation.


## Reproduction

The main exact commands are

```text
magma -b code/contact6_m612_affine_dual_exact_leading_mod5.m
magma -b code/contact6_m612_t0_local5_second_blowup.m
magma -b code/contact6_m612_weighted_R3_p16_independent.m
magma -b code/contact6_m612_weighted_R3_geometry_p8_param.m
magma -b code/contact6_m612_weighted_R3_geometry_p8_audit.m
magma -b height:=30 code/contact6_m612_weighted_E9_P8_core_intersection.m
magma -b code/contact6_m612_relative3_exact_reconstruct.m
magma -b mode:=search height:=1000 prime_bound:=71 max_exact:=500 \
  code/contact6_m612_p8_extra3_residue_sieve.m
magma -b code/contact6_m612_p8_p7_bad_disks.m
magma -b code/contact6_m612_p8_p17_bad_disks.m
magma -b code/m612_hlp_GA_monodromy_geometry.m
```

Detailed derivations are in

```text
notes/contact6_m612_core_audit_2026_07_10.md
notes/contact6_m612_boundary_charts_2026_07_10.md
notes/contact6_m612_affine_dual_exact_2026_07_10.md
notes/contact6_m612_t0_local5_2026_07_10.md
notes/contact6_m612_weighted_R3_geometry_2026_07_10.md
notes/contact6_m612_weighted_R3_p16_obstruction_2026_07_10.md
notes/contact6_m612_weighted_E9_P8_core_intersection_2026_07_10.md
notes/contact6_m612_p8_relative3_exact_2026_07_10.md
notes/contact6_m612_p8_relative3_status_2026_07_11.md
notes/m612_hlp_GA_probe_2026_07_10.md
```
