# High-torsion bounded sprint: handoff (2026-07-13)

## Outcome and resource policy

This sprint produced no new high-torsion genus-2 Jacobian, but it made four
candidate lanes substantially more precise.  The computations were organized
with at most one substantive process per lane.  Every substantive process had
a hard address-space cap of at most `8 GiB` and a wall-time cap of at most five
minutes; full trigonal replays used at most `4 GiB`.  Smaller validation and
finite-field jobs used tighter caps where possible.

All negative conclusions below are bounded or restricted to explicitly stated
opens.  None is a global nonexistence theorem.

## 1. Cyclic order `44` from the marked order-`22` families

No order-`44` example was found.  The exact marked-class cover has open local
points at several good primes, but every good open residue is killed at both
`p=5` and `p=17`.  A rational half must therefore enter degenerate boundary
disks at both primes.  The cancelled boundary member is a genuine smooth curve
and its marked class is not divisible by `2`; the other three regular charts
specialize to square, non-genus-2 fibers and are only formally live through
first order.  Those formal seeds are not smooth Jacobian or Hensel points.

The conservative exact search through finite parameter height `1000` tested
`1,216,767` reduced rational parameters and found no marked half after `278`
smooth exact divisibility tests.  The parameter at infinity was audited
separately as a degenerate square fiber.

Artifacts:

- [Detailed note](order44_from_order22_scout_2026_07_13.md)
- [Finite-field scout](../code/order44_from_order22_finite_scout.m)
- [Boundary audit](../code/order44_from_order22_boundary_audit.m)
- [Exact halving equations](../code/order44_from_order22_halving_cover.m)
- [Rational sieve](../code/order44_from_order22_rational_sieve.m)
- [Family-equivalence check](../code/order44_from_order22_family_equivalence.m)
- [Recorded output](../data/order44_from_order22_scout_2026_07_13.txt)

## 2. Cyclic order `60` on the `M(12)+5` linear-`B` cover

No order-`60` example was found.  On the explicitly squarefree/coprime open of
the generic linear-`B` quotient, the exact finite-field census at
`p=7,11,13,17,19` found `1,2,3,5,10` open quotient points and
`2,2,2,4,10` signed points.  Every one of the `21` quotient points has
equation-Jacobian rank `4`, hence is a smooth point of local dimension one.
The signed cover is locally alive at every tested prime.

The raw global saturation remained unsuitable under the cap.  A fixed
modulo-`7` base fiber decomposes to one reduced linear point, but this is not a
characteristic-zero component.  The separately treated constant-`B` open has
no point at `p=7,13` and one reduced isolated quotient point at `p=11`, so it
was deprioritized.

The census deliberately excludes repeated-support and shared-Weierstrass
boundaries.  Those can encode valid divisors and remain separate undecided
tasks; the generic-open counts do not rule them out.

Artifacts:

- [Detailed note](m12_general5_b2zero_modular_triage_2026_07_13.md)
- [Finite-field enumerator](../code/m12_general5_b2zero_modular_points.py)
- [Bounded Magma triage](../code/m12_general5_b2zero_modular_triage.m)
- [Root-quotient derivation and boundary caveats](m12_general5_b2zero_rootquotient.md)
- [Recorded output](../data/m12_general5_b2zero_modular_triage_2026_07_13.log)

## 3. The contact-`30` C3-root route to `[2,60]`

No `[2,60]` candidate was found.  The rational-root condition is an
irreducible trigonal cover of the original parameter line with genus `12` and
different degree `28`.  It has a rational involution whose exact irreducible
trigonal quotient has genus `6` and different degree `16`.  The quotient is a
substantial reduction from the previous total-degree-`32` plane projection.

Exact searches through height `5000` on both the original `R`-line and the
quotient `v`-line tested `30,401,831` reduced finite rational parameters in
each box.  The original search found six rational-root points and the quotient
search found five; all original points are contact boundary, and every
quotient point either lifts only to boundary or has nonsquare lift
discriminant.  The infinity cubics have no rational root.

These searches do not determine the rational points of either the genus-`12`
curve or its genus-`6` quotient.  The rejected additional transformation is
only the displayed base-Mobius candidate, not a proof that no other geometric
map exists.

Artifacts:

- [Detailed note](contact30_c3root_trigonal_geometry_2026_07_13.md)
- [Genus and discriminant geometry](../code/contact30_c3root_trigonal_bounded_audit.m)
- [Quotient and bounded searches](../code/contact30_c3root_trigonal_quotient_search.m)
- [Recorded output](../data/contact30_c3root_trigonal_geometry_2026_07_13.txt)

## 4. Clebsch--Klein orbit-`12` source halving toward order `160`

No order-`160` cover point was found in the selected CRT coset.  After an
integral-chart correction, one product of certified open local disks modulo
`11^3`, `19^2`, and `23^2` becomes two numerator-denominator congruence
lattices modulo `254,179,739`.

At Euclidean radius `10^6`, the exact search combined `3,152` short `t`
vectors with `3,138` short `m` vectors.  All `9,890,976` resulting labelled
chart points are smooth rational Clebsch--Klein points, and none has even one
of the four required exact rational-square radicands.

This is a result for one selected product of three Hensel disks only.  It does
not enumerate the many other compatible local states and is not a global
obstruction to the order-`160` cover.

Artifacts:

- [Detailed note](elkies22210_orbit12_ck_crt_lattice_2026_07_13.md)
- [Exact CRT/lattice search](../code/elkies22210_orbit12_ck_crt_lattice.py)
- [Recorded radius-`10^6` output](../data/elkies22210_orbit12_ck_crt_lattice_r1000000.txt)
- [Prior local-disk certification](elkies22210_source_halving_local_and_search_2026_07_11.md)

## Ranked next steps

1. **Generic `M(12)+5` open:** replace monolithic saturation by structural
   projection, function-field elimination, or modular interpolation one
   equation at a time.  This remains the strongest live route to an exact
   cyclic order-`60` example.
2. **Genus-`6` trigonal quotient:** compute its Jacobian rank and decomposition
   and search for maps to lower-genus curves; if the rank permits, apply
   Chabauty or a Mordell--Weil sieve.
3. **Clebsch--Klein order `160`:** organize all compatible local disks into a
   finite collection of integral chart cosets, or derive global component
   geometry, before any further radius search.
4. **Order `44`:** revisit only through weighted higher-order saturation and
   normalization of the three square boundary disks, not a larger undirected
   height box.

The constant-`B` `M(12)+5` open is lower priority than these four directions;
its excluded collision boundaries should be treated as their own
divisor-arithmetic problem rather than inferred from the open census.
