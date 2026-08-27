# Top-ten frontier smoke tests: `[6,12]` and `[2,2,4,8]`

Date: 2026-07-18

These are bounded reproducibility tests for two targets in the pre-test
top-ten queue.  They are not nonexistence searches.

## `[6,12]`

The Prym dependency in the `E8 -> E4` route was rerun with

```text
magma -b code/contact6_m612_prym_rank_verifier.m
```

It returned

```text
RankBounds(J(D)) = 1..1
torsion invariants = []
generator (Mumford): (x^2 + 2*x + 4, 5*x + 5)
PRYM_RANK_VERIFIED
```

for

```text
D: y^2 = -3*x^6 + 24*x^3 - 75.
```

Thus the advertised Prym has rank one and trivial rational torsion, so the
genus-four cover remains in Prym-Chabauty range once the exact quotient and
bigonal maps are certified.

The near-miss verifier was rerun with

```text
magma -b code/contact6_m612_nearmiss_recertify.m
```

It returned exact torsion `[2,2,6]`, not the `[2,12]` stated in the script
header and some prose, together with a valid root-power witness

```text
p = 7
chi = T^4 + T^3 + 5*T^2 + 7*T + 49
deg minpoly(pi^n) = 4, 2 <= n <= 12.
```

The simplicity evidence is sound, but the near-miss description should be
corrected before it is used as structural evidence for `[6,12]`.

The next high-leverage step is not another height box: certify the exact
characteristic-zero degree-24-to-`E8` quotient and the bigonal divisor
transport, then implement the precision-safe Coleman plus Abel-Prym
Mordell-Weil sieve described in
`notes/m612_review_and_top3_plan_2026_07_13.md`.

## `[2,2,4,8]`

The HLP normalization test was rerun with

```text
magma -b code/m2248_hpl_normalization_check.m
```

Both known split HLP examples are recognized by the full cover equations.
For each example the square interpretation gives six normalizations, four
full sources, and 32 full witnesses; the nonsquare interpretation gives none.
This is a useful positive control for the equations, but all of these points
are geometrically split.

Two bounded searches were then rerun:

```text
magma -b height:=6 max_hits:=3 code/m2248_rst_direct_search.m
magma -b height:=8 max_hits:=3 code/m18_m14_full_split_search.m
magma -b height:=8 max_hits:=3 code/m18_m2228_curve_extra_halving_search.m
```

Their summaries were

```text
rho-sigma-tau: checked 3908, cover_square 2, cd_square 0, hits 0
full-split tangent: checked 7224, split_quadratics 0, hits 0
known [2,2,2,8] curve: checked 84, exact base 84, extra-half hits 0.
```

The current HLP-centered chart is therefore a verified negative control, not
a discovery lane.  The best next move is the simultaneous pair-fiber
geometry: compute ranks and rational points on the three genus-one fibers
through the normalized HLP point, or move to the Richelot graph of a simple
`[2,4,8]` seed and impose the extra rational 2-direction there.

## Promoted target `[8,8]`

After the first three entries in the pre-test queue were realized, `[8,8]`
was promoted into the revised open top ten.  The reduced-cover positive
control was rerun with

```text
magma -b primes:="7,11" samples:=2 \
    code/m18_m14_88_reduced_finite.m
```

The reduced equations agree exactly with the intrinsic condition
`T_x in 4*J(F_p)`:

```text
p=7:  good bases 16, intrinsic 8, reduced 8, mismatches 0
p=11: good bases 56, intrinsic 24, reduced 24, mismatches 0.
```

A fresh rational smoke test

```text
magma -b height:=8 max_hits:=2 \
    code/m18_m14_88_reduced_search.m
```

gave

```text
checked 7569, open_good 22, first_bases 6,
first_solutions 24, second_solutions 0, hits 0.
```

Thus `[8,8]` remains locally plausible, but a larger rational height box is
not the next move: eliminate the reduced second-cover variables or determine
the rational points on its genus-one pair fibers.
