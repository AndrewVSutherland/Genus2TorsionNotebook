# Target `Z/35` (top-10 item 10): full `A_1(5)` threefold — local go/no-go and first CRT-guided search

Goal: a genus-2 curve over `Q` with geometrically simple Jacobian and a rational
torsion point of order `35`.  Cyclic `35` is one of the four open cyclic gaps
`{31, 35, 37, 38}` and the only one with strong positive evidence nearby:
HLP (Forum Math. 12, 2000) realize `C35` on SPLIT genus-2 Jacobians by a
positive-rank elliptic-curve family, and Howe's split order-`70` curve contains
a `35`-torsion point.  No simple example is known.

## Dossier summary

The in-repo M-type (point-contact) charts are locally obstructed and searched
empty (`notes/simple_35_attempt.md`):

```text
contact-7 plus 5: open surface obstructed mod 3; all p=3 boundary branches dead
                  ((1,1) branch has limiting generalized Jacobian of order 14,
                  no 5-part; the branches with a 5-part degenerate the 7-point)
contact-5 plus 7: open surface obstructed mod 3 AND mod 5; searches empty
```

But both charts only cover classes of the form `n(P - infinity)` (torsion
supported at a single curve point).  Elkies (Contemp. Math. 796, 2024, LuCaNT)
gives the FULL rational moduli threefold `A_1(5)` of pairs `(C, T5)`:

```text
universal curve  y^2 + (L'Q' - LQ) y = Q^2 Q',
L = x, L' = 1, Q = q2 x^2 + q1 x + q0, Q' = Q - x,
coordinates (q0,q1,q2); involution T -> 2T is (q0,q1,q2) -> (q2, 1-q1, q0);
the class T = [{Q=0,y=0} - K] has order 5.
```

A `C35` curve needs a rational 5-torsion class AND a rational 7-torsion class;
the locus `A_1(5,7)` is a threefold cover of moduli (codimension 0), so
rational points are not dimension-obstructed — the question is whether the
non-split components have (small-height) rational points.

## Strategy (ranked)

1. **Route A (executed below): sieve the full `A_1(5)` chart for `7 | #J`.**
   On `(q0,q1,q2)` in `Q^3`, the 5-torsion is free; a rational 7-torsion class
   forces `35 | #J(F_p)` at every good `p != 5,7` (prime-to-`p` part at
   `p = 3,5,7`).  First establish the local landscape (does the mod-3
   obstruction of the contact charts persist on the full threefold?), then run
   a CRT-guided height search with conservative bad-reduction passes, then
   exact torsion in Magma for survivors.
2. **Route B: reconstruct the `A_1(7)` threefold ("NDE 2003", cited by Elkies)
   and sieve it for `5 | #J`.**  Mirror of Route A; the chart guarantees the
   rarer 7-structure exactly, and the imposed condition (5-divisibility,
   local pass rate ~0.23) is easier to satisfy, so the rational points of
   `A_1(5,7)` may sit at lower height in the 7-chart coordinates.  Requires
   digging the explicit equations out of Elkies' g2_tors page / NDE 2003.
3. **Route C (fallback, pessimistic): RM slices.**  Elkies' `Delta = 0`
   surface parametrizes RM5 curves with a sqrt5-torsion point; imposing 7
   there is a codim-1 condition on a rational surface.  But all RM abelian
   surfaces over `Q` are GL(2)-type hence modular, and Cowan's `J_0(N)` scan
   (N prime up to 2*10^6) found essentially no 35-torsion — treat this route
   as a long shot only.

## Validation (calibration before anything new)

Script: `code/claude_c35_validate.m` (run `nice -n 15 magma -b claude_c35_validate.m`).

```text
Elkies identity: (q0,q1,q2) in {(1,1,1),(2,-1,3),(1/2,3,-2),(-3,2/5,7)}:
  genus 2, Order(J![Qmonic,0]) = 5 in all four cases.        [PASS]
contact-7 [28] curve y^2 = 4x^5+21x^4-70x^3+79x^2-42x+9:
  torsion [28], marked class J![x-1,-1] has order 7,
  L_5 = 25T^4+2T^2+1 irreducible.                            [PASS]
  CAVEAT: L_5 is biquadratic, so its 12th-power transform is REDUCIBLE:
  p=5 certifies Q-simplicity only, NOT the repo geometric-simplicity
  certificate.  Use a different prime for geometric simplicity here.
Howe split C70 curve (contains C35): 35 | #J(F_q) at all 19 good
  q in [11,97] (bad model at 11,13 skipped conservatively)   [PASS]
  — i.e. a genuine 35-curve sails through the sieve semantics.
```

## Main test 1: the `p=3` landscape (the decisive prime)

Scripts: `code/claude_c35_p3landscape.gp`, `code/claude_c35_f3curves.gp`.

Exhaustive enumeration of ALL genus-2 models `y^2 = f`, `f` squarefree of
degree 5 or 6 over `F_3` (1944 models, 1296 squarefree):

```text
#J(F_3) attains 35 on exactly 12 models  (Weil range: #J <= 55, max attained 36)
all 12 have the SAME Frobenius charpoly
    x^4 + 4x^3 + 9x^2 + 12x + 9 = (x^2 + x + 3)(x^2 + 3x + 3)
    (#J = 5 * 7 = (1+1+3)(1+3+3); split Frobenius, ordinary)
example: y^2 = x^6 + 2x^2 + 1 over F_3.
```

Consequences:

```text
(a) The mod-3 obstruction of the contact charts does NOT persist globally:
    a rational C35 curve MAY have good reduction at 3.
(b) It is extremely rigid: good reduction at 3 forces reduction isogenous to
    E5 x E7 with traces (-1, -3) — a single Weil polynomial.
(c) Any curve with a rational point of order 70 must have bad reduction at 3
    (70 > 55.7 = Weil bound), e.g. Howe's C70 curve (3-adically bad model).
```

On the Elkies chart mod 3 (27 triples): 20 are boundary (non-genus-2), and
ALL 7 smooth residues have `7 !| #J` (orders: `5 | #J` in all seven, none
divisible by 35).  So a rational `C35` point of `A_1(5)` is forced to be
**3-adically boundary in this chart** (nonintegral `q_i` or singular chart
reduction) — but NOT forced to have bad reduction: the 12 good `F_3` curves
above are `F_3`-points of `A_1(5)` lying outside the smooth chart locus.

## Main test 2: local densities on the full threefold

Script: `code/claude_c35_sweep.c` (C, single-threaded; validated against gp
exactly at `p = 3, 5, 7`).  `./c35_sweep density p [nsample]`:

```text
p    total    boundary  smooth   5|J      7|J     35|J    frac35
3    27       20        7        7        0       0       0.000
5    125      56        69       69       8       8       0.116   (7|#J filter)
7    343      128       215      215      32      32      0.149   (5|#J vacuous)
11   1331     298       1033     1033     156     156     0.151
13   2197     470       1727     1727     300     300     0.174
29   24389    2332      22057    22057    3162    3162    0.143
41   68921    4768      64153    64153    10494   10494   0.164
71   60000-sample  2482 57518    57518    8067    8067    0.140
```

Findings: `5 | #J` holds at EVERY smooth chart point at every tested prime
(the universal 5-torsion never dies), so `35 | #J` is exactly `7 | #J`, with
healthy densities `~ 1/7 + O(1/p)` everywhere except the empty `p = 3` chart.
**Local verdict: GO — no local obstruction on the full threefold; the only
constraint is the rigid 3-adic boundary condition.**  This definitively
explains the contact-chart failures as chart artifacts.

## Main test 3: CRT-guided rational search, height <= 16

`./c35_sweep search H` sieves all `(q0,q1,q2)` with `|num|, den <= H` through
per-prime lookup tables at `p in {3,5,7,11,13,17,19,23,29}` (kill only when
the reduction is a smooth genus-2 chart point failing the prime-to-`p` part of
35; boundary/nonintegral = conservative pass).  Post-processing
(`code/claude_c35_post.gp`) drops globally degenerate triples (planes
`q0 = 0`, `q2 = 0` — they ride the conservative pass) and extends the exact
sieve `35 | #J(F_q)` through `q <= 97`.

```text
H=8 :  triples 658503,   genuine survivors of C-stage 32,   final 0
H=12:  triples 6128487,  genuine survivors 705,             final 0
H=16:  triples 32461759, genuine survivors 6551,            final 0
H=16 kill distribution beyond 29:
   31:5317  37:1025  41:174  43:30  47:3  59:2  (then empty)
```

The kill curve decays geometrically at the expected `~0.15-0.18` per prime
with no anomalous tail: **no rational `C35` point on the `A_1(5)` chart with
coordinate height <= 16**, and no sign of a low-height positive-dimensional
35-locus (the split HLP component evidently does not meet this chart at small
height in these coordinates).  Wall time: the whole pipeline (validation +
landscape + densities + three searches) is ~7 CPU-minutes.

Artifacts (scratchpad `/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/`):
`c35_p3_full.txt`, `c35_search_h{8,12,16}.txt` (+ `.log`).

## Verdict

Upgraded from "locally obstructed in all known charts" to **locally
unobstructed on the full 5-torsion moduli threefold**.  The mod-3 wall was an
artifact of the `n(P - infinity)` contact charts; on `A_1(5)` the necessary
conditions have positive density at every prime except the (non-obstructing,
rigid) 3-adic boundary condition.  No candidate through chart height 16; since
`A_1(5,7)` is a high-degree cover (presumably general type), expect finitely
many non-split rational points, plausibly at moderate height — a production
search is justified and cheap.

## Next steps

1. **Production search, height 40–60**: move the `31..43` tables into the C
   stage, quotient by the involution `(q0,q1,q2) -> (q2,1-q1,q0)` (halves the
   space), and parallelize by `q2`-slice.  Cost scales like `H^6`; `H=40` is
   ~ 2*10^9 triples of table lookups — hours, not days.
2. **Sharpen the 3-adic filter** (biggest lever): currently every 3-boundary
   triple passes.  Blow up the chart mod 3 to determine which nonintegral
   valuation patterns of `(q0,q1,q2)` can reduce to the 12 good `F_3` curves
   (all with charpoly `(x^2+x+3)(x^2+3x+3)`) vs genuinely bad reduction; kill
   the patterns that force good reduction with the wrong Weil polynomial.
   The analogous 5-adic and 7-adic boundary analysis is also unexploited.
3. **Bad-prime strengthening**: add component-group/generalized-Jacobian
   necessary conditions at bad sieve primes (the Howe control passed via bad
   11, 13 — real candidates will do the same; a Néron-model-order filter at
   bad primes would cut the conservative-pass tail).
4. **Route B in parallel**: reconstruct the NDE 2003 `A_1(7)` threefold and
   run the mirror sieve (`5 | #J`); compare which chart exposes `A_1(5,7)`'s
   rational points at lower height.
5. For any final survivor: exact `TorsionSubgroup` in Magma, then geometric
   simplicity via the 12th-power-transform certificate — at more than one
   prime (see the biquadratic `L_5` caveat above), plus a split check (a hit
   on the HLP component would show split Frobenius factorizations at most
   primes, like Howe's curve at 17, 19).
