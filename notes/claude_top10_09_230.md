# top10 #9: (2,30) via Elkies' A_1(5) threefold — strategy + test run (GO)

Target: `[2,30]` (order 60) on a geometrically simple genus-2 Jacobian over Q.
Status quo: products only — not even a split genus-2 Jacobian is known. C30
simple IS known (Nicholls 2018; plus this repo's infinite contact-5/contact-6
`[30]` family). The in-family extra-2 route is CLOSED unconditionally
(`notes/claude_tier2_230_order30_extra2.md`: quadratic branch empty by
Chabauty0 theorems, cubic branch = genus-12 curve clean to height 5000), so a
first `[2,30]` needs a different order-30 construction. The team's designated
fresh route (June 26 todo item 3, never previously set up): combine the
A(6)/M(6) 6-torsion construction with Elkies' LuCaNT-2024 5-torsion moduli
threefold — the engine that produced the (2,2,2,10) record.

## Headline results of this test run

1. **Calibration passed.** Repo fact reproduced: the (u,s)=(125,5415)
   contact-5/6 member has Order(D5)=5, Order(D6)=6, Order(D5+D6)=30, exact
   torsion `[30]`, simplicity cert at p=13. Literature fact reproduced:
   Elkies' curve `y^2 = x(x+1)(x-1)(3x-7)(8x-13)(24x+25)` has torsion
   `[2,2,2,10]`. Elkies' universal 5-torsion class verified: for 5 random
   `(q0,q1,q2)` the class `[{Q=0,y=0} - K]` has exact order 5, and the
   involution `(q0,q1,q2) -> (q2,1-q1,q0)` preserves G2-invariants.
2. **No local obstruction (the key go/no-go).** Full F_p^3 scans of the Elkies
   chart at p = 7..23: the (2,30) necessary condition {2-rank(J(F_p)) >= 2 and
   3 | #J(F_p)} holds on a positive fraction of good fibers at EVERY prime
   (0.9% at p=7 rising to 4.6% at p=23). Unlike the repo's contact charts
   (dead at small primes), the full A_1(5) space has good-reduction room.
   Self-test: 5 | #J(F_p) on 100.0% of good fibers at all six primes
   (23,614/23,614) — validates both the built-in 5-torsion and the
   point-count engine (also cross-checked against gp `hyperellcharpoly` on 10
   random fibers, exact #J agreement).
3. **The (2,10) sub-stratum is populated by SIMPLE members** (only the 3-part
   is missing for (2,30)): a 14-prime sieve over rational `(q0,q1,q2)` of small
   height found **five distinct geometrically simple curves with exact torsion
   `[2,10]`**, all on the odd-degree boundary `q2 = -1/4`
   (`data/claude_230_simple_210_curves.txt`), e.g. `q = (-4, 3, -1/4)`:
   `y^2 = 6x^5 - 167x^4 + 1520x^3 - 5152x^2 + 7424x - 3840`, cert p=7.
   (2,10) is already a realized group — the point is pipeline validation: the
   Elkies chart + extra-2 machinery reaches simple members at tiny height.
4. **The 3-part will not come from sieving alone:** 0 survivors of the
   {3 | #J everywhere} level at height <= 15 on the 3-fold (1.2M triples) and
   0 on the `q2=-1/4` surface at height <= 60 (386K pairs) — consistent with
   the repo lesson that torsion must be FORCED by identities. Production step
   = algebraic 3-torsion imposition (Route 1 below).
5. Structural facts learned (useful for production): (i) survivors on the
   involution-fixed line `(t, 1/2, t)` never certify simple (extra-automorphism
   locus, as expected) — includes a curious exact-`[5,10]` member at t=3/5
   (flag for the (5,10) target; likely split, no cert p <= 97); (ii) the line
   `q1 = -(3/4) q0` on the `q2=-1/4` surface is a systematic sieve
   false-positive: quartic cofactor with biquadratic-type Galois (local 2-rank
   2 at every p, global 2-rank 1, torsion `[10]`).

## Setup

Elkies (LuCaNT 2024, Contemp. Math. 796): universal curve for a genus-2
Jacobian with a marked rational 5-torsion class:

```text
y^2 + (L'Q' - LQ) y = Q^2 Q',   L = x, L' = 1,
Q  = q2 x^2 + q1 x + q0,        Q' = Q - L L' = Q - x.
```

Coordinates (q0,q1,q2) birational on the moduli threefold A_1(5); T -> 2T is
(q0,q1,q2) -> (q2, 1-q1, q0); the 5-torsion class is [{Q=0, y=0} - K].
Completing the square:

```text
y^2 = F(x),   F = (Q' - xQ)^2 + 4 Q^2 Q',   deg F = 6,
lc(F) = q2^2 (1 + 4 q2)   [so q2 = -1/4 gives the deg-5 boundary chart].
```

Since (2,30) = Z/2 x Z/2 x Z/3 x Z/5, on top of the built-in 5 we need:
(a) 2-rank(J(Q)) >= 2, i.e. F has >= the right rational factorization
(deg-6: three stable even orbits; deg-5: quintic with >= 3 factors), and
(b) a rational 3-torsion class. Necessary at every odd good prime p:
2-rank(J(F_p)) >= 2 and (p != 3) 3 | #J(F_p); 5 | #J(F_p) is automatic.

## Strategy (ranked routes)

**Route 1 (main, the designated 6x5 combination).** Intersect A_1(5) with the
3-torsion (equivalently 6-torsion) locus, staying in the chart validated
today. Concretely, impose a cubic-contact 3-torsion class on `y^2 = F(x)`:
exists cubic h3, quadratic k, scalar lambda with

```text
h3(x)^2 - F(x) = lambda * k(x)^3,
```

(the repo's standard plus-3 contact condition), a system in
(q0,q1,q2; h3; k; lambda) to be reduced by elimination. Counting: 3 (chart) +
4 (h3) + 3 (k) + 1 (lambda) - 6 (coefficient equations of the sextic identity)
- 1 (scaling) - 1 (x-translation redundancy... none, chart is rigid) leaves an
expected ~2-dim locus {(C, T5, T3)} — the '(15)-locus' — matching the moduli
count (finite cover of M2 minus nothing: pairs (C,T15) form a 3-fold; the
contact form h3 restricts to classes [P1+P2-K], a coindex-0 open part).
Then cut by the extra-2 condition (rational quadratic factor of F, validated
today on the `q2=-1/4` boundary where one Weierstrass point is free) to reach
(2,30) = 2 x 2 x 3 x 5. Recommended implementation: work on the `q2 = -1/4`
surface first (all five simple (2,10) hits live there): 2 chart params,
so {3-torsion contact} should give a CURVE in (q0,q1); compute its genus; if
<= 2, run the rank/Chabauty pipeline; independently scan its F_p fibers
for the extra-2 condition (same engine as today).

**Route 2 (mirror chart: A(6) + forced extra 2, then 5).** A(6) model
(NotesAndTodo.tex): `y^2 = R (R x^2 + 4(R+x-1)(R-1))`, R = a x^2 + b x + c;
[R,0]-type class of order 6 built in (3-part AND one 2-torsion class free).
The cofactor is quadratic in R: `G = 4R^2 + R(x^2+4x-8) - 4(x-1)`. Extra
2-torsion class <=> G = G2 * G2' with rational quadratics <=> the quartic

```text
disc_R(x) = (x^2+4x-8)^2 + 16*G0(x) + 64(x-1)
```

is a perfect square of a quadratic when G0 = G2*G2' is prescribed — i.e.
choose monic G2, G2' (4 free coefficients), impose 2 square conditions: a
2-dim (2,6)-forced surface with explicit equations. Then impose 5 by sieve
(5 | #J at all good p) or, better, by GL2-equivalence of `R*G` with
`F(q0,q1,q2)` (full 6x5 moduli intersection — same 3-fold as Route 1 seen
from the other chart). Caution: the M(12)+5 analog died at p=7
(`notes/order60_attempts.md`); run the F_7 room scan on this surface FIRST
(one scan230.c variant, minutes).

**Route 3 (overshoot, cheapest big win).** Elkies' Clebsch–Klein cubic
surface family already forces (2,2,2,10); a member with rational 3-torsion
would have torsion >= (2,2,30) [invariant factors], containing (2,30) and
realizing an even larger new group. 2-parameter rational family, so a
3-divisibility residue sieve + exact checks is a direct reuse of today's
engine. Lesson-based expectation is low (sieved 3-part), but the cost is one
afternoon and any hit is spectacular.

## Test run: exact commands and results

Workdir: `/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/work230`
(code copied to `code/claude_230_*`).

```text
magma -b validate230.m            # calibration (a),(b),(c) above     [~40 s]
gcc -O2 -march=native -o scan230 scan230.c
./scan230 scan                    # full F_p^3 scans, p = 7..23       [~2 s]
./scan230 base                    # random-sextic baseline, 20K/prime
./scan230 one p q0 q1 q2          # single-fiber #J/rank (gp cross-check)
gp -q crosscheck.gp               # independent #J via hyperellcharpoly
./scan230 sieve 15 5              # 3-fold sieve: |n|<=15, den<=5, 14 primes (3..47)
./scan230 sieve5 60 8             # q2=-1/4 surface sieve: |n|<=60, den<=8
magma -b checksurv.m / checksurv5.m   # exact torsion + simplicity certs
```

Scan table (Elkies chart, all of F_p^3; good = smooth deg-5/6 fiber):

```text
p   total  good   5|#J   3|#J   rank>=2  BOTH  frac   baseline-BOTH
7    343    215    215     69      17      2   0.93%   3.7%
11  1331   1033   1033    382     115     30   2.90%   5.4%
13  2197   1727   1727    584     171     57   3.30%   4.6%
17  4913   4095   4095   1658     481    174   4.25%   5.7%
19  6859   5897   5897   2097     739    223   3.78%   5.2%
23 12167  10647  10647   4438    1329    494   4.64%   6.2%
```

(baseline = 20000 random sextics per prime; baseline 5|#J ~ 20-24% vs 100%
on-locus.) Mild suppression vs baseline at small p, but no obstruction.

Sieve results: 3-fold height-15: 21 level-A survivors (rank>=2 at all good
p in {3,...,47}), 0 level-B (+ 3|#J). Exact checks: 6 simple `[2,10]` hits =
3 distinct curves (involution pairs matched exactly, e.g. (-1/4,-2,-9) =
image of (-9,3,-1/4)). Surface sieve (q2=-1/4, height 60): 28 level-A,
0 level-B; exact checks give 5 distinct simple `[2,10]` curves total
(certs p = 7,7,11,11,17), 12 further `[10]`-only members on the
`q1 = -(3/4) q0` false-positive line, 4 exact-`[5]`, plus the non-certifying
involution-fixed-line members incl. the `[5,10]` curiosity at (3/5,1/2,3/5).

CPU: ~15 min total, single-threaded, nice -n 15 throughout.

## Verdict

**GO.** The designated 6x5 route survives its first contact with reality on
every axis tested: the chart is computationally friendly, the necessary local
conditions have room at all small primes (the killer of previous (2,30)
attempts), the (2,10)-with-simplicity stratum is demonstrably populated at
tiny height on an identifiable boundary surface (q2 = -1/4), and the
remaining gap is exactly one algebraically-imposable condition (the 3-part),
for which the repo already owns the contact machinery. The realistic hard
step is the elimination for the {T5, T3} locus and the geometry (genus) of
its `q2=-1/4` slice.

## Next steps

1. Derive the plus-3 cubic-contact system `h3^2 - F = lambda k^3` on the
   `q2 = -1/4` surface (2 params): eliminate (Groebner/resultants, PARI or
   Magma) to a plane curve in (q0,q1); compute genus + rational points. If a
   rational curve: instant infinite `[15]/[30]`-type family inside A_1(5),
   then extra-2 cut for (2,30). If genus 1-2: rank/Chabauty pipeline.
2. Same elimination on the full 3-fold (expected surface): even a mod-p
   fiber count of the {T5,T3} locus at p = 7,11,13 (extend scan230.c: for
   each F_p-triple, test existence of the contact (h3,k,lambda) over F_p by
   brute force over k) would locate the right chart before any symbolic work.
3. Route 2 opener: F_7/F_11 room scan of the explicit (2,6)-forced surface
   {G2, G2' : disc_R square} for 5 | #J — one scan230.c variant. If p=7 is
   dead as in the M(12)+5 attempt, kill Route 2 early.
4. Deeper level-B sieve with the existing engine (3-fold to height ~40,
   surface to height ~500, early-abort makes this minutes) — cheap lottery
   tickets while the elimination is developed.
5. Housekeeping flags: (i) hand the `[5,10]` member (3/5,1/2,3/5) to the
   (5,10) effort for a proper split/simple decision; (ii) the five simple
   `[2,10]` curves are new explicit small examples on Elkies' threefold —
   possibly worth an LMFDB conductor check.

## Files

```text
code/claude_230_validate_elkies5.m   calibration: repo [30] member, Elkies (2,2,2,10), order-5 class
code/claude_230_elkies5_scan.c       F_p^3 scans, baseline, single-fiber, sieve + sieve5 modes
code/claude_230_checksurv.m          exact torsion + simplicity certs for survivors
data/claude_230_simple_210_curves.txt  the five simple [2,10] curves (models + certs)
(scratchpad work230/: validate230.out, sieve_h15.out, sieve5.out, crosscheck.gp)
```

## Caveats

* Level-A/B sieve conditions are NECESSARY only; 3 | #J at all tested primes
  does not imply rational 3-torsion (exact Magma check remains the arbiter).
* cert_p = 0 (no certificate p <= 97) is evidence of, not proof of,
  non-simplicity; the involution-fixed-line members were not further analyzed.
* Bad-reduction fibers are always "allowed" by the sieve; globally singular
  triples are only weeded out at the exact-check stage.
* The moduli-dimension counts for the {T5,T3} locus are heuristic until the
  elimination in Next step 1/2 is actually run.
