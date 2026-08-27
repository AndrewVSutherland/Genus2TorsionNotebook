# prod 01: (3,12) — FOURTH [3,12] hit found (split); two simple Q(zeta_3) near-misses; both fiber quotients rank 2; no low-genus fibers

Date: 2026-07-18.  Production follow-up to notes/claude_top10_01_312.md (which
established: S12 carrier surface has a genus-4 fibration pi_z over the z-line;
z=-5/3 fiber holds 6 known rational points = two split-hit triples).  All new
scripts/logs in scratchpad `t312prod/` (see Files below).  ~55 CPU-min total.

## Strategy recap (3 lines)

(A) Reconstruct the z=-5/3 fiber of S12 as an explicit curve over Q and verify.
(B) Exploit its (newly found) degree-3 etale map to a genus-2 curve: Mordell-Weil
+ Chabauty on the quotient to classify all rational points of the fiber.
(C) Interpolate F12(U;z,r) mod p in both variables; scan ALL z in F_p for
degenerate/low-genus fibers of pi_z; boundary checks.

## (A) The fiber over Q — DONE, verified

- `data/claude_prod_01_312_fiber_z53.m`: G(W;r), monic degree 12, W = (r+1)U,
  coefficient numerators degree <= 10, denominators (r+1)^(12-j).  Built from
  ~58 per-point Q-level eliminants (qbatch_a.log) + interpolation, 6 held-out
  points verified.
- Verified: full match mod 10007 against the direct function-field slice
  (modpcmp.m: true); divisibility of the degree-40 minpoly at 6 points mod
  31013 (selfcmp31013.m: all true); all 6 known hit points satisfy G = 0
  exactly (verify_q.log).
- Genus over Q: **4** (FunctionField/Genus, 376 s, verify2.log) — matches mod-p.

## (B) Genus-2 quotient, Mordell-Weil, new point — MAJOR RESULTS

1. The degree-12 function field has a UNIQUE nontrivial subfield over Q(r):
   degree 4, genus 2 (verify2.log; predicted mod p in subfields_p.log).  So the
   genus-4 fiber C4 has a **degree-3 map to a genus-2 curve E**, etale by
   Riemann-Hurwitz (6 = 3*2).  Defining quartic g4(W;r) in
   `data/claude_prod_01_312_quotient.m`.
2. Hyperelliptic model over Q (IsHyperelliptic + minimization, q2curve.log):
   **Hm: y^2 + (x^2+1) y = x^5 + 61x^4 - 86x^3 - 624x^2 - 5657x + 18400**,
   simplified Hs: y^2 = 4x^5 + 245x^4 - 344x^3 - 2494x^2 - 22628x + 73601.
3. Each known split triple maps to ONE point: r=-3/5 -> (5,-144),
   r=-9/35 -> (5,144) (hyperelliptic-involution partners).
4. Rational points found on Hs (bound 1e5, re-checked to 2e6 — see hsearch.log):
   infinity, (5,±144), (11/4,0), **(-16,±3615)** (new pair).
5. **MW(Jac(Hs)) = Z/2 x Z x Z proved** (MordellWeilGroupGenus2, flags
   true/true, 143 s, q2mw.log).  Rank exactly 2 = genus => classical Chabauty
   NOT applicable.  Height-pairing det of the two visible generators 9.456.
6. (r,x)-correspondence Phi of bidegree (2,4) (phi.m): over x=-16 the places
   have r = 10161/6025 and r = infinity; over x=11/4: r = -81/800 (its fiber
   cubic has NO rational root); over x=inf: r = inf.
7. **NEW RATIONAL POINT of the z=-5/3 fiber** (first non-split rational point
   of S12): r = 10161/6025, W = -1549/3615 (U = -7745/48558).  Contact data
   M = -16767142144/16875 < 0 (NOT a square), V = 227028025/37726069824.
8. Jackpot protocol on the corresponding curve at (z,r)=(-5/3, 10161/6025)
   (newpt.m): torsion = **Z/12** (not (3,12) — the extra 3-class is only
   Q(sqrt(M))-rational); **geometrically simple** (L_p irreducible + 12th-power
   transform irreducible at p=29 AND p=53); disc primes {2,3,5,19,241,379,739,
   1129,8093}; G2-invariants differ from both split hits.  Integral model in
   `data/claude_prod_01_312_quotient.m`.
   **M = -3 (129488/225)^2**: the extra 3-class is rational exactly over
   Q(zeta_3) — a geometrically simple Jacobian /Q with Z/12 torsion growing to
   contain (3,12) over the cyclotomic field (closest simple near-miss so far).
9. **Structural correction to the campaign picture**: rational points of S12
   do NOT automatically give [3,12] — the extra 3-torsion class is rational iff
   additionally M is a square.  The true carrier is the etale double cover
   m^2 = M of S12.  M is square at all 6 split points, nonsquare at the new
   point; M is NOT constant along the degree-3 fibers (256, 256/81, 256/289 at
   the r=-3/5 triple), so the cover does not descend to E.
10. CONDITIONAL CLASSIFICATION: if Hs(Q) = {the 6 points in 4}, then the
    z=-5/3 fiber has exactly 7 rational points: 6 split (3,12) points + the new
    non-split Z/12 point, and hence NO new (3,12) curve on this fiber.  Only a
    Mordell-Weil sieve (rank 2, MW known exactly, both generators explicit) is
    missing to make this unconditional.

## (B'') Uniformity: a genus-2-fibered quotient surface of S12

The degree-4 genus-2 subfield exists in EVERY tested fiber: z = -5/3 (over Q)
and z = 17, 50, 101, 137 (mod 10007, subf_family.m), each time unique.  So S12
admits a global degree-3 fiberwise-etale quotient map to a surface T with a
**genus-2 fibration over P^1_z** (spreading out from the uniform fiberwise
structure).  Every rational point of S12 maps into T(Q); the generic fiber of
T is a genus-2 curve over Q(z), opening a family-Mordell-Weil / fiberwise-
Chabauty attack as the systematic continuation (the z=-5/3 and z=-5/4 fibers
of T are exactly the quotients computed in (B) and (B')).

## (B') Second fiber z=-5/4 — COMPLETED, with a FOURTH [3,12] HIT

`t54/`: Q-level qbatch (54 points), interpolation verified (6 held-out), fiber
at r=-32/35 has the expected 3 rational W-roots.  Genus over Q: 4 (644 s).
Unique degree-4 genus-2 subfield again; quotient minimal model
**y^2 = x^5 + 80x^4 - 11x^3 - 90x^2 - 675x - 675**; points found (bound 1e5):
infinity, (-1,0), (3,+-54), (-79/4,+-194430/64);
**MW = Z/2 x Z x Z, rank exactly 2 proved AGAIN** (pipe54.log) — Chabauty
blocked at both computed fibers by rank = genus = 2.
Correspondence Phi (bidegree (2,4), phi54.m) lift results:
- x=3 covers r=-32/35 (the known 3rd split hit, 3 roots, M squares) AND
  **r=-32/65 — a previously unknown point with 3 rational W-roots, ALL M
  square** => FOURTH EXACT [3,12] REALIZATION (newpt54.m):
  curve y^2 = 52x^6 + 156x^5 - 1043x^4 - 2346x^3 - 629x^2 + 570x + 225,
  TorsionSubgroup = [3,12] verified; G2-invariants new; but SPLIT
  (aut order 4; elliptic quotients both conductor 4290, torsions Z/12, Z/6;
  split54.m).  Data: `data/claude_prod_01_312_newhits.m`.  NOTE: height 65
  parameter — beyond all previous residue-scan reach; the quotient method
  finds hits the direct scans missed.
- x=-79/4 covers r=811171/648100: SECOND SIMPLE NEAR-MISS —
  M = -3 (1459271/50)^2 (Q(zeta_3) again!), torsion Z/12, geometrically
  simple (certificate at p=29,67,71,83,97,101,103).
- x=-1 (Weierstrass) covers r=-128/225 (double): fiber [3,9], no rational
  root; r=inf places only over x=-79/4 and x=inf; no rational x over r=-1.
**Pattern/conjecture**: every non-split rational point of S12 seen so far has
M in -3 (Q^*)^2, i.e. 3-class rational exactly over Q(zeta_3) — suggests a
mu_3/Weil-pairing constraint forcing square class of M into {1, -3} on S12;
if true, the (3,12) double cover of S12 is governed by a single -3 twist.

## (C) Degenerate-fiber scan — COMPLETE at p=10007

Method: 19880-point mod-p grid (140 z-columns x 142 r-samples, 4 ms/point,
sampler40.m), per-column rational-function reconstruction of the degree-40
minpoly + factorization over F_p(t) (recon_cols.m; nullspace-based Cauchy
interpolation, 8 held-out points per function), then 2-variable reconstruction
of F12(U;z,t) (twovar.m; max z-degree 24) and a disc_U-pattern scan over ALL
z0 in F_p (27 s).

Results (twovar.log, recon_a.log, colssp.txt):
- All 140 sampled columns: factor pattern [1,12,27], fiber genus 4, disc
  pattern <150,67,132>.
- Full scan of 10007 z-values: generic disc pattern <294,226> at 9975 values;
  deviations only at 29 finite z (gcd degree +1..+3) and z = 0, 1, -1
  (identically degenerate: the 40-dim contact algebra collapses there —
  gridsp.txt shows all-142 DIM/EMPTY flags at z=0,1,-1).
- Genus recheck at ALL strongest deviations (z=2641,7366 <+3>, 4693,5766 <+2>,
  and three <+1> samples): **fiber genus is STILL 4** in every case
  (branch-point collisions, not genus drops).
- CONCLUSION: **pi_z has NO fiber of genus <= 3 over any z in F_10007 outside
  the boundary {0, +-1}**.  In particular no genus <= 1 fiber => no sweepable
  subfamily in this pencil.
- Second-prime cross-run at p=31013 (grid_b/cols_b/twovar_b): same generic
  pattern <294,226> (30998 of 31013), same poles z=0,1, same identical
  degeneration at z=-1, deviations at 12 finite z-values.  Rational
  reconstruction of deviation values at BOTH primes: **empty intersection** —
  every height <= 70 "recognition" was a coincidence.  So the degeneracy
  divisor of pi_z has NO rational point of height <= 70, and all its fibers
  (any prime tested) keep genus 4 anyway.
- Boundary of the z=-5/3 fiber itself (phib.m, unconditional): Phi(-1,x) has
  no rational root => no rational place over r=-1; all r=inf places lie over
  x=-16 and x=inf.

## Verdict

- The quotient method WORKS as a discovery engine: one session produced a
  FOURTH exact [3,12] realization (split, new curve, parameter height beyond
  all previous scans) and TWO new geometrically simple curves whose extra
  3-class is rational exactly over Q(zeta_3).
- Simple-(3,12)-over-Q remains open, but the obstruction is now sharply
  localized: on both computed fibers, every rational point is either split
  with M square, or simple with M in -3 (Q^*)^2.  If the {1,-3} square-class
  dichotomy on S12 is a theorem, the hunt becomes precisely "kill the -3
  twist" — the clearest structural target the campaign has had.
- Both fiber quotients have MW rank exactly 2 = genus (proved), so per-fiber
  classification needs a Mordell-Weil sieve, not Chabauty.
- No low-genus fiber shortcut exists in pi_z (settled at 2 primes).  Route 2
  (BFS sqrt3-division, degree-9 carrier over rational P^2) remains the best
  structural hedge.

## Next steps

1. PROVE/refute the M square-class dichotomy {1, -3} on S12 (check it first
   at many mod-p points of S12; then attempt the mu_3/Weil-pairing argument).
   If true, study the single quadratic twist obstruction globally.
2. Mordell-Weil sieve on both quotients (rank-2 MW groups fully known,
   generators explicit): make the fiber classifications unconditional.
3. Iterate the quotient pipeline over other small-height z (each fiber ~20
   min: qbatch -> interp -> subfield -> quotient -> MW).  Any quotient with
   rank <= 1 gives an unconditional Chabauty fiber; any new x-pair gives new
   S12 points to protocol.  The genus-2-fibered quotient surface T (B'')
   also invites a family/MW-lattice attack over Q(z).
4. Study the M-square double cover C8 -> C4 / S12~ -> S12 (the true (3,12)
   carrier); with the -3 dichotomy this is a single twist family.
5. The two simple Z/12 curves with Q(zeta_3) 3-class are independently
   bankable (torsion growth over cyclotomic fields); in
   `data/claude_prod_01_312_quotient.m` and `data/claude_prod_01_312_newhits.m`.

## Files (scratchpad `t312prod/`)

- Task A: ppelim.m/ppelim2.m (per-point eliminants), interp.m, f12_q_raw.m,
  modpcmp.m, selfcmp31013.m, verify_q.log, verify2.m/.log (genus+subfields /Q).
- Task B: subf2model_p.m (mod-p prototype), g4_q.m, q2curve.m/.log,
  q2rank.m/.log, q2mw.m/.log, q2lift.m, q2lift2.m (fiber lifts), phi.m
  (correspondence), newpt.m (jackpot protocol), hsearch.m/.log.
- Task B': t54/ (mksamples.py, pipe54.m/.log, qbatch54.log, g54_q_raw.m,
  g54sub.m, phi54.m (lifts; found the 4th hit), newpt54.m (jackpot protocol,
  torsion [3,12] verified), split54.m (split structure + 2nd near-miss)).
- Task C: sampler40.m, recon_cols.m, degtest.m, twovar.m/.log, grid_a.txt,
  cols_a.txt, gridsp.txt, colssp.txt; second prime: twovar2.m, grid_b.txt,
  cols_b.txt, recon_b.log, twovar_b.log.
- Data: data/claude_prod_01_312_fiber_z53.m (fiber model G),
  data/claude_prod_01_312_quotient.m (quotient, MW, new point, new curve).
