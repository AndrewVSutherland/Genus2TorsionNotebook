# Top-10 campaign summary — 2026-07-17

Companion to `claude_top10_ranking.md`. One strategy+test session per target (files
`claude_top10_01..10_*.md`); every test ran single-threaded/niced, ~30 CPU-min budget each.
Headline outcomes, ranked as in the ranking file:

| # | Target | Test executed | Outcome / verdict |
|---|--------|---------------|-------------------|
| 1 | (3,12) | S12 slice-genus computation + chart-free F_2/F_5 enumerations | **S12 is a genus-4-fibered surface** (z-fibers genus 4, not general type); hit line z=-5/3 is a nondegenerate genus-4 fiber /Q with >=6 known (split) points — per-fiber Chabauty now concrete. NEW THEOREM: any genus-2 Jacobian /Q with (3,3) ⊆ torsion has bad reduction at 2, and with (3,12) also at 5 (exhaustive F_2/F_5 enumerations). Outlook IMPROVED. |
| 2 | (2,24) | **Mechanism B on G1 — first-ever run** | Reduces to a genus-13 cover with genus-3 non-hyperelliptic quotient X_B (plane quartic computed); 7 rational points found, all degenerate; rank(Jac X_B)>=1, tors|3, simple — **no rank-0 shortcut; G1 almost certainly sterile but needs plane-quartic descent + Chabauty–Coleman to close**. Realization hope now rests on the unfitted THIRD component (>=10 unexplained points on the z=25/7 fiber). |
| 3 | (4,16) | **Per-R elliptic solve built and validated** — the O(H^4) wall is down | For fixed R the membership condition is two explicit genus-1 quartics; 0.1 s/fiber to w-height 1e5. Re-found all known points, and found a **NEW certified-simple [4,8] point at height 5093** (6.4x beyond the old frontier, ~4000x cheaper than blind sweep). ellrank closes rank-0 fibers unconditionally (R=-8,2,3/2,5 closed). GO for production (R-height 100, w-height 1e6). |
| 4 | (2,2,2,12) | H=800 sweep reproduced (64 near-misses exact) + rank computation on the T3 pencil | **Spec'd rank-0 kill REFUTED**: T3 near-miss locus fibers over the conic pencil λB3+μB5; on each member the 3 passing conditions collapse to ONE isomorphism class of elliptic curves with rank exactly 2 resp. 3 — the algebraic-correlation mechanism that hits require exists in 3 of 4 conditions. Whether a pencil member exists where the 4th condition joins is now a FINITE symbolic question (rho-scan of resultants). Major structural upgrade. |
| 5 | (2,6,6) | Chart scan extended H=60→150 (1.5e9 pairs) | 4 distinct [2,6,6] classes now known (3 new) — **all bielliptic/split; certified: no simple [2,6,6] on this chart to height 150**. Thin-set confinement to the split locus (same phenomenon as [3,12]/S12). Escape route: BFS sqrt3-RM (3,3)-first construction with conic 2-rank conditions. Also: the (6,6) curve's strict certificate prime is p=37, not 23. |
| 6 | (2,2,4,8) | **First-ever enumeration of the twisted second component** (d<=3000) + genus-3 split-locus search | 1947 twisted tuples, delta(T_0) trivial 1947/1947 (consistency perfect), 0 (2,2,4,8) hits; best-square stats mirror component I (4/5 parity-blocked, 3/5 at ~11%) — log-divergent expected count justifies 2-3 more orders of magnitude. Genus-3 split curve empty to height 2000; E1,E2 rank bounds (1,1) → elliptic-Chabauty over Q(sqrt 2) is the decidable next step. |
| 7 | (2,2,16) | Fresh-chart divisibility test (25 twisted curves × all 32 order-8 classes) + d=-2 fiber descent | Chart 2: 0 halvings, uniform 1/5 delta ceiling (~800 curves) → symbolic parity audit required before scaling. Chart 1: the three d=-2 fibers have **unconditional rank 2 → Y(Q) is infinite**, yet 0/14 new points pass the field condition — Brauer–Manin hypothesis strengthened and now certifiable fiber-by-fiber. Exact (2,2,16) can only come from the 7 never-enumerated twisted norm surfaces Y_t. |
| 8 | (8,8) | Nicholls Prop 5.9.6 (4,4)-family implemented + validated + box scans | **Substrate GO**: 10/10 members have rational (4,4), generically simple — a simple-(4,4) production line on a rational 3-fold (strictly richer than the dead M_1(8,4) chart). No (8,4)+ at height <=4; halving loci are thin. Decisive next step: geometry (genus profile) of the stage-1 halving cover before any search. |
| 9 | (2,30) | 6x5-combination local go/no-go on Elkies' A_1(5) threefold + sieve | **GO — the contact-chart local death does NOT recur**: (2,30)-room at every prime tested (0.9%→4.6% densities). Bonus: **five new certified-simple [2,10] curves** on the q2=-1/4 boundary (data/claude_230_simple_210_curves.txt) + an exact [5,10] member flagged for the (5,10) target. Next: plus-3 cubic contact on the q2=-1/4 surface. |
| 10 | (35) | Full A_1(5)+7|#J landscape at p=3..71 + CRT search to height 16 | **Upgraded to GO**: F_3 supports #J=35 (12 models, unique Weil poly (x^2+x+3)(x^2+3x+3)) — the contact-chart mod-3 obstruction is a chart artifact. Positive densities at every prime; search clean to height 16 (all 6551 candidates killed by q<=59). Production height 40-60 sieve justified (hours). Caveat found: biquadratic L-polys (e.g. p=5 on the [28] curve) do NOT certify simplicity — scan primes. |

## Corrections to project records found during testing
- The simplicity certificate must SCAN primes: p=13 fails for (29,121,125,145) (works at 61);
  p=11 fails for Nicholls' [30] curve (works at 13); p=23 gives only the weak certificate for
  the (6,6) curve (first strict prime 37); p=5 is biquadratic for the [28] curve.
- The two "distinct" [2,6,6] examples in contact6_m36.md were already known Q-isomorphic;
  now 4 distinct classes exist to H=150, all split.

## New banked by-products (verify + add to tables)
- Simple [4,8] at (R,w)=(1/11,-1843/5093) on M_1(8,4), cert p=43 (top10_03 notes).
- Five simple [2,10] curves on Elkies' A_1(5), q2=-1/4 (data/claude_230_simple_210_curves.txt).
- Exact [5,10] member at q=(3/5,1/2,3/5) — split/simple status UNDECIDED (no cert p<=97);
  potential first-ever (5,10) if simple: top-priority follow-up.
- Three new split [2,6,6] classes (top10_05 notes) for the split census.
- Universal theorem: (3,3) ⊆ J(Q)_tors forces bad reduction at 2; (3,12) forces bad at 2 and 5.
