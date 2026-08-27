# Production campaign summary — 2026-07-18

Ten parallel production lanes (notes/claude_prod_01..10_*.md), executed on the revised top-10
(notes/claude_review_gpt56_plan.md). Wall ~1.9h, ~16 threads. Hits first.

## HEADLINE: (2,2,2,12) REALIZED — new group, new record (order 96)

Curve (M(2,2,2,6) chart, (s,m,n) = (336396, -689185, -166464)):
  y^2 = prod_{i=1..5} (A_i + B_i x),  A = [1,1,1,2,2],
  B = [282322361376, -8243383980, -64241207724, -114724491840, 561915878400]
Reduced minimal model:
  y^2 + (x^2+1)y = 756900x^6 + 737595570x^5 + 150572203590x^4 - 15854483576121x^3
                   - 530648977741620x^2 + 32014154874551031x + 830742747091037849
- TorsionSubgroup EXACTLY [2,2,2,12] (order 96). INDEPENDENTLY RE-VERIFIED from scratch
  (orchestrator Magma run): torsion + certificates at p = 37, 73, 113, 149 (Frobenius charpoly
  irreducible AND 12th-power transform irreducible of degree 4 at each).
- Geometric simplicity: 4-prime certificate excludes geometric splitting; QM is excluded
  unconditionally by Laga-Schembri-Shnidman-Voight (O-PQM Jacobians have #tors <= 16 < 96).
- Significance: (i) first (2,2,2,12) on ANY genus-2 Jacobian /Q (previously known only for
  products of elliptic curves, not even split Jacobians); (ii) largest torsion group ever
  recorded on a geometrically simple genus-2 Jacobian /Q — beats Elkies' (2,2,2,10), order 80;
  (iii) found structurally (T5-pencil point u=-97/48 on member rho'=-49/240), not blind search;
  second isomorphic representation at (u,rho')=(133/145, 289/240).
- Data: data/claude_prod_02_22212_hit.txt; full derivation notes/claude_prod_02_22212.md.

## Other lane outcomes (one line each; details in claude_prod_NN files)
- (3,12): 4th exact realization found at (z,r)=(-5/4,-32/65) — split, beyond all prior scan
  ranges — plus the FIRST non-split rational points of the carrier S12 (two certified-simple
  curves whose (3,12) is rational exactly over Q(zeta_3)); thin-set fear refuted; new
  M-square-class dichotomy {1,-3} conjectured = the localized obstruction.
- (4,16): no hit to w-height 1e5 over all 12,172 R (h<=100); 627 fibers closed unconditionally
  (rank-0 certificates); NEW two-coset descent theorem confines any hit to two explicit strata
  (one governed by congruent-number curves, empty on 2,108 fibers); two record-height simple
  [4,8] points banked (heights 4489, 7459).
- Z/35: A_1(5) chart clean to height 88 (geometric kill decay, no near-miss tail); contact lane
  CLOSED: it collapses to the irreducible genus-7 plane curve Phi38 whose 4 accessible rational
  points are all degenerate. Next: the A_1(7) mirror threefold.
- (5,5): THEOREM — the smooth F_7/F_11 slice charts are reductions of irrational points (the
  saturated slice over Q has zero rational points; orbits of degrees 28+6); h1=0 sub-branch
  closed; lane reduced to one explicit square-condition curve over a genus-0 base (genus TBD).
- (2,24): both required loci are individually DENSE (new: 2-parameter certified-simple [2,12]
  family; 37 new exact-[24] curves) yet their intersection is empty at 6,426 exact checks with
  no local obstruction — clean target for a global 2-descent incompatibility theorem.
- (2,2,4,8): sign-reduction THEOREM — every route/twist reduces to the tor2244 list, so no
  (2,2,4,8) exists in the all-squares model for ANY d <= 65535; the genus-3 split-locus curve
  is completely closed (third elliptic quotient, rank 0 unconditional); twisted slice extended
  to d<=30000 (0 hits), B=100000 resumed and running.
- (2,30): (1,1,3) route closed unconditionally (rank-0 conductor-24 quotient); target pinned to
  one explicit genus-5 gate curve PC2. Flagged [5,10] curve resolved: SPLIT (66.c3 x 66.c4).
  Banked: infinite certified-simple exact-[15] family with proven Mumford generators, and a new
  infinite certified-simple exact-[30] family (genus-0 parametrization).
- (8,8): GENUS GATE PASSED (stage-1 cover fibers genus 0 and 1, vs 21-31 in the dead chart);
  route GO. Banked: new certified-simple (2,2,4,4) member (order 64) in a provably infinite
  structured family. Remaining obstacle: one conic-type lift layer (symbolic route mapped).
- (2,6,6): split confinement now MECHANISTIC — two proven infinite rational families of split
  [2,6,6] explain the entire census; BFS sqrt3-RM escape closed (mu_3 part never rational).
  Recommend demoting below live lanes.

## Still running (detached, niced; markers RESUME_ALL_DONE in scratchpad)
- t2248prod: twisted-family enumeration to B=100000 + finalize sweep.
- t416/prod: (4,16) tier-2 w-height 1e6 sweep (all ranges redone cleanly) + postproc.

## Suggested immediate human actions
1. Add the (2,2,2,12) curve to the paper's tables (it also supersedes the '(2,2,2,12) not even
   split' line) and circulate to the team for independent confirmation.
2. Optional definitive endomorphism check (Sage geometric_endomorphism_ring_is_ZZ) — though
   LSSV already excludes QM at this torsion order.
3. The T5-pencil empirical laws (rho'-numerator = +-odd^2, denominator 240) suggest siblings:
   the bounded-descent sibling hunt is spec'd in notes/claude_prod_02_22212.md.
