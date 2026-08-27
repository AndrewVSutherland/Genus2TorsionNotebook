# Top-10 ranked list: torsion subgroups not yet realized by a geometrically simple genus 2 Jacobian /Q

Date: 2026-07-17. Synthesized from: NotesAndTodo.tex, main.tex, paper/g2torsion.tex (Sutherland's
table, PDF dated 2026-02-18), all notes/*.md, code/ + data/ artifacts, an LMFDB census
(g2c_curves split by is_simple_geom), and a literature survey (Elkies 2002/2024, HLP 2000,
Howe 2015, Nicholls 2018, Platonov et al., LSSV QM-torsion bound, 2024-26 arXiv scan).

## Baseline (known geometrically simple, incl. project banked wins)
Cyclic 1-30, 32, 33, 34, 36, 39, 40; (2,2n) n <= 11, 13, 14; (2,2,2n) n <= 7 plus (2,2,20);
(2,2,2,2n) n <= 5 (incl. (2,2,2,8) infinitely, (2,2,2,10) Elkies 2024); (3,3), (3,6), (3,9);
(4,4), (2,2,4,4) infinitely, (2,4,4) infinitely, (4,8), (2,4,8), (6,6).

## Ranking criteria
A. Existence evidence: realized by SPLIT genus-2 Jacobians (moduli has rational points;
   only geometric simplicity missing) > realized only by products of elliptic curves > nowhere.
B. No proven/suspected obstruction in the explored charts (local go/no-go verdicts).
C. Concrete unfinished machinery in this repo (tractability of the next experiment).
D. Prize value (record order, conspicuous gap in a realized chain).

## The ranked ten

1. **(3,12)** [order 36 - smallest open group anywhere]. Split realizations abound (HLP P^2
   family; three split in-project hits e.g. M(2,12) at (z,r,a)=(-5/3,-3/5,-10/9), splits as
   90c3 x 510g1). (3,3),(3,6),(3,9) all simple-realized. No obstruction theorem; the carrier
   surface S12 over the M(2,12) chart is irreducible with UNDETERMINED geometry. Fresh routes:
   S12 type determination; M(12)+cubic-contact-3 (pure, without the extra 2); Bruin-Flynn-
   Shnidman sqrt3-RM full-level rational surface (arXiv:2102.04319).
2. **(2,24)** [order 48 - conspicuous gap: (2,22),(2,26),(2,28) known]. Split (2,24) = HLP P^2
   family. Project has TWO rational one-parameter simple [24] families (main M(24), proven
   [2,24]-free; new G1 family z=(1+t^2)/(2t), r=(t^2-1)^2/(8(1+t^2)), a=-(t^2-1)^2(1+t^2)/
   (2t^2(t^2+3)^2)) plus 25 unfitted off-main halving points = more components. G1 Mechanism B
   (rational quadratic factor of Q4) NEVER RUN - the single most concrete unfinished step.
3. **(4,16)** [order 64]. Split realization known (Sutherland split-not-product list). Route
   revived: NO local obstruction at 7 (218 ALIVE 7-adic reps); rigorous negative only to
   height 800 on M_1(8,4); 52 in-family near-misses "locally compatible, globally killed".
   Unlocks: spec'd-but-unbuilt per-R elliptic solve for Delta=square (replaces O(H^4) loop);
   also the unfinished rational-elliptic-surface parameterization inside M_1(8,2^w).
4. **(2,2,2,12)** [order 96 - not known even split; products only]. But the project's
   best-diagnosed target: forced bad reduction at 3,5,7,11,13 (15015 | N explains database
   absence); 13-adic GO (3552 Hensel-liftable branches); no obstruction p >= 17; production
   enumerator code/tor22212.c ready, NEVER run at scale (~260 expected local survivors at
   H=10^4). Caveat: independence heuristic says a hit needs an algebraic subfamily.
   Complementary: M(2,2,2,6) subsumption sweep clean to H=1200 with 64 T3/T5 near-misses
   awaiting a rank-0 kill.
5. **(2,6,6)** [order 72]. Split family = HLP positive-rank elliptic surface. The contact-6
   chart ALREADY produces [2,6,6] curves - just all Q-isomorphic to one non-simple curve up to
   height 60. Torsion structure achieved; escaping the split locus is the residual problem.
   Height-150 chart scan costs only ~10 CPU-min. (6,6) simple was just realized on this chart.
6. **(2,2,4,8)** [order 128 - the flagship; would tie the all-time genus-2 record and beat
   Elkies' simple record of 80]. Split = HLP positive-rank elliptic surface (10 models stored,
   heights 10^19+). Simple locus: exhaustive negative over all 30,387 (2,2,4,4) primitives
   d <= 65535 (all components/twists); K3 route killed by signed local analysis at 11/23 +
   real place. LIVE routes: the HPL (Z/2)^4-cover threefold (fibration-finiteness does NOT
   apply there); the never-enumerated twisted-(2,2,2,8) second component intersected with
   (2,2,4,4) conditions; the genus-3 split-locus curve {y^2=(u-3)(u+1)(u^2-6u+1),
   z^2=-(u-1)(u^2-6u+1)} (rank-1 quotients, undecided).
7. **(2,2,16)** [order 64]. Split (non-product) realization known. Norm-surface route shows a
   0/69 global-only failure pattern (suspected Brauer-Manin on the 2-cover; d=-2 least
   obstructed) - worth settling either way. UNEXPLORED alternative: second-halving tower over
   the newly discovered twisted (2,2,2,8) family (52 seeds, e.g. (2,3,12,18); simple member
   (29,121,125,145)).
8. **(8,8)** [order 64]. Split = HLP P^2 family. NO local obstruction through p=43 in the
   M_1(8,4) chart, but the second-halving locus there is a general-type surface (fiber genera
   23/31) - dry by structure, not by accident. Needs a construction with two independent
   8-chains: Nicholls' 3-parameter pointwise-rational (4,4)-kernel family (Prop 5.9.6) +
   double 8-divisibility is the untried substrate.
9. **(2,30)** [order 60 - products only; not even a split genus-2 Jacobian known]. C30 simple
   is known (Nicholls; plus the project's infinite contact-5/6 [30] family, whose extra-2
   branch is closed by theorems). The team's own June 26 idea - combine the A(6)/M(6)
   6-torsion construction with Elkies' 5-torsion threefold (the (2,2,2,10) engine) - is the
   designated fresh route and has never been set up.
10. **(35)** [cyclic - the most accessible of the cyclic gaps 31, 35, 37, 38]. Split C35
    exists (HLP positive-rank elliptic curve family). Contact-5+7 / contact-7+5 charts are
    obstructed mod 3 (and 5), but those are M-type (point-on-curve) charts only; Elkies'
    LuCaNT A_1(5) threefold (all 5-torsion classes, not just P-infinity ones) intersected
    with a 7-torsion condition (NDE-2003 7-torsion family) is strictly bigger and untried.

Near-list (honorable mentions): (2,2,24) [split-realized, but both in-repo routes strongly
obstructed], (6,12) [strictly above (3,12) - do (3,12) first], (2,2,18) [not even products],
cyclic 31/37/38 [no known construction anywhere], (5,10), (4,12), (2,2,2,14).

## Cross-cutting strategy lessons (from the repo's own history)
- Every realized group came from constructions where ALL the torsion is forced by algebraic
  identities; every "family + necessary-condition sieve" attempt died at a small-prime
  boundary that is Hensel-live locally but globally dry.
- Richelot/2-power isogeny neighborhoods always LOSE torsion.
- The QM locus is capped at torsion order <= 18 (LSSV) - all hunts must live in RM or End=Z.
- Palindromic sextics are never simple; check simplicity certificates early
  (irreducible L_p + irreducible 12th-power transform; Nicholls Prop 2.4.2).
