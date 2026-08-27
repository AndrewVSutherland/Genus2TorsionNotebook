# Session 2026-08-14: [2,2,6,6] via the five non-identity sigma-surfaces

Executes `notes/claude_split_2266_sigma_plan_2026_08_14.md` (successor to the
08-13 session).  Environment: local Magma (license works again), <= 3
concurrent jobs, external `timeout` wall caps.

## 0. Validation gates (all PASS)

- `lane_2266_derive.m` rerun: log byte-identical to 08-13
  (`derive2266_rerun.log` = `derive2266.log`); all five anchors verified
  (ANCHOR2 prints `[ 1 ]`, not the word PASS - gate coded accordingly).
- New join code (`lane_2266_sigma.m`, Census:=1) reproduces the lane_misc2
  H=40 census EXACTLY - all 18 (sigma, drop) counts.
- Negative control (Scramble:=1, swapped probe coordinates): drop=2 rows
  collapse to 2-8 accidental matches (vs 297-1983); drop=1/3 rows land
  exactly on OTHER sigmas' counts - explained: swapping the kept key
  coordinates composes sigma with a pair-slot relabeling, which for keep
  {C1,C2} is induced by conjugation (sign-free -> another true surface),
  while for keep {C1,C3} the induced relabeling carries a sign twist
  (kernel-incompatible -> noise).  The scramble is simultaneously a
  negative control and a consistency check of the KO mechanics.

## 1. The join stage (lane_2266_sigma.m, H=150)

27425 parameter values; for every (si, drop), si in {2..6}:

- surface points (2-of-3): 3183-5328 per (si,drop), written to
  `product/data/sigma_si<si>_d<drop>_pts.txt` (t u deck columns);
- full 3-of-3 matches: 3027-5092 per (si,drop) - ALL deck pairs (j-equal);
  **genuinefull = 0 across the board**.  Consistent with the 08-13 H=200
  all-sigma sweep (17711 matches, all deck).
- swap cross-check: si5/si6 lists are (u,t)-mirrors (PASS, scanner sec. 4).

## 2. STRUCTURAL DISCOVERY: universal deck sections

`deck_derive.m`: the j-equal correspondence j(E26(t)) = j(E26(u)) factors
into the diagonal, FIVE (1,1)-components (Mobius maps), and six
(3,3)-components.  `deck_sigma_match.m`: the Mobius components biject with
the sigma-systems - each sigma-surface has exactly one as a UNIVERSAL
SECTION (all three conditions identically satisfied over Q(u)):

| sigma | universal section |
|---|---|
| id [1,2,3]      | t = u (HLP)          |
| (12) [2,1,3]    | t = (u+15)/(u-1)     |
| (13) [3,2,1]    | t = 6-u              |
| (23) [1,3,2]    | t = (5u-9)/(u-5)     |
| (123) [2,3,1]   | t = (u-21)/(u-5)     |
| (132) [3,1,2]   | t = (5u-21)/(u-1)    |

This REFINES the plan's premise: the non-id surfaces are NOT sectionless
torsor fibrations - they are elliptic surfaces with a universal section
whose points are all deck pairs (excluded as gluing candidates), but whose
section-MULTIPLES roam the fibers genuinely.  Empirically every processed
fiber has rank 1-2 (see sec. 4) - the deck section appears to be of
infinite order on essentially every fiber.

- Explains full=deckfull in sec. 1: the deck graphs lie inside the full
  3-condition locus.
- Explains the fiber population: every fiber over u0 has the deck point
  as a free rational base point.

## 3. Symbolic section-multiple machinery (lane_2266_sigma_sections.m)

The exact analogue of `lane_2266_sections.m` (which closed sigma=id), run
per sigma with base = the universal deck section, sections P1 = (-s0, Y0),
P2 = (s0, -Y0), multiples a*P1 + b*P2, a in [0..4], b in [-2..2]:

- ALL FIVE surfaces mirror the id-case exactly: the multiples produce
  condition curves of genus 7, 23, 47, 79, 119 (deg 16-240), each with
  only 0-12 small points, all degenerate; **no parametric family, no
  SECSURV**.  The distinct-multiple t(u)-values collapse to ~6 per sigma
  (effective rank 1 + torsion of the visible subgroup).
- => the deck-section-multiple route is CLOSED on all six sigma-surfaces
  (id on 08-13, the five non-id today), by the same genus-growth wall.

## 4. Per-fiber Mordell-Weil (lane_2266_sigma_fib.m)

Fibration (uniform in sigma, derived from the always-C2-linear t-side):
t = 3 + s^2/d, d = 2*A(sg[1],sg[3])(u0); org-A fiber
Y^2 = cA*(s^2+6d)(s^2-2d).  Populated by the in-process join (H=150) +
direct Points sweep (HU=40) + IsSquare(cA) infinity fibers.

DISCIPLINE (the 08-13 stall, root-caused today): `RankBounds` on the RAW
curve from `EllipticCurve(C, P0)` stalls indefinitely (probe: the same
fiber's MINIMAL model takes 0.00 s).  Also Magma's `Alarm` does NOT
interrupt deep descent - wall caps must be EXTERNAL (`timeout -s KILL`).
Final discipline: MinimalModel ALWAYS; RankBounds on the minimal model;
rank 0 -> TorsionSubgroup only; rank >= 1 -> bounded TwoDescent-cover
point search (NO unbounded MordellWeilGroup anywhere); unequal bounds ->
MWSKIP + Points(C : Bound := 1e5) fallback; coefficient box on found
generators + torsion; per-fiber side-channel; TimeCap + external timeout.

Smoke (Sig=2, 40 fibers): 0.25 s/fiber, 0 MWSKIP, 0 GENFAIL, all fibers
rank 1-2, and the third-condition survivor set = exactly the deck points
(SKIPISO) - an end-to-end positive control (the pipeline finds and
correctly dispositions the known full-condition points).

Production (H=150, org-A, all five sigmas; fibers processed multi-point
first, then by ascending height; TimeCap truncation records the frontier):

| shard | fibers (multi) | processed | rank>=1 | MWSKIP | survivors | genuine |
|---|---|---|---|---|---|---|
| si3 (13) | 4612 (151) | 2193 (TimeCap 7200s) | 2103 (96%) | 90 | 2181 = all deck | 0 |
| si2 (12) | 4114 (152) | 2397 (TimeCap 7200s) | 2314 (96%) | 83 | 2386 = all deck | 0 |
| si4 (23) | 5188 (190) | 1959 (TimeCap 7200s) | 1890 (96%) | 69 | 1943 = all deck | 0 |
| si5 (123) | 3072 | 1831 (TimeCap 3600s) | 1802 (98%) | 29 | 1818 = all deck | 0 |
| si6 (132) | 3094 (147) | 1884 (TimeCap 3600s) | 1839 (98%) | 45 | 1866 = all deck | 0 |

TOTALS: 10 264 fibers processed with the bounded-MW machinery across the
five surfaces (multi-point fibers all covered - they were processed first;
the TimeCap frontiers cut only the single-deck-point tails at large
height(u0)); 96-98% of fibers have rank >= 1 and rank0 = 0 everywhere -
the deck section has infinite order on essentially every fiber; 316 MWSKIP
(unequal 2-descent bounds, ~3%; Points fallback applied), 0 GENFAIL;
10 194 third-condition survivors and EVERY ONE is a deck point; genuine
candidates: ZERO.  (Org-B (fold C3) pass deferred: with org-A uniformly
deck-confined on all five surfaces, the independent org-B boxes are the
cheapest unexplored coverage but were left as future work - see sec. 10.)

COMPLETION PASS (codex round 3, `lane_2266_sigma_fib2.m`, 2026-08-15): the
original lane capped generator collection at TWO unverified points, so the
6444 rank >= 2 fibers had potentially incomplete boxes (certainly so for
the 1462 with rank 3-4).  All 6444 were reprocessed with generators drawn
from ALL TwoDescent covers and a VERIFIED-INDEPENDENT subset selected by
height-pairing regulator (> 1e-6), boxes over up to min(rank, 4)
generators (NB 10/10/5/4):

| shard | fibers | complete | GENSHORT | survivors | genuine |
|---|---|---|---|---|---|
| si2 | 1673 | 1398 | 274 | 1658 = all deck | 0 |
| si3 | 1295 | 1164 | 131 | 1286 = all deck | 0 |
| si4 | 1141 | 1006 | 123 | 1117 = all deck | 0 |
| si5 | 1264 | 1118 | 145 | 1253 = all deck | 0 |
| si6 | 1071 |  992 |  66 | 1035 = all deck | 0 |
| total | 6444 | 5678 (88%) | 739 | 6349 = all deck | 0 |

GENSHORT = the cover-point search (Bound 3000, <= 12 candidates) could not
supply rank-many independent directions - those fibers are honestly
INCOMPLETE (u-lists in the logs), not silently truncated.  Survivors remain
100% deck; the deck-confinement statement now rests on verified-complete
boxes for 88% of the rank >= 2 fibers and >= 2 independent directions on
the rest.

NOBASE redo (codex round 4): 25 fibers (12 si4, 13 si6) whose only points
came from the original direct sweep had no base point in the completion
pass (which only rebuilt the join) and were skipped while still counted as
processed.  Fixed (bounded Points fallback for the base) and reprocessed
via `fib2_usets_nobase.m`.

SATURATED RERUN (codex round 5, `sigfib3_si*.log`): the round-3 pass used
verified-INDEPENDENT generators but never verified SATURATION - a
finite-index sublattice's box misses whole cosets.  Final pass over all
6444 rank >= 2 fibers with Saturation(gens, 100) after the greedy
selection, plus conclusive handling of unequal rank bounds (rank pinned
to rhi when rhi independent points are exhibited, else RBOPEN):

| class | count | meaning |
|---|---|---|
| complete | 5701 (88.5%) | full box over a p<=100-SATURATED full-rank basis (see caveat) |
| GENSHORT | 740 | fewer independent points than the (conclusive) rank - incomplete |
| RBOPEN | 3 | rank bounds [1,2] unresolved (si2 u=-47/73, u=-65/4; si5 u=-65/4) |
| nobase/mwskip/crashskip | 0 | - |

All 6444 accounted (`fib3_account.py` stitches per-fiber dispositions
across crash-resumed segments, deduping survivors by (t,u) since resumed
segments overlap a few fibers).  ~~KEY EMPIRICAL FACT: satfix = 0 on every
shard~~ **RETRACTED (2026-08-16, full-saturation certificate pass): the
satfix = 0 outcome was an ARTIFACT of a bug in the saturation-adoption
guard - whenever Magma's Saturation returned extra points (it habitually
appends a torsion generator), `#gsat eq #gens` failed and the result was
silently discarded, so saturation was effectively never verified on most
fibers.  The certificate pass (below) found REAL enlargements: e.g. si5
u=2 has exact index 2 (regulator ratio 4.000..., old gens = (1,0), (-1,2)
in the saturated basis).  Affected fibers' boxes enumerated a sublattice
and their box claims are void until the corrected rerun (SATFIXBOX list)
lands - see the certificate-pass block in sec. 12.**  Saturation caveat
(codex round 6): Saturation(gens, 100) certifies p-saturation only for
p <= 100; the full-saturation certificate (Minkowski + CPS/Silverman +
bounded search + Saturation-to-Bidx) is the 2026-08-16 pass of sec. 12.

**CORRECTED BOX FLEET (2026-08-17, `sigfib4_si*.log`) - the missed-coset
hole CLOSED**: all 6444 rank >= 2 fibers re-enumerated with
verified-saturated bases (SatAdopt at p <= 100 with ADOPTION):

| shard | fibers | complete | genshort | rbopen | satfix | survivors | genuine |
|---|---|---|---|---|---|---|---|
| si2 | 1673 | 1398 | 272 | 3 | 920 | 1673 = all deck | 0 |
| si3 | 1295 | 1162 | 133 | 0 | 830 | 1295 = all deck | 0 |
| si4 | 1141 | 1018 | 123 | 0 | 769 | 1141 = all deck | 0 |
| si5 | 1264 | 1120 | 144 | 0 | 786 | 1264 = all deck | 0 |
| si6 | 1071 | 1005 |  66 | 0 | 765 | 1071 = all deck | 0 |
| total | 6444 | 5703 (88.5%) | 738 | 3 | **4070 (63%)** | 6444 - 100% deck | **0** |

63% of the rank >= 2 fibers had genuinely enlarged (typically index-2)
lattices - the round-5 boxes were sublattice-short at scale - and with
every previously-missed coset now enumerated, the survivors are STILL
100% deck points and there are ZERO genuine candidates.  Precisely
scoped (codex): the deck-confinement statement is RESTORED on corrected,
verified-saturated lattices FOR THE 5703 VERIFIED-COMPLETE FIBERS; the
738 GENSHORT + 3 RBOPEN fibers ran partial boxes and remain flagged
incomplete as tabulated.  (`sat_usets.m` regenerated from the sigfib4
dispositions - the certificate fleet targets exactly this corrected
5703-fiber set.)  Survivors 6384 distinct (t,u) pairs -
all deck; genuine instances 0.  Operational: Magma segfaulted
state-dependently once (si2, fiber 737, ~700 descents into the process;
the same fiber passes in a fresh process) - handled by `fib3_driver.sh`
(crash-resuming wrapper: retry once fresh, CRASHSKIP on a double crash;
none needed a skip).

## 5. Empirical rational-curve structure + line sweeps

`sigma_section_scan.py` (exact-arithmetic scanner; selftest with planted
sections PASS) on the H=150 lists found, per sigma, dozens of rational
curves carrying 5-56 surface points (counts as corrected in codex round 2,
sec. 11 - the original scan let common num/den roots inflate deg-2 support):

- si3 carries t = u+4 and t = u-4 (56 points each): on t = u+4 condition
  C3 is IDENTICALLY square ( -(t-1)(t-9)*-(u+3)(u-5) = (u+3)^2(u-5)^2 ),
  and the remaining system collapses to two CONIC conditions - the full
  [2,2,6,6] system on this line is a genus-1 problem in u alone.
- si3 also shows the lambda-ruling t = 3 - k^2(u-3) (C2 identically
  square for every slope -k^2) - the second elliptic fibration of the
  surface; top representatives carry 13-16 points.
- si3/si5/si6 carry deg-2/deg-2 curves with 38-56 points; si2/si4 carry
  a swap-pair of lines t = (7/15)u + 8/5 / t = (15/7)u - 24/7 (6 each).
- si5/si6 share their lines (swap symmetry), as do si2/si4.

`lane_2266_sigma_lines.m`: for 48 such curves, symbolic restriction of
the three conditions (RESIDUAL classes printed; PARAMETRIC-FULL-FAMILY
detector armed) + exact u-sweep to height HSweep=2000 (13x the join
height) + funnel of every full pass.  RESULTS (458 s):

- Magma's exact simplification exposed several scanner deg-2 candidates as
  Mobius curves in disguise: si3's three deg-2 curves are t = u+4 (again),
  t = (6u-9)/u, t = -9/(u-6); si5's two are both t = (21u-9)/(5u-9).
- No PARAMETRIC full family; no new full-condition points: across all 48
  curves the sweep found exactly 26 LINEFULL passes, and EVERY one lies on
  the ambient sigma's universal deck section (verified pointwise, e.g.
  si2 (7/3, 13) = ((13+15)/(13-1), 13); si5 (125/49, -101/19) on
  t=(u-21)/(u-5)) - all SKIPISO, funneled=0.
- Residual structure recorded: on the marquee curves one condition is
  identically square and the residual pair has degrees [2,2], [1,2] or
  [2,1] - genus <= 1 reductions (closure certificates: sec. 5b below);
  the lambda-ruling lines have residual degrees [4,4]; generic scanner
  lines have [4,3,3] (no trivial condition - not special, just populous).

**Pattern (multi-method, now strong): every full-condition rational point
found by ANY method - join to H=200 (08-13), sigma-join to H=150, fiber MW
boxes, 48-curve sweeps to u-height 2000 - lies on one of the six universal
deck sections.  Conjecturally the rational points of the full [2,2,6,6]
locus are EXACTLY the deck graphs.**

### 5b. Closure certificates for the reduced lines: ALL EMPTY at p=3
`lane_2266_line_close.m`: for each special curve (one condition identically
square, residual pair of degrees <= 2), parametrize the first residual conic
and substitute into the second -> an integral quartic W^2 = g(U).  The nine
curve entries collapse to THREE distinct quartics:

    t=u+4 (si3):                       W^2 = -U^4-4U^3-118U^2-228U-4209
    t=u-4, -9/(u-6) (si3);
      (5u+11)/(u-1) (si6, from all
      three deg-2 scanner curves):     W^2 = -15U^4-34U^2-15
    (6u-9)/u (si3);
      (21u-9)/(5u-9) (si5, both
      deg-2 scanner curves):           W^2 = -393216U^4-...-3840

and ALL THREE are LOCALLY INSOLUBLE AT p=3 (checked 2 <= p <= 100; no
rational points to 1e4 as a sanity cross-check).  **Every special Mobius
curve on every sigma-surface is provably empty - by a pure p=3 local
obstruction** - while the ambient surfaces are everywhere-locally-soluble
(sec. 7).  The [2,2,6,6] condition is 3-adically obstructed on exactly the
strata where it is "almost" satisfied; note p = 3 is the degeneration prime
of the E26 family (t = +-3 excluded), which is presumably the source.
Together with sec. 3 (deck-multiple condition curves genus >= 7, dry) and
the fiber MW sweeps, every identified rational curve on the surfaces is
now either the deck section or provably [2,2,6,6]-empty.

## 6. Secondary lane: congruence sieve (CLOSED - dry)

`lane_cong_sieve.m` with Fam:="26", H:=150: pool = 21658 E26 curves
(j-deduped), all ~2.3*10^8 pairs sieved mod 5 AND mod 7 (MINCMP 14):
**SIEVE_DONE cands=0** (523.9 s) - not a single congruent pair, isogenous
or not, in the whole E26 family to height 150.  The genuine-5/7-congruent
[2,6]x[2,6] route has NO instances at family scale; combined with the
08-13 closure of the (N,N)-congruence route this lane is exhausted.
(`lane_sigma_strictglue.m` with the strict rationality ladder was built
and is ready, but has no targets; kept for future use.)

## 7. Local analysis (thin vs obstructed): ELS

`lane_sigma_local.m`: for every sigma and every place p <= 47 and R, the
full 3-condition variety has a Q_p-witness (16/16 places soluble per
sigma, instant).  No local obstruction anywhere: a dry global outcome
means the rational points are GENUINELY THIN (boundary/deck-confined),
not locally impossible.

## 8. Cleanup theorem (twisted diagonals): UNCONDITIONAL

Three-stage finish of session-note-08-13 item 2.4 (`lane_tw_cleanup.m` /
`_cleanup2.m` / `_cleanup4.m`):

1. **Rank 0 without GRH**: RankBound(J) = 0 for both genus-3 Jacobians with
   default (unconditional) class-group bounds - the GRH caveat of
   lane_misc2.m is GONE.
2. **Exact 2-torsion** from the odd-degree-model factorizations:
   tw12 factors [1,2,2,2] -> J(Q)[2] = (Z/2)^3; tw10 factors [1,1,1,2,2]
   -> J(Q)[2] = (Z/2)^4.  (gcd bounds over 93 primes stall at 128 resp. 32
   with min v2 = 7 resp. 5 - systematic 2-adic slack, so gcd alone cannot
   finish; DivisionPoints is not implemented in genus 3.)
3. **Class-order sieve**: reduction J(Q)_tors -> J(F_p) is injective for
   odd good p, so ord[P - infinity] is preserved exactly.  Over 6 primes
   per curve, the set of pure-2-power orders realized by points of C(F_p)
   intersects to {1, 2} for BOTH curves => every rational class has order
   <= 2.  The Mumford representative of [P - infinity] for affine
   non-Weierstrass P is (x - x_P, y_P), y_P != 0, which is never
   self-inverse - so no such P exists, and
   **C(Q) = {infinity} + rational Weierstrass points** for both twisted
   diagonals: the known degenerate points EXHAUST C(Q), unconditionally.
   The [10,10]/[12,12] conjugation-self-gluing route is closed as a
   THEOREM (no GRH, no height cutoff).

## 9. Operational lessons recorded

- MinimalModel before ANY descent call (RankBounds on raw glued models
  stalls; minimal model instant).
- Magma Alarm() does not fire inside deep C descent code - use external
  `timeout --signal=KILL`.
- Never call unbounded MordellWeilGroup in a fiber loop - TorsionSubgroup
  + RankBounds + bounded TwoDescent-cover search covers every case with
  strictly bounded work.
- The scramble negative control doubling as a KO-mechanics consistency
  check (sign-twisted vs sign-free relabelings).

## 10. Frontier / next session

**Closed this session (do not repeat):**
- The five non-id sigma-surface fibrations: deck-section multiples (genus
  wall), per-fiber MW org-A boxes (H=150 joins + HU=40 sweeps + NB=10
  boxes), 48-curve sweeps to u-height 2000.
- The nine special Mobius curves - provably empty (p=3).
- E26-family internal 5/7-congruence (0 pairs to H=150).
- Twisted diagonals - unconditional theorem (sec. 8).

**Structural facts to build on:**
- Universal deck sections delta_sigma (sec. 2 table) - the surfaces are
  positive-generic-rank elliptic surfaces; deck confinement of the full
  locus is uniform across every probe.
- Each surface has (at least) three elliptic fibrations: the u-fibration
  (C2-parametrized, swept here), the lambda-ruling (C2-identical lines,
  residual degrees [4,4] - NOT systematically MW-attacked), and the
  condition-identical Mobius pencils (the 1-parameter families through the
  special curves; only the rational-scalar members were closed).
- The deg-(3,3) deck correspondences (deck_derive.m) were not exploited
  (non-rational sections; their rational points are deck pairs anyway).

## 12. Full-saturation certificate pass (2026-08-16, owner-approved)

`lane_2266_sigma_sat.m` over the 5701 "complete" fibers.  Per fiber:
verified saturation-adoption at p <= 100 (SatAdopt - the fix of the
silent-skip bug retracted in sec. 4), then a capped exhaustive search
below Teff := log(Bsearch) - ub (ub := min(CPS, Silverman) upper bound
for h - hhat; conventions calibrated in probe_ht.m: Magma hhat ~ log H(x),
CPS bounds h - hhat, Points' Bound caps H(x)); Minkowski turns the proven
m1 lower bound into an index bound Bidx := sqrt(R_sub gamma_r^r / m1^r),
and Saturation(gens, Bidx) then excludes every possible prime factor of
the index: idx = 1 PROVEN (SATCERT / SATCERT2).  Enlargements are adopted
and flagged SATFIXBOX (their round-5 boxes ran on a proper sublattice -
rerun with lane_2266_sigma_fib2.m, which now also adopts saturation before
enumerating).  Fibers where even the hard search cap cannot reach
Teff > 0 (CPS/Silverman ub too large) are SATOPENB-final with a plain
Saturation pushed to 20000 recorded.

RESULTS (CLOSED OUT 2026-08-18 by owner decision, fleet stopped partway):
532 of the 5703 corrected-complete fibers carry full idx = 1 certificates
(21 SATCERT + 511 SATCERT2; si2 183/1398, si3 205/1162, si4 144/1018
processed prefixes; si5/si6 not started), 26 SATOPENB, 1 skip, and -
importantly - ZERO SATVIOLATIONs: on every fiber checked, the certificate
machinery confirmed the box lattices rather than contradicting them.

**FINAL STATUS OF THE [2,2,6,6]-VIA-E26xE26 ROUTE (owner's framing, for
the record): the paper will make positive claims only, so this chart's
closure is recorded as evidence, not as a theorem.  We have NOT completely
ruled out this path to a realization - the residual possibilities are
(i) an index with all prime factors > 100 on the ~91% of complete fibers
without an idx = 1 certificate, (ii) the 738 GENSHORT + 3 RBOPEN partial
fibers, (iii) box reach beyond NB = 10, heights beyond the sweeps, and
other fibrations/organizations - but after the joins, the corrected
coset-complete boxes, the sweeps, the section machinery, the sieve, and
the local analysis all came back deck-confined, a realization through
this chart now seems LESS LIKELY THAN OTHER POSSIBILITIES (the sec. 5 gap
routes: (3,3)-gluings, section-gain constructions, other charts).**

**Plausible next moves (in decreasing structure):**
1. The lambda-ruling as a FIBRATION: fibers over lambda are genus-3 (two
   quartic conditions) - but special lambda values may degenerate; a
   discriminant analysis could find the lambda where the fiber drops genus.
2. Org-B (Org:=3) fiber pass - independent boxes on the same surfaces
   (deferred; org-A uniform deck-confinement makes a different outcome
   unlikely but it is the cheapest unexplored coverage).
3. Prove deck-confinement: the full 3-condition variety is a
   (Z/2)^2-cover of each sigma-surface; the p=3 mechanism that kills the
   special curves suggests a uniform 3-adic argument on the complement of
   the deck sections - if it works, [2,2,6,6] via E26xE26 (2,2)-gluing
   would be IMPOSSIBLE outright, a clean theorem.
4. Leave the chart: [2,6,12] / [12,12] via (3,3)- or higher gluings (the
   gap classification in split_torsion_table.md sec. 5), or entirely
   different 144-decompositions.

## 11. Codex review rounds (PR #16, 2026-08-15)

Round 1 (5858ff9): [8,8] family count corrected 86 -> 59 distinct (the
H=120 instances are a subset of the H=250 log's 59 - verified by t-value
comparison); scanner selftest fixed for the supported no-numpy mode.

Round 2: three findings, all confirmed against the tree, one of which
unraveled a false belief:

- **The analytic-glue "validation" of 08-13 was itself a weak-test
  artifact.**  Investigating the P1 on `analytic_glue.m`'s weak RatApprox:
  the old `validate_glue.log` GLUE line (prec 120) printed a height-47
  "rational invariant" for M=[0,1;1,0], and the log ends right there (the
  run stalled in Mestre before ever comparing with the algebraic
  gluings).  Probes at prec 120/200/300/400/500 show the pipeline is fully
  precise (stage drift 10^-198..10^-201 between prec 200 and 500) and the
  analytic j1 sits EXACTLY 10^-80 from that height-47 rational at every
  precision - the printed value was a best-approximant artifact.  The
  recognition-free comparison then showed the truth: the six analytic
  anti-isometries match the six algebraic gluings BIJECTIVELY to
  10^-175..10^-191, with true invariant heights 100-182 digits.
  Fixes: `StrictRatApprox` (the glue_window height-relative ladder,
  err < 10^-(2*height+15), prec-driven rungs) is now the only recognition
  path in GlueScan/GlueCandidates (weak RatApprox kept only as a warning
  exhibit); `validate_glue.m` rewritten as a two-gate validation -
  (1) recognition-free bijective match at prec 200 (6/6 PASS at
  10^-175..191), (2) strict GlueScan at prec 600 returns exactly the six
  algebraic rational triples (6/6 PASS, heights to 182).  The analytic
  gluer is hereby PROPERLY validated for the first time; every
  weak-recognition output in older logs (isoglue2 SCANINV lines, the old
  validate_glue GLUE line) is superseded.  Note the practical corollary:
  genuine glued invariants at these conductors have heights 100-182, so
  strict recognition needs prec >= ~400 and Mestre at such heights is the
  documented stall - GlueScan + invariant-identity is the right
  verification pattern, not curve reconstruction.
- Scanner deg-2 support counting let common num/den roots (poles) accept
  arbitrary t at that u; fixed with a den != 0 guard; supports drop
  58->56, 50->48, 49->48, 40->38, 39->38, 54->52 (candidate set and all
  conclusions unchanged; selftest PASS in both numpy modes).
- Four 08-13 scripts (classsweep/degsweep/pin90/verify600) hardcoded
  ChangeDirectory("/home/claude/torsion_jac/product/code"); replaced by
  the house run-from-product/code convention with an optional
  TORSION_JAC_ROOT override.

Round 3: two P1s.
- **Fiber generator completeness (confirmed - the sharpest catch of the
  review).**  The original fiber lane capped cover-point collection at two
  unverified generators; 1462 fibers with rank 3-4 (and 4982 rank-2 fibers
  with unverified independence) therefore had potentially incomplete
  coefficient boxes.  Remedy: `lane_2266_sigma_fib2.m` completion pass
  over all 6444 rank >= 2 fibers (u-lists parsed from the committed logs
  into `product/data/fib2_usets.m`) with verified-independent generator
  sets (height-pairing regulator) and rank-adaptive boxes; results in the
  sec. 4 completion table - 88% verified-complete, 739 GENSHORT flagged
  incomplete, survivors still 100% deck, zero genuine.  The
  deck-confinement statement was re-qualified accordingly.
- **StrictRatApprox hq (refuted - already max(num, den) as committed);
  but the described defect class existed mirror-imaged in two 08-13
  scripts** (glue_window.m and pin90.m used denominator-only hq, so an
  integer-like q would face only a 10^-15 bar); both fixed to the max.

Round 4: two P2s, both confirmed - the 25 sweep-populated NOBASE fibers
(recovered, see sec. 4) and the stale GRH wording in the 08-13 report's
twisted-variant row (updated to the unconditional closure).

Round 5: two comments, both acted on - saturation of the generator sets
(P1, confirmed in principle; the "satfix = 0" outcome later RETRACTED as
a verification-guard artifact - see sec. 12) and the two unaccounted
MWSKIP fibers (P2,
confirmed; the rerun's conclusive-rank handling resolves one class and
the final accounting lists 3 RBOPEN fibers explicitly - no fiber is
silently absent from the totals).

