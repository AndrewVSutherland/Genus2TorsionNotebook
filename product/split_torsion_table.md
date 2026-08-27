# Torsion subgroups of split Jacobians of genus-2 curves over Q

**Master table, assembled 2026-08-12.** "Split" means the Jacobian is *not
geometrically simple*, i.e. Qbar-isogenous to a product of two elliptic curves.
A group is listed when an explicit genus-2 curve C/Q is known with
J = Jac(C) geometrically split and J(Q)_tors *exactly* the given group.

**Sources.**
- **DB**: the LMFDB alpha genus-2 database `g2c_curves_new` (6,216,959 curves,
  of which 3,415,569 geometrically split, 3,377,845 split over Q), queried via
  the LMFDB MCP on 2026-08-12.  Per-group counts and minimal-conductor split
  witnesses below come from these queries.
- **Literature**: verified sweep in `product/literature_split_torsion.md`
  (every claim re-checked against the source paper and re-verified in Magma).
  Key: HLP = Howe–Leprévost–Poonen, Forum Math. 12 (2000); Howe15 = Howe,
  Bull. LMS 47 (2015); Lep95 = Leprévost, JTNB 7 (1995); Ogg73 = Ogg (PSPM 24).
- **repo**: exact split witnesses recorded earlier in this repository
  (lane-8 anchor bank `results/claude_ov_lane8_verify.log`, 2026-07-25; the
  `[2,2,4,4]` source-51 audit; contact-6/M(2,12) split banks).
- **new**: witnesses constructed by this session's gluing lanes
  (`product/code/lane_*.m`, logs in `product/logs/`, verification in
  `product/logs/verify_witnesses*.log`, curves in
  `product/data/new_split_witnesses.txt`).

## 1. The 75 groups known to occur (plus the trivial group)

Sorted by order then invariants.  "DB n" = number of geometrically split
curves in the alpha DB with exactly this torsion; witness = minimal-conductor
DB example (LMFDB label) or pointer to a construction.

| group | order | DB n | witness | notes / literature |
|---|---|---|---|---|
| [] | 1 | 50559 | 961.a2 |  |
| [2] | 2 | 674171 | 336.a1 |  |
| [3] | 3 | 26703 | 324.a2 | X0(37) is a split [3] (lit) |
| [4] | 4 | 174341 | 289.a2 |  |
| [2,2] | 4 | 1727196 | 256.a2 |  |
| [5] | 5 | 6517 | 121.a2 |  |
| [6] | 6 | 119916 | 196.a3 |  |
| [7] | 7 | 1083 | 676.b3 |  |
| [8] | 8 | 26187 | 300.a2 |  |
| [2,4] | 8 | 471270 | 225.a1 |  |
| [2,2,2] | 8 | 13681 | 315.a2 |  |
| [9] | 9 | 386 | 722.a1 |  |
| [3,3] | 9 | 735 | 361.a1 |  |
| [10] | 10 | 4913 | 363.a2 |  |
| [12] | 12 | 7387 | 294.a1 |  |
| [2,6] | 12 | 76505 | 196.a2 |  |
| [14] | 14 | 19 | 1352.a1 |  |
| [15] | 15 | 159 | 484.a1 | X0(50) (lit) |
| [16] | 16 | 154 | 2700.d1 | **extra** (not ⊂ any product) |
| [2,8] | 16 | 12523 | 225.a8 |  |
| [4,4] | 16 | 8350 | 225.a6 |  |
| [2,2,4] | 16 | 5544 | 225.a3 |  |
| [2,2,2,2] | 16 | 609 | 784.c1 | HLP [2,2,4,8]-family contains it |
| [18] | 18 | 27 | 490.a1 |  |
| [3,6] | 18 | 1039 | 196.a4 |  |
| [19] | 19 | 1 | 169.a1 = X1(13) | **Q-simple only**; Mazur–Tate 1973; extra |
| [20] | 20 | 267 | 363.a3 | HLP (P²-family) |
| [2,10] | 20 | 79 | 256.a1 = X1(16) | min ex. Q-simple (Ogg73); Q-split exist |
| [21] | 21 | 39 | 324.a1 = X1(18) | min ex. Q-simple; Q-split: Lep95 C21 ~ J0(26) |
| [24] | 24 | 295 | 300.a1 | Lep95 C24,2 split |
| [2,12] | 24 | 3516 | 450.a2 |  |
| [2,2,6] | 24 | 274 | 600.a1 |  |
| [25] | 25 | 1 | 42336.dw1 | **Q-simple only** (Lep95 C25, split per PZP 2013); extra |
| [5,5] | 25 | 7 | 121.a1 | X0(22) (HLP Cor. 5) |
| [27] | 27 | 5 | 7776.e2 | Lep95 C27 split/Q; extra |
| [3,9] | 27 | 21 | 2916.e2 | HLP (P²) |
| [28] | 28 | 3 | 4732.c1 | Howe15 (3,3)-gluings, 5 curves |
| [2,14] | 28 | 3 | 86436.z3 |  |
| [30] | 30 | 88 | 4356.j1 | HLP §3.5 |
| [2,16] | 32 | 30 | 4608.h5 | **extra** |
| [4,8] | 32 | 366 | 225.a2 |  |
| [2,2,8] | 32 | 227 | 360.a1 |  |
| [2,4,4] | 32 | 21 | 225.a5 |  |
| [2,2,2,4] | 32 | 27 | 882.a5 |  |
| [35] | 35 | 6 | 7436.a2 | HLP §3.6; repo lane8 exact anchor |
| [36] | 36 | 3 | 4860.f1 | = P&P 2012's f36 (models: Platonov RMS 69:1 (2014) §6; non-simplicity proven there §7; Magma-verified 2026-08-13, `logs/verify_pp_curves.log`) |
| [2,18] | 36 | 1 | 64980.bo1 |  |
| [3,12] | 36 | 18 | 2700.c1 | HLP §3.5; repo M(2,12) hits |
| [6,6] | 36 | 176 | 196.a1 = X0(28) | HLP §3.4; repo split bank |
| [40] | 40 | 9 | 3168.i1 | HLP §3.5 |
| [2,20] | 40 | 3 | 30492.q2 |  |
| [2,2,10] | 40 | 1 | 30492.q1 |  |
| [45] | 45 | 0 | repo lane8 anchor | HLP §3.6 (pos.-rank elliptic curve) |
| [48] | 48 | 5 | 1764.a1 | = P&P 2012's f48,1 = (x²−4)(x⁴−10x²−3) (earliest [48]); f48,2 = 5292.c2; models Platonov RMS 69:1 (2014) §6, split §7 (Magma-verified); Howe15 infinite family |
| [2,24] | 48 | 76 | 1440.c2 | HLP (P²); repo lane8 |
| [4,12] | 48 | 2 | 8100.j1 |  |
| [2,2,12] | 48 | 4 | 1350.d1 |  |
| [2,2,2,6] | 48 | 7 | 900.b2 |  |
| [7,7] | 49 | 0 | HLP single curve | exact (HLP §3.6); repo lane8 |
| [5,10] | 50 | 1 | 85500.cj1 | HLP §3.5; repo lane8 |
| [60] | 60 | 1 | 13500.r4 | HLP §3.5 |
| [2,30] | 60 | 2 | 13068.f2 |  |
| [63] | 63 | 0 | HLP single curve | exact (HLP eq. (4)); repo lane8 |
| [4,16] | 64 | 1 | 50400.dj1 | repo palindromic curve; extra |
| [8,8] | 64 | 0 | **new** (W1, W2) | HLP §3.4 containment; **first exact witnesses this session** |
| [2,2,16] | 64 | 4 | 69300.bu1 | **extra** |
| [2,4,8] | 64 | 1 | 573300.ir1 | HLP §3.4; repo tor2244 rows |
| [2,2,2,8] | 64 | 1 | 4410.d2 | repo bielliptic tuples |
| [2,2,4,4] | 64 | 0 | repo (source-51 audit) | exact, split by explicit Q-isogeny |
| [70] | 70 | 0 | Howe15 f70 | exact; repo lane8 |
| [3,24] | 72 | 1 | 8100.e4 | **extra** |
| [6,12] | 72 | 0 | **new** (V1, V2) | HLP §3.5 containment; repo C_HLP (unlogged); **exact witnesses this session** |
| [2,6,6] | 72 | 1 | 132300.eo2 | HLP §3.7; repo contact-6 bank (7 classes) |
| [2,48] | 96 | 1 | 18900.e2 | **extra** |
| [2,2,24] | 96 | 0 | **new** (HLP §3.7 pts) | HLP §3.7 containment; **first exact witnesses this session** |
| [2,2,4,8] | 128 | 0 | **new** (HLP §3.7 pt) | HLP §3.7 containment (+10 huge repo anchors); **first exact witness this session** |

Largest known split order: **128** = [2,2,4,8] (the maximum of HLP's
Table 1); largest cyclic: **70** (Howe15).  On the geometrically simple side
the current maximum order is 96 = [2,2,2,12] (this project), so split
Jacobians still hold the order record.  (Corrected 2026-08-13: an earlier
revision of this file misstated the orders of nine multi-factor rows,
including listing [2,6,6] as order 144; orders are now computed from the
invariant factors.)

## 2. Groups not a subgroup of any product E1(Q)tors × E2(Q)tors ("isogeny extras")

**[16], [19], [25], [27], [48], [2,16], [2,48], [3,24], [4,16], [2,2,16]** —
these can only arise through the isogeny J → E1×E2 having nontrivial effect on
torsion.

**Added 2026-08-26: [11]** (isogeny-extra of conjugate-pair type; 11 is not a
factor of any Mazur-group order).  First split realization — previously
simple-only (Ogg 353.a).  Witness in `product/data/new_split_witnesses.txt`:
the (2,2)-gluing of the X1(11) quadratic-point curve E/Q(√11) (s0 = 4/5 on
r² − r(s³−3s²+4s) + s = 0; E(K)tors = [11]) with its conjugate, Weil-descended
to Q; J(Q)tors = [11] exactly, J is Q-simple (E ≁ E^σ) and geometrically
split.  With this, **76 nontrivial groups** have exact split witnesses.
Construction lane `product/code/lane_qglue.m`; the same machinery targets
[13] (X1(13) quadratic points) next.

**Sanity check vs `g2torsion.pdf`**: the slides list C27, C48, C2×C48, C3×C24,
C4×C16, C2×C2×C16 and correctly flag C19, C25 as Q-simple-only.  Two
omissions: **C16 and C2×C16** are also not subgroups of any product (a product
has 2-power cyclic part at most C8) yet occur for Q-split Jacobians (2700.d1,
4608.h5).  Everything else on the slides agrees with the DB + literature.
(Also: HLP contains no order-65 example — 63 is their largest cyclic; 65 is not
listed on the slides either, so no discrepancy, just a note.)

## 3. Q-simple but geometrically split

[19] (J1(13), unique DB example) and [25] (42336.dw1, unique) occur **only**
for Q-simple, geometrically split Jacobians — matching the slides.  Many small
groups also occur Q-simple geometrically split (CM / M2(Q) types, e.g. X1(16),
X1(18), 784.c1).

## 4. New witnesses from this session (all verified: exact torsion + split)

Full data in `product/data/new_split_witnesses.txt`; verification logs
`product/logs/verify_witnesses{,2}.log`.  All from Genus2Elliptic2 (2-torsion
gluing, Bröker–Howe–Lauter–Stevenhagen / HLP implementation in `genus2.m`).

- **[8,8]** (order 64): W1: y² = 836x⁶+88596x⁵+88597x⁴+1800118x³−4045487x²−4535664x+84285504
  (factors: [2,8]-curves 82110.bs5 × 210.e6).  W2 (even model):
  y² = 814016x⁶+1077715393x⁴+246368752516x²−3741868748800 (46410.ck6 × 210.e6).
  ~120 further examples in the logs; every pair of [2,8]-curves with matching
  2-torsion field yields [8,8] under the "cross" gluings.
- **[6,12]** (order 72): V1: y² = 132x⁶+396x⁵−6347x⁴−13354x³+75207x²+81950x+88825
  (factors [2,6]: cond 630 × 30.a6).  V2 (even):
  y² = 19299418112x⁶−178475623260x⁴+413375730009x²−15993254144
  ([12]-curves, cond 3131310 × 274170).  11 examples found.
- **[2,2,24]** (order 96): from HLP §3.7's own parameter point (y,u)=(2/9,1/3),
  i.e. (t,u)=(241/81,1/3) on E^t_{2,6} × E^u_{2,8} with condition (6):
  y² = 581449680x⁶−4134794160x⁵−574778279x⁴+10516435114x³−15881777387x²+11223443268x−1729978236
  (factors 249690 [2,8] × 825330 [2,6]).  Second instance at (t,u)=(−95,3/7).
- **[2,2,4,8]** (order 128): instance (t,u)=(1/6,4/17) of HLP condition (5) on
  E^t_{2,8} × E^u_{2,8}:
  y² = −144295356865660x⁶+289009358554092x⁵+860505249465645x⁴−1006367755763986x³−1568053362370059x²+812225828599020x+949701503960100
  (factors [2,8]-curves, cond 19291662390 × 82110).  Far smaller than the ten
  height-10^19..10^108 anchors in `paper/scripts_and_data/ten2248models_abcd.txt`.

With these, **every group in HLP's Table 1 now has an explicit exact-torsion
witness**.

## 5. The product benchmark and the remaining gaps

There are 96 distinct products of two Mazur groups; **66 of the 96** now have
exact split-Jacobian witnesses.  The remaining **30 product gaps** (occur for
E1×E2, unknown for Jacobians), with route classification:

- **(3,3)-glue could give the full product** (needs an anti-3-congruent pair
  with the stated torsions; pairs on Frengley's Z(3,2) surface):
  [56]=[7]×[8], [2,28]=[7]×[2,4], [2,40]=[8]×[10], [2,56]=[7]×[2,8],
  [10,10]=[10]×[10], [2,2,20]=[10]×[2,4], [2,2,40]=[10]×[2,8],
  [2,8,8]=[8]×[2,8], [2,2,8,8]=[2,8]×[2,8].
- **2-glue image loses index 2; needs a §3.7-style gain**: [2,60]=[10]×[12],
  [4,24]=[8]×[12], [9,9] (also needs a 2-congruent pair of 9-torsion curves —
  none exists in the LMFDB nor on X1(9)² to parameter height 20), [12,12],
  [2,2,2,12], [2,2,2,24], [2,2,6,6].
- **no 2-/3-gluing route from Mazur pairs** (every product decomposition has
  mismatched rational 2- and 3-structures): [42], [72], [84], [90], [2,36],
  [2,42], [2,72], [3,18], [3,36], [6,18], [2,2,30], [2,4,12], [2,4,24],
  [2,6,12].  (Some of these may still arise via gains on top of other gluings,
  (N,N)-gluings for N ≥ 4, or entirely different constructions.)

**[2,2,6,6] status (2026-08-14 session, `notes/claude_split_2266_sigma_session_2026_08_14.md`):**
the (2,2)-gluing route through E26×E26 is now closed on every structured
stratum.  Each of the six sigma-surfaces (root-matching class systems) has a
*universal deck section* — a Möbius map t = δ(u) with j(E26(t)) = j(E26(u))
satisfying all three square conditions identically (id: t=u; (12):
(u+15)/(u−1); (13): 6−u; (23): (5u−9)/(u−5); (123): (u−21)/(u−5); (132):
(5u−21)/(u−1)) — and every full-condition rational point found by any method
(joins to height 200, per-fiber Mordell–Weil boxes over ~4–5k fibers per
surface, 48 rational-curve sweeps to u-height 2000) lies **on** a deck
section, where the gluing degenerates (j-equal pair).  The deck-section
multiples hit the same genus-7/23/47/79/119 condition-curve wall as the
sigma=id fibration; the special Möbius curves with one condition identically
satisfied (t = u±4, (6u−9)/u, −9/(u−6) on (13); (21u−9)/(5u−9) on (123);
(5u+11)/(u−1) on (132)) are all **provably empty — locally insoluble at
p = 3** — while the ambient surfaces are everywhere-locally-soluble.  The
E26-family trace sieve (21 658 curves, all pairs, mod 5 and mod 7) found no
congruent pair to height 150, closing the product-injection alternative at
family scale.  Verdict: the locus is genuinely thin (deck-confined), not
locally obstructed; [2,2,6,6] remains open but every known structured attack
on this chart is exhausted.  (Post-review addendum, 2026-08-16/18: review
uncovered that most fiber boxes had enumerated index-2 sublattices; the
corrected coset-complete rerun over verified-saturated lattices — 4070 of
6444 rank-2+ fibers enlarged — still found ONLY deck points, and a partial
idx=1 certificate pass (532 fibers, zero violations) confirmed the lattices
everywhere it looked.  This path to a realization is not completely ruled
out — no negative claim is made — but it now seems less likely than the
other open routes above.)

## 6. Session construction summary

| lane | method | outcome |
|---|---|---|
| glue2 | 2337 field-matched LMFDB pairs, all Genus2Elliptic2 gluings funneled | 63 hits: [8,8] (from the six [2,8]-curves), [6,12] |
| x1fam | X1(N) family sweeps N∈{7,8,9,10,12} (h≤20), field-bucketed pairs | 68 hits: [8,8], [6,12] beyond LMFDB conductors |
| hlp37 | HLP §3.7 conditions (5)/(6) swept, instances glued | exact [2,2,24] ×2, [2,2,4,8] ×1 |
| siblings | Richelot webs + isogeny-class (2,2)/(3,3)-gluings around 40 anchors | 342 exact-known, nothing new (DB thorough in its box) |
| gainhunt | non-generic-factorization prefilter, N∈{6,8,10} h≤24-32 | gains land on known groups only |
| biell | y²=c_t(x²) sweeps over LMFDB seeds | nothing beyond DB |
| richelot | Richelot closures of 42 split anchors | nothing new |
| z3 | Frengley Z(3,2) sweep + (3,3)-gluing (targets [2,40], [10,10], ...) | 7306 surface points, 0 usable moduli pairs (ZNrModuli interface mismatch — parked; the (3,3)-route classification in §5 stands as future work) |
| 2266sigma (08-14) | five non-id sigma-surface 2-of-3 joins H=150 + third-condition upgrade | 62 079 surface points listed; every full match a deck pair; 0 genuine |
| 2266sigfib (08-14) | per-fiber MW over each sigma-surface (deck-section base points, RankBounds + bounded TwoDescent) | ~4–5k fibers/surface, rank 1–2 generic; survivors = deck points only |
| 2266lines (08-14) | 48 rational curves on the surfaces, symbolic reduction + u-sweep to height 2000 | all full points on deck sections; the 9 special Möbius curves provably empty (p=3 obstruction) |
| cong_sieve E26 (08-14) | all 21 658 E26 curves pairwise mod-5/mod-7 trace sieve (MINCMP 14) | 0 congruent pairs to H=150 |
| tw cleanup (08-14) | genus-3 twisted diagonals: unconditional RankBound 0 + class-order sieve | C(Q) exhausted by degenerate points — [10,10]/[12,12] twisted-diagonal impossibility now unconditional |

Funnel: 12 good odd primes p, J(F_p) invariants, sound abelian-group
compatibility test against the known set (odd part pinned exactly for
2-gluings), exact `TorsionSubgroup` only for survivors.  Files:
`product/code/split_lab.m` (shared), `lane_*.m`, logs in `product/logs/`.
