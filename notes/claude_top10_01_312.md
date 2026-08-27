# top10 item 1: (3,12) — strategy + test run

Date: 2026-07-17.  Target: torsion (3,12), order 36 — the smallest group not yet
realized by a geometrically simple genus 2 Jacobian /Q (rank #1 on
notes/claude_top10_ranking.md).  All new work in the session scratchpad
(`scratchpad/t312/`); repo untouched except this note.
Total compute this session: ~25 CPU-min, single-threaded, nice -15.

## Dossier recap (inputs)

- Split realizations abound: HLP P^2 family; three exact split [3,12] hits on the
  M(2,12) chart at (z,r) = (-5/3,-3/5), (-5/3,-9/35), (-5/4,-32/35); the first has
  odd model Y^2 = 5668704X^5 - 22143375X^4 + 36098622X^3 - 30305259X^2
  + 12990780X - 2259900 and splits as Cremona 90c3 x 510g1
  (notes/m212_three_torsion.md).
- Proven (notes/claude_tier2_266_312_contact6.md): over the M(2,12) chart the
  residual 3-torsion cover splits over Q(z,r) into exactly two irreducible
  components S12 (degree 12, Weil-isotropic w.r.t. the built-in class 4D) and S27
  (degree 27); every rational [3,12] point lies on S12; the split hits are a thin
  set on it; the geometry of S12 (rational / elliptic-fibered / general type) was
  UNDETERMINED ("genus expected large") — flagged as the decisive open question.
- Open-chart obstruction at p=5 on M(2,12): all good residue pairs have
  J(F_5)[3] rank <= 1, so rational [3,12] candidates live over the p=5 boundary.
- Untried routes: pure M(12) + independent cubic-contact-3 (drop the extra-2);
  Bruin-Flynn-Shnidman sqrt3-RM full-level rational surface (arXiv:2102.04319)
  as a (3,3)-substrate on which to impose a 4-torsion class.
- Machinery reused: the 3-torsion cubic-contact system E2z,E1z,E0z in
  Q[z,r,M,U,V] (h3^2 - f5 = M q^3, q = x^2 + Ux + V, M = m^2) and the
  multiplication-matrix eliminant method, both from the earlier tier2 session
  (recovered from `tier2_266/b312_equations.m`, copied into `t312/`).

## Strategy: three ranked routes

### Route 1 (primary): settle the geometry of S12 and exploit its low-genus pencil

The decisive question from the previous session was whether S12 is of general
type.  THIS SESSION'S TEST ANSWERS THE FIBER-GENUS PART: the z-pencil of S12 has
fibers of GENUS 4 (see test run below) — much smaller than feared, and the
special slice z=-5/3 through two of the three known hits has the same genus as a
generic slice, i.e. the hits do NOT sit on a degenerate fiber.  A genus-4-fibered
surface is far more workable than a general-type black box:

- Each z-slice is a Faltings-finite genus-4 curve; the z=-5/3 slice already
  carries >= 6 known rational points (3 extra-3-torsion points over each of the
  two hit values r=-3/5, r=-9/35, all with M a rational square).  A genus-4 curve
  with 6+ points is rich; computing its full Mordell-Weil-quotient/Chabauty
  structure is feasible with current tools if a rank <= 3 quotient exists.
- Concrete exploitation plan: (a) interpolate the degree-12 U-eliminant factor
  F12(U; z, r) over Q(z, r) along more slices (cheap now: ~165 s per slice mod
  p); (b) hunt for SPECIAL slices — values of z (or r, or other pencils, e.g.
  z = c*r + d lines through the p=5-boundary residue classes) where the fiber
  genus drops to <= 1; the degeneration locus of a genus-4 fibration is a
  divisor, and any rational point of it with fiber genus <= 1 gives an
  infinite-sweep parametrization; (c) if no genus-drop slice exists, run the
  known genus-4 slices through the simplicity certificate at their known points
  and finish z=-5/3 by explicit Chabauty (its points = all [3,12] curves on that
  line; deciding whether all are split would convert the empirical confinement
  into a proved statement for this line).
- Simplicity screen for any new point: L_p irreducible + 12th-power transform
  irreducible of degree 4 (validated below on the split hit — it FAILS there, as
  it must).

### Route 2 (new substrate, best structural leverage): BFS sqrt3-RM surface + sqrt3-division of D1 + one 4-torsion class

(3,12) = Z/3 x Z/3 x Z/4.  The Bruin-Flynn-Shnidman full-level-sqrt3 family
(arXiv:2102.04319; cached at scratchpad/full2.txt) is the RATIONAL surface
P^2_{(a:b:c)} \ Delta with universal curve C_{a,b,c}: y^2 = f_{a,b,c}(x),
disc_x(f) = 2^12 a^3 b^9 c^14 (a-c)^12 q3^3 t3^3 q4^4 t4^4, where
G1^2 + lambda1 H1^3 = q4^2 f_{a,b,c} (their eq. (9); explicit G1, H1,
lambda1 = 4ac t3 in their Sec. 3).  CORRECTION to the dossier framing: full
sqrt3-level here means J[sqrt3] = Z/3 x mu_3, so only ONE pointwise-rational
3-torsion class D1 (Mumford (H1/lc, G1 mod H1)) is free; the mu_3 part has no
rational points.  The right way to get the second rational 3-class w (needed
since (3,12) has 3-rank 2) is NEW and sharp: sqrt3 is a rational endomorphism
and sqrt3*w is rational, so w rational and independent forces

    sqrt3 * w = +-D1   (w in J[sqrt3](Q) = <D1> is excluded),

i.e. w lies in the sqrt3-DIVISION TORSOR of D1.  As unordered pairs {w,-w}
this is a degree-9 cover T9 -> P^2_{(a:b:c)} (18 elements {w : sqrt3 w = +-D1},
mod +-), STRICTLY smaller than the M(2,12) carrier S12 (degree 12), living over
a rational base, and carrying torsor structure under J[sqrt3]: whether a fiber
has a rational point is a class in H^1(Q, J[sqrt3]) — attackable by the
explicit 3-isogeny descent machinery BFS themselves develop (their
A = J/<D1>, B = J/<D2>, Cor. 3.7: A_{-3} ~ B at the involuted parameter
(q3 : 3b(b+a)-q3 : q4)).  Equivalently, w exists iff D1 is hit by A(Q) under
the complementary isogeny psi: A -> J with kernel J[sqrt3]/<D1> = mu_3.
Plan: (i) build T9 with the same cubic-contact eliminant machinery
(h^2 - f_{a,b,c} = M q^3 over Q(a,b,c); the residual 39 pairs must split
FURTHER here than the generic [12,27] because RM constrains the monodromy —
expect a component of degree <= 9 defined by sqrt3 w = +-D1); verify the
predicted degree-9 component and compute its slice genera exactly as done for
S12 below; (ii) rational points on T9 give (3,3)-RM curves in bulk; on those,
filter for a rational linear factor of the sextic and run the repo's validated
x-T 4-divisibility criteria (C presieve ports in ~1 h) to force the Z/4;
(iii) certificate-screen.  The local facts below (bad reduction at 2 and 5
forced) are cheap pre-filters.  Bonus: every T9 point is also a fresh
simple-(3,3)-with-RM candidate independent of the (3,12) goal.

### Route 3: pure M(12) + independent contact-3 over A^2_{a,r}

Mirror of the b312 system without the z-level-2 structure: on
y^2 + (x-r)(T+1)y = a x^2 T(T+1), T = ax^2 - x + r (marked order-12 class),
impose an independent 3-torsion class via cubic contact on the completed square
W = h^2 + 4a x^2 T(T+1) (even model — needs the weight-2 contact system or an
odd model after locating a rational root of W; the T+1 root trick used on
M(2,12) is not available without z, so use the generic even-model contact
h^2 - W = M q^3 with deg h = 3).  The residual cover should again split by
Weil-pairing into isotropic/anisotropic parts; the isotropic surface S12' over
A^2_{a,r} pulls back to S12 under the double cover z^2 = 1 - 4a(r+1).  S12'
could a priori have better geometry than its double cover.  LOCAL VERDICT FROM
THIS SESSION: the p=5 obstruction persists on the full M(12) chart (0 of 9
good-reduction pairs admit 3-rank >= 2 — see test run), so this route inherits
exactly the same p=5-boundary confinement as M(2,12); p=7,11,13 are alive with
healthy densities (6/24, 12/78, 39/116 good pairs).  Rank this third: same
thin-set risk as Route 1 with a new eliminant to build; do it only if Route 1
finds S12' -relevant structure (e.g. if S12 turns out to be the pullback of a
smaller surface, the descent to S12' is the natural place to see it).

## Test run (executed this session)

Scripts in `scratchpad/t312/`: `val_split_hit.m`, `m12_local_feas.m`,
`f5_universal.m`, `f2_universal.m`, `slice_genus.m` (+ `b312_equations.m`
copied from the tier2_266 scratch area).  All runs `nice -n 15 magma -b`,
single-threaded.

### 0. Validation (calibration against the dossier)

`val_split_hit.m`: the dossier's odd model at (z,r,a)=(-5/3,-3/5,-10/9):

    torsion invariants: [ 3, 12 ]
    p=7,11,13,19,23: chi_p REDUCIBLE (factor degrees [2,2] or [2]),
    simplicity certificate FAILS at every tested prime (p=17 bad).

Reproduces the known facts exactly (torsion [3,12]; split as 90c x 510g, so
every L_p factors).  ~30 s.

### 1. M(12)+extra-3 local feasibility (`m12_local_feas.m`)

For p in {5,7,11,13}, all (a,r) in F_p^2 with good genus-2 reduction of W;
count pairs with J(F_p)[3] rank >= 2 and with J(F_p) >= Z/3 x Z/12
(necessary conditions for reduction of a rational [3,12] point):

    p  5: total  25  good   9  rank3ge2  0  contains_Z3xZ12  0
    p  7: total  49  good  24  rank3ge2  6  contains_Z3xZ12  6
    p 11: total 121  good  78  rank3ge2 12  contains_Z3xZ12 12
    p 13: total 169  good 116  rank3ge2 39  contains_Z3xZ12 39

The p=5 obstruction PERSISTS after dropping the extra-2 condition; p >= 7 all
alive.  ~1 min.

### 2. Universal local facts at p=2 and p=5 (new small theorems, chart-free)

`f2_universal.m`: ALL 768 genus-2 curves /F_2 (y^2 + h y = f, deg h <= 3,
deg f <= 6): NONE has (3,3) contained in J(F_2).

> Consequence: EVERY genus-2 Jacobian /Q with torsion containing (3,3) — in
> particular every (3,12) — has bad reduction at 2.

`f5_universal.m`: all 27,500 squarefree models y^2 = f over F_5 (deg 5 monic
and deg 6 with lc in {1, nonsquare} — this covers every isomorphism class,
since odd-degree leading coefficients normalize to 1 and even-degree ones to
their square class):

    total models 27500   with 36 | #J(F_5): 367   with J(F_5) >= Z/3 x Z/12: 0

> Consequence: EVERY genus-2 Jacobian /Q with torsion containing Z/3 x Z/12
> also has bad reduction at 5.  The p=5 open-chart obstruction found on
> M(2,12) (and reconfirmed on M(12) above) is UNIVERSAL, not a chart defect.

Combined lemma (new, chart-free, same genre as the repo's forced-bad-reduction
theorems for (2,2,4,8) and (2,2,2,12)): any genus-2 Jacobian /Q with a (3,12)
subgroup has conductor divisible by 10 (bad reduction at 2 AND 5); the (3,3)
part alone already forces the 2.

Cross-checks: the split hits' elliptic factors have conductors 90 = 2*3^2*5 and
510 = 2*3*5*17 (both divisible by 2 and 5); all LMFDB (3,3)-containing curves
have even conductor (e.g. 26244 = 2^2*3^8, 196 = 2^2*7^2).

### 3. S12 slice geometry (`slice_genus.m`) — the decisive question, TEST scale

Method: fix z = z0 (resp. r = r0), work over k = F_10007(t); the contact system
is a 0-dim degree-40 algebra over k; multiplication-matrix minimal polynomial of
U factors over k; the degree-12 factor generates the function field of the S12
slice; Magma `FunctionField` + `Genus`.  ~165 s per slice (GB 5 s, minpoly 73 s,
genus 1-86 s).

    slice z = -5/3 (hit line):    factors [1, 12, 27];
        S12 slice: degree 12, GENUS 4      S27 slice: degree 27, genus 13
    slice z = 2/7 (generic):      factors [1, 12, 27];
        S12 slice: degree 12, GENUS 4      S27 slice: degree 27, genus 13
    slice r = 1/2 (cross pencil): factors [1, 12, 27];
        S12 slice: degree 12, GENUS 16
    cross-check at p = 31013 (z = -5/3): factors [1, 12, 27];
        S12 slice: degree 12, GENUS 4  (matches p = 10007)

So S12 carries a GENUS-4 FIBRATION over the z-line (pi_z: S12 -> P^1_z), the
same genus at the hit line as at a generic z (the hits do NOT sit on a
degenerate fiber), while the cross pencil pi_r has genus-16 fibers.  The old
working assumption "genus expected large" is refuted for the z-pencil: S12 is
far more structured than feared.  The z = -5/3 fiber is a genus-4 curve over Q
carrying >= 6 known rational points (the three extra (U,V,M) triples over each
of r = -3/5, r = -9/35 listed in claude_tier2_266_312_contact6.md).

## Results summary and verdict

New facts established this session (all at TEST scale, ~25 CPU-min total):

1. **S12 carries a genus-4 fibration** pi_z: S12 -> P^1_z (fiber genus 4 at
   both a generic z and the hit line z = -5/3; cross pencil pi_r has genus-16
   fibers; the S27 sibling has genus-13 z-fibers).  The previous session's
   "genus expected large" is refuted; the decisive S12-geometry question is
   now reduced to concrete, cheap computations on a genus-4 pencil.  The
   z = -5/3 fiber is a genus-4 curve /Q with >= 6 known rational points.
2. **Universal forced bad reduction at 2 and 5** for any genus-2 Jacobian /Q
   with (3,12)-torsion (exhaustive F_2 and F_5 enumerations; the (3,3)
   subgroup alone forces p=2).  Explains the p=5 chart obstructions as
   universal; gives a rigorous conductor pre-filter (10 | N) for all searches;
   consistent with the split hits (conductors 90, 510) and every
   (3,3)-containing curve in LMFDB (all have even conductor).
3. **The p=5 obstruction is NOT an artifact of the extra-2 level structure**:
   on the full M(12) chart it persists (0/9 good pairs at p=5), while p = 7,
   11, 13 admit healthy densities of (a,r) with J(F_p) >= Z/3 x Z/12
   (6/24, 12/78, 39/116) — no new local obstruction for the pure route.
4. Route 2 sharpened: on the BFS sqrt3-RM surface the second rational 3-class
   is forced into the sqrt3-division torsor of D1 — a degree-9 carrier over a
   rational P^2, with torsor structure attackable by the 3-isogeny descent
   machinery already in BFS.  (Paper analysis only; no computation yet.)

VERDICT: realizability still open but the outlook is BETTER than at ranking
time.  No new obstruction appeared; the two local obstructions that exist are
now understood as universal (hence harmless to route choice — they just force
conductors divisible by 10); and the carrier surface S12, previously a
black box, is a genus-4-fibered surface — structured enough for a real
geometry determination (rational / elliptic / general type) at production
scale.  The split hits remain the only known rational points of S12, so the
thin-set risk of Route 1 is unchanged; Route 2 (BFS + sqrt3-division) is the
best hedge because its base is rational and its carrier degree (9) is smaller.

## Next steps (production scale, in order)

1. [~1 CPU-h, mod p] Interpolate F12(U; z, r) (the degree-12 U-eliminant
   factor) over F_p in BOTH variables via pure per-point eliminant
   computations (~5 ms/point, ~3-4k points) + 2D rational reconstruction —
   NOT via per-slice function-field minpolys (177 s each).  Then: disc_U along
   the z-pencil, locate ALL degenerate fibers of pi_z, compute their genera;
   any fiber with genus <= 1 and a rational point gives a sweepable family.
   Check also the boundary values z = +-1, r = -1 excluded by the chart and
   the p=5-boundary residue directions (Rinf+Z0 class).
2. [~30 min, Q-level] Reconstruct the z = -5/3 fiber as an explicit genus-4
   curve /Q (interpolate the degree-12 factor over Q(r) from ~100 per-point
   Q-eliminants at 0.5-3 s each; b312_qlevel.m already computes single points).
   Feed it the 6 known points; compute automorphisms/quotient maps (bielliptic
   or low-genus quotients would enable rank bounds + Chabauty).  Outcome
   either: a proof that all its rational points are split (first PROVED
   negative statement for a whole 1-parameter slice), or new points = new
   [3,12] curves to certificate-screen.
3. [~1-2 h build + short runs] Route 2 build-out: derive the extra-3 contact
   system over the BFS P^2 (same eliminant machinery, even-model contact
   h^2 - f_{a,b,c} = M q^3 with deg h = 3), verify the predicted degree-<=9
   sqrt3-division component, slice-genus it as done here for S12.  If its
   slices have genus <= 2, this becomes the primary route.
4. [cheap, opportunistic] Extend the M(2,12) residue-table search
   (code/m212_extra3_residue_search.m) from height 30 to height ~120 with the
   simplicity certificate added as prescribed in notes/m212_three_torsion.md,
   restricted by the new conductor filter (only (z,r) whose curve has bad
   reduction at 2 and 5 — most do, so the real gain is the certificate screen).
5. [writeup] The p=2/p=5 universal lemma is paper-ready (statement + two
   exhaustive finite enumerations); also record that (3,3), (3,6), (3,9),
   (6,6), (3,12), (6,12) all inherit the p=2 part.

## Caveats

- All slice genera computed mod p = 10007 (z = -5/3 additionally cross-checked
  at p = 31013, same answer — table in section 3); char-0 genus could in
  principle differ at bad interpolation primes, but two independent slices and
  two primes agreeing makes that negligible.
- The degree-27 r-slice genus was not finished (killed as irrelevant — S27
  cannot carry rational (3,12) points).
- The sqrt3-division degree-9 claim for Route 2 is a structural argument
  (rational endomorphism sqrt3 maps a second rational 3-point to +-D1), not
  yet verified by an explicit eliminant factorization.
- The F_5/F_2 enumerations cover all isomorphism classes via standard model
  normalizations; they were run once, not independently re-derived.

## Files (this session's scratchpad, `scratchpad/t312/`)

- `val_split_hit.m` (validation), `m12_local_feas.m` (M(12) local diagnostic),
  `f2_universal.m` / `f5_universal.m` (+ .log) (universal local lemmas),
  `slice_genus.m` (parametrized slice-genus engine; run logs `slice_z53.log`,
  `slice_z27.log`, `slice_r12.log`, `slice_z53_p2.log`),
  `b312_equations.m` (copied from the tier2_266 scratch area — the 3-torsion
  contact system over Q[z,r,M,U,V]).
