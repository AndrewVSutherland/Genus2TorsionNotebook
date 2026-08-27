# Lane A: T5 (u,g)-surface production sweep — SECOND (2,2,2,12) CURVE FOUND

Date: 2026-07-18.  Sibling lane A of the (2,2,2,12) campaign (parent note:
notes/claude_prod_02_22212.md).  Scripts: scratchpad `sib_A_t5surf/`
(surf_id.gp, t5surf.c, sibverify.m, fingerprint.gp, reduce2b.m, logs).
Banked: data/claude_sib_t5surf_curve2.txt (the NEW curve, full certificates),
data/claude_sib_t5_nearmisses.txt (2-of-3 dataset).  No repo file modified.

## 0. HEADLINE

**A SECOND geometrically simple genus-2 Jacobian /Q with torsion exactly
(Z/2)^3 x Z/12 was found and fully verified** (torsion [2,2,2,12] exact;
simplicity certificates p=71,103,127,137 chi & chi^12 irreducible; G2INV
distinct from curve #1).  Smallest model (s,m,n) = (2208,-8303,-7200),
minimal model y^2 + (x^2+x)y = 36x^6 + 36750x^5 - 462983772x^4
- 301623595823x^3 + 1518598238654317x^2 + 397058962729817115x
- 1282993930035013443975.  See data/claude_sib_t5surf_curve2.txt.

Additionally: **8 new representations of curve #1** and **8 of curve #2**
(each curve now has 12 known pencil representations), a THIRD point on member
-49/240 at u = 3637/7105 (beyond the old u-height-4000 scan), and two new
structural laws (Sections 2-3).

## 1. The chart (validated derivation, surf_id.gp)

All T5 hits lie on S: z^2 = (q rho'-1)(q rho'-2u+1) (their product must be a
square).  Pythagorean parameterization z^2 + (u-1)^2 = (q rho'-u)^2 with
rational g:  **rho' = (2gu + (u-1)(g^2+1)) / (2g q(u))**, q = 4u^2-6u+3;
g <-> 1/g gives the same point.  On this chart (proven by exact issquare in
Q(u,g), validated at 60 random points and by 13152-point C selftest):

* W1 = (u-1)  q [(u-1)g^2 + 2ug + (u-1)]        == V'1 mod squares
* W2 =        q [g^2 + (6-8u)g + 1]             == V'2 mod squares
* W4 = (2u-1) q [(2u-1)g^2 + (6-4u)g + (2u-1)]  == V'4 mod squares
* W3 == W1 W2 W4 mod squares (V'3 DEPENDENT on S, signs included)

so a hit = THREE positive-square conditions, each quadratic in g (all three
palindromic in g).  Integer-cleared (u=p/q0, g=a/b, QQ=4p^2-6pq0+3q0^2):
Y1 = (p-q0) QQ ((p-q0)a^2+2pab+(p-q0)b^2);  Y2 = QQ q0 (q0(a^2+b^2)+(6q0-8p)ab);
Y4 = (2p-q0) QQ ((2p-q0)a^2+(6q0-4p)ab+(2p-q0)b^2)  — sizes ~200 N^4 G^2,
i128-safe to N,G ~ 10^5.  The known hit = (u,g) = (-97/48, 725/288).
The 6 known T5 near-misses lie OFF S (z^2 classes 10,6,3,3,6,-1) — on S,
near-miss <=> hit, so every 3-of-3 surface point is automatically a full hit.

## 2. g-RIGIDITY (the key empirical discovery)

Every verified hit representation (24 points known by end of session) lies on
one of FOUR g-fibers: **g = +-725/288 (all curve #1) and g = +-459/23 (all
curve #2)** — dedupe by F_p-fingerprint (fingerprint.gp, 8 primes) confirmed
by exact Magma G2INV on representatives (candG=curve1, candO=curve2 hashes).
Deep scans --gfiber to u-height 30000: SIX points per fiber, ALL the same
curve per +-fiber pair.  The fibration is NOT isotrivial in moduli
(isotriv.m: generic u on the fiber gives different G2) — the hit points of a
special fiber are nevertheless permuted representations of ONE curve.

**Consequence: one special |g| = one (2,2,2,12) curve. New curves <=> new
special g values.**  The pencil involution acts as (u,g) -> (u',-g).

## 3. Empirical laws (24-point dataset)

* rho'-numerator = +-(perfect square) in EVERY representation: 49, 289, 169,
  34969=187^2, 5929=77^2, 9025=95^2, 43681=209^2, 7744=88^2, 33124=182^2,
  4225=65^2, 12544=112^2.  The old "+-odd^2" refines to any square once
  even-square members appeared (88, 182, 112).  The old "denominator 240"
  law is DEAD (rd in {240, 3519, 5760, 10557, 34800}).
* members come in same-rd sets; each special fiber carries 3 members x 2
  points (u-heights <= 30000).
* curve #1 fiber: |g|=725/288 (725=5^2*29, 288=2^5 3^2, 725^2-288^2=437*1013,
  437=19*23); curve #2 fiber: |g|=459/23 (459=3^3*17; 459-23=4*109,
  459+23=2*241; note 109, 241 divide disc(curve2)).

## 4. Searches run (final tallies; ALL hits dedupe to curves 1-2: 12+12)

* Calibration N=150 G=800 (11 s, 6 threads): re-found curve1 reps 1-2 AND
  discovered curve2 (reps 1-2) + curve1 reps 3-4.  3077 NEAR2.
* Production P1 N=300 G=3000 (5 threads, 722 s): 10 hits (all on the 4
  special fibers), 7662 NEAR2, 6.6e5 1-cond passes.
* Production P2 N=150 G=8000 (5 threads, 1394 s): 5 hits (all known),
  5687 NEAR2.  **NO new special fiber to g-height 8000 at u-height <= 150,
  nor to g-height 3000 at u-height <= 300.**
* ULTRA fibers +-725/288, +-459/23 to u-height 100000: **exactly 6 points
  per fiber, no more** (near2 0/0/8/8) — each special fiber's point set is
  complete at 6 to height 1e5; 12 representations per curve total.
* Member scans (t5sweep --member) -9025/3519 to u-height 20000 and
  43681/10557 to 18000: second points found (curve2 reps 3-4).
* 183 "warm" candidate fibers (NEAR2 multiplicity >= 5 in P1) scanned to
  u-height 4000: 0 hits, 1511 NEAR2 — warm never converts, special fibers
  have FEW near2 (they convert to hits): near2-density is NOT the signal,
  E2-rank is (Section 4b).
* Fingerprint dedupe (fpfile.gp / fpclassify.txt): all 24 hit points from
  every sweep = 12 CURVE1 + 12 CURVE2, ZERO new classes.
* **T3 deep-w sweep (t3sweep 150 2500)**: 0 hits — extends the old T3
  w-locus kill from w-height 150 to 2500 (T1-passes 172338, 142 s).  If a
  T3 special fiber analog existed at |w| <= 2500 with a point at u-height
  <= 150 (both T5 curves' first points had u-height <= 121), it would have
  been seen.

## 4b. Rank criterion for special fibers (specialg.m, rankscan.m)

For fixed g the W2 condition is the elliptic curve E2(g): y^2 =
q(u)((g^2+6g+1)-8gu), and E2(-g) ~ E2(g) via u -> 3/2-u (q(3/2-u)=q(u))
— this explains the +-g fiber pairing exactly.  RankBounds:

* SPECIAL fibers: g=725/288 rank **6**; g=+-459/23 rank **5** (tors Z/2).
* Warm controls (many NEAR2, 0 hits): 32/9, -75/11, 66/5, -32/9 rank 4;
  116/9 rank 5 (rank 5 alone is not sufficient — E1/E4 must cooperate).
* Random g: rank 1-2.  Exhaustive |a|<=40, b<=16 scan: NO fiber with rank
  >= 5 exists at small height (28 fibers reach rank 4, list in rankscan.log).

So high E2-rank (>=5) is a strong necessary-condition filter for special
fibers, and the next special g has height beyond {40/16}-box — consistent
with the observed 459/23 and 725/288.

Extended scan (rankscan2.m, b=17..30, |a|<=80, completed; stalled at b=31
on a hard class group — use SetClassGroupBounds("GRH") next time): the ONLY
rank>=5 fiber is g=+-67/22 (rank [5,5]).  Deep scan --gfiber +-67/22 to
u-height 30000: **0 hits** (3 near2 each sign).  Together with 116/9 (rank
5, hitless to 4000): E2-rank >= 5 is NOT sufficient — the E1/E4 covers must
cooperate.  Next refinement: compute rank bounds of all three covers; a
triple-high-rank criterion likely pins special fibers exactly.

## 4c. Fiber orbit structure (observation)

The 6 points of fiber +725/288 pair with the 6 of -725/288 by EQUAL
rho'-denominator: members {-49/240 ~ 289/240}, {34969/34800 ~ -169/34800},
{-169/5760 ~ 5929/5760}.  Under this pairing the u-numerators form two
3-CYCLES: {97, 133, 169} (arithmetic progression, difference 36!) and
{3637, 6767, 14041}.  Same shape on the 459/23 fibers: {23, 265, 271, 317,
553} with members {-9025/3519, -4225/3519, 7744/3519, 12544/3519,
-33124/10557, 43681/10557}.  Suggests a rank-1-orbit/translation structure
on each special fiber worth deriving exactly.

## 5. Verification protocol used (jackpot)

magma -b cls:=5 pp:=P qq:=Q rn:=RN rd:=RD sibverify.m — exact integral model,
TorsionSubgroup = [2,2,2,12], >=4 simplicity certificates (chi irreducible
deg 4 AND chi^12 irreducible: p=71,103,127,137 for curve2; 37,73,113,139 for
curve1 reps), G2Invariants hash dedupe (curve1 d2e1e8c1..., curve2 abe99337...).

## 5b. Resume state (all processes stopped at session end)

* Scratchpad `sib_A_t5surf/`: t5surf.c (+binary t5surf/t5surf2 with --gfiber,
  --ufiber, --test, --selftest, q0-chunking), surf_id.gp, fpfile.gp/fpin.txt,
  sibverify.m, reduce2b.m, specialg.m, rankscan(2).m, all logs/outputs
  (p1_300_3000.*, p2_150_8000.*, gfiber_ultra.*, gcand_scan.txt,
  t3_deepw_150_2500.*, ver_cand[A-G,O].log).
* Resume box: `./t5surf 150 16000 6 [q0lo q0hi]` (next g-height doubling;
  ~90 min at 6 threads) or N=300 G=8000 (~2x P2).
* Resume rank scan: rankscan2.m pattern with SetClassGroupBounds("GRH"),
  b=31..40 + the uncovered b<=16, 40<|a|<=80 strip; probe E1/E4 ranks at
  725/288, 459/23, 67/22, 116/9 to sharpen the criterion.
* Banked: data/claude_sib_t5surf_curve2.txt, data/claude_sib_t5_nearmisses.txt
  (14903 lines), data/claude_sib_t5surf_hitclassify.txt.

## 6. Next steps (ranked)

1. **Hunt new special g fibers**: extend the (u,g) box in g-height (P2/P3
   chunks; q0-chunking supported in t5surf).  Each new special |g| is a new
   curve.  The N=150 u-height suffices empirically (all first-discovery
   points had u-height <= 121).
2. **Characterize special g**: for fixed g the three conditions are elliptic
   quartics in u; compute the rank of the E2(g) cover (Magma) at the two
   special g vs random/warm g — a rank criterion would turn the hunt into
   root-finding.  Warm fibers (many NEAR2, 0 hits: +-32/9, 116/9, -75/11,
   66/5, ...) are the natural controls.
3. Deep-extend --gfiber on the 4 special fibers (u-height 10^5, i128-safe)
   to grow the representation sets and find the translation structure.
4. The two curves' 12-member tables (rd sets {240,5760,34800} and
   {3519,10557}) suggest a finite orbit structure on members — derive the
   member-level involution formulas exactly (gp, 30 min).
5. T3-side analog: the T3 w-locus is ALSO a surface with a Pythagorean-type
   chart; rerun this exact machinery there (the identity is already proven).

## Addendum (Codex review, PR #4): checked-in scripts
The exact-commands sections above reference the session scratchpad. The key scripts are
now checked in under code/claude_sib_lanes/A/ (same filenames); sweep binaries rebuild
with the gcc lines given in the commands. Scratchpad paths remain valid only on the
discovery machine.
