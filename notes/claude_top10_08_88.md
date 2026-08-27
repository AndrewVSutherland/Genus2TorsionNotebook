# top10 #8: (8,8), order 64 — Nicholls (4,4)-substrate validated; first (8,4)-halving box scan

Date: 2026-07-17.  Target: J(Q)_tors = Z/8 x Z/8 on a geometrically simple genus 2
Jacobian over Q.  Fresh substrate (per dossier): Nicholls' thesis (Oxford 2018),
Prop 5.9.6 — a 3-parameter family (s,t,v) with pointwise-Q-rational (4,4)-kernel.
Working directory (scripts + raw logs):
/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/t88/

## Dossier summary

Split (8,8) exists (HLP gluing).  For simple: the only in-repo route was halving
the second generator on the M_1(8,4) = [4,8] family; no local obstruction through
p = 43, but the second-halving locus H_eta is a general-type surface (conjugate
z-degree-8 components of fiber genus 23/31, min observed 21), CRT/height searches
to 2000 structurally dry, and the 2-step Richelot neighborhood of all 82 seeds
only LOSES 2-power torsion (notes/claude_tier2_88_richelot_sweep.md).  Verdict on
file: need a construction with TWO INDEPENDENT 8-chains over a bigger base.
Nicholls Prop 5.9.6 supplies a 3-parameter base with (4,4) pointwise rational —
one more degree of freedom than the dead 2-parameter chart, and both order-4
generators visible from the start.

## The substrate (exact construction, from the thesis text)

Source: scratchpad/nicholls.txt lines 5747-5990 (Example 5.9.4, Sections
5.9.2-5.9.3, Prop 5.9.6, Tables 5.1-5.3).

- Parameters s,t,v in Q.  A := s^2 - t^4 + t^2,
  u := (-s^2 A v^2 - 2 A v - 1) / (-s^2 t A v^2 + t),
  (a,b,c) := ( A/(1-t^2), A/(u^2 s^2 + 1 - t^2), t^2 ),
  d2 := A (s^2 u^2 + t^4 - 2t^2 + 1)(s^4 u^2 - s^2 t^2 u^2 + s^2 u^2 - t^6 + 3t^4 - 3t^2 + 1).
- Quintic-side curve C2: y^2 = x(x-1)(x-a)(x-b)(x-c) (all 6 Weierstrass points
  rational).  Richelot splitting (x, (x-1)(x-a), (x-b)(x-c)) gives the codomain
  C1: y^2 = L1 L2 L3 with (Example 5.9.4, common denominator a-bc dropped)
  l1 = (-a+b+c-1)x^2 + (2a-2bc)x + (abc-ab-ac+bc), l2 = -x^2+bc, l3 = x^2-a;
  f1 := l1 l2 l3 / lead.  Then C1: y^2 = d2 f1(x) has a (4,4)-kernel Sigma
  pointwise rational: (Z/4)^2 <= J1(Q)_tors.  (Twist class of d2 is what matters;
  scripts use the squarefree part.)
- Provenance of the family: c = t^2, (1-t^2)(a-t^2) = s^2, b(a-b)(a-t^2) = z^2
  (row Lambda_334 of Table 5.1 discs), conic chain in (b,z) then in (u,v)
  (5.123)-(5.128); the twist d2 equalizes the a29 square classes so both
  order-4 generators lift to J1.  Table 5.2/5.3 list the OTHER Lambda_ijk
  branches (more parametrized (4,4) families — see route 2).
- Structure: J1 --phi1--> J2 = Jac(C2) --phi2--> J3, Sigma = ker(phi2 phi1) with
  Sigma[2] = ker phi1.  Generators D1, D2 (order 4, 2D_i = T_i in Sigma[2]).

Key push-forward fact used below: if E_i in J1(Q) with 2E_i = D_i (the (8,8)
condition), then phi1(E_i) has order 4 in J2(Q) (since 4E_i = T_i in ker phi1)
and 2 phi1(E_i) = phi1(D_i), a NONZERO rational 2-torsion class of the quintic
Jacobian J2.  So (8,8) on J1 forces two independent rational 4-torsion points on
the full-2-level quintic side, above two prescribed 2-torsion classes — each a
classical repo-style x-T tangent/halving condition — plus one more a29-type
square condition per generator to lift back to J1.

## Test run (all single-threaded, nice -n 15; total approx 14 CPU-min)

### 0. Calibration (reproduce a dossier fact)

calibrate.m: Nicholls' geometrically simple order-30 curve
y^2 = x^6 - 16/3 x^5 + 70/9 x^4 + 131/27 x^2 + 16/27 x + 64/81.
Command: `magma -b calibrate.m` in the t88 directory.
Result: minimal model y^2 = 4x^6-20x^5+5x^4+70x^3+167x^2+142x+73,
TorsionSubgroup = [30]; simplicity cert (Frobenius charpoly irreducible AND
12th-power transform irreducible of degree 4) PASSES at p = 13
(charpoly x^4-4x^3+6x^2-52x+169).  Note: at p = 11 the charpoly is reducible —
the cert must search primes, one reducible prime proves nothing.  PASS.

### 1. Family implementation spot check

family.m: 12 small tuples, 10 nondegenerate.  Results (exact TorsionSubgroup
after minimization):

| (s,t,v) | torsion |
|---|---|
| (2,3,1), (2,3,-1), (5,2,1), (2,3,1/2), (1/2,1/3,1), (5,3,-2), (4,3,1), (7,2,1), (2,5,1) | [4,4] |
| (3,1/2,1) | [2,4,4] |

10/10 contain (4,4) — Prop 5.9.6 verified as implemented (twist convention:
squarefree part of d2 times monic f1).  One member picks up an extra rational
2-torsion point.

cert.m: member (2,3,1) → a = 17/2, b = 670761/69169, c = 9; minimal model
y^2 = 363482437500x^5 + 30563334341100x^4 + 45986276929439x^3
- 124217445278528x^2 + 62679319617852x - 6682982158656; torsion [4,4];
SIMPLICITY CERT PASS at p = 31 (charpoly x^4-4x^3-2x^2-124x+961).
So the family generically gives GEOMETRICALLY SIMPLE Jacobians with rational
(4,4) — already a clean simple-(4,4) production line as a side product.

### 2. Box scans with local (8,4) filter

Filter: Z/8 x Z/4 must embed in J1(F_p) (2-adic valuations of the two largest
invariant factors >= (2,3)) for the first 6 good primes >= 11 — necessary for
one order-4 generator to halve over Q.  scan.m (heights <= 3, 588 tuples) and
scan4.m (heights <= 4, 2420 tuples, survivors re-tested at 10 primes; the whole
2344-curve filter runs in 12 s).

- h <= 3: 568 valid curves; members passing first k primes, k = 1..6:
  [234, 102, 41, 18, 10, 3] — per-prime survival approx 0.42, MEMORYLESS
  (pure-chance profile).  All 3 six-prime survivors have exact torsion [4,4].
- h <= 4: 2344 valid; 21 six-prime survivors; 3 ten-prime survivors
  (1/3,3/2,1/3), (2,3/4,-4), (4,4,1/2) — ALL exactly [4,4] (last.m).

ZERO (8,4)/(4,8)/(8,8) members through height 4.

### 3. Quintic-side stage-1 diagnostic

scan_j2.m (h <= 3 box, same 568 members, on J2 = Jac(C2) which has full rational
2-torsion): local embeddability of (2,2,2,4) into J2(F_p) for k = 1..6 primes:
[359, 261, 196, 162, 129, 114] — flattens at 20%, NOT memoryless.  But j2exact.m
(exact torsion of the first 4 six-prime survivors) gives [2,2,2,2] every time.
Explanation: the (2,2,2,4)-anywhere filter is a union over the 15 Weierstrass-pair
classes; per-class mod-p 2-divisibility is approx 1/2 per prime, so
E[some class survives k primes] approx 15/2^k → 15/64 = 23% at k = 6, matching
the observed 20%.  The flattening is class-multiplicity, not a rational
stage-1 sublocus.  The both-chains condition (2,2,4,4) decays to 0 by k = 6.

## Interpretation / verdict

- Substrate VALIDATED: pointwise-rational (4,4) for free on a 3-parameter
  rational base, generically geometrically simple.  Structurally strictly
  richer than the dead M_1(8,4) chart (both generators present, +1 dimension,
  quintic side with full rational 2-torsion where the repo's x-T machinery
  applies verbatim).
- But NO free lunch: no halving hits at height <= 4, and both local filters
  behave exactly like chance models.  The halving loci are thin covers not hit
  by naive boxes — same lesson as the old chart, now established cheaply here.
- Go/no-go: GO for the symbolic first-halving locus derivation (route 1);
  NO-GO for larger naive (s,t,v) box scans — they measure nothing but the
  chance model until the locus is derived and parametrized.

## Strategy (ranked)

1. **Two-stage symbolic halving on the Lambda_334 family (primary).**
   (i) Identify the prescribed classes phi1(D_1), phi1(D_2) in J2[2](Q)
   explicitly (one Kummer/Richelot image computation on a sample member, then
   symbolically — they are Weierstrass-pair classes of {0,1,a,b,c,inf}).
   (ii) Stage-1 necessary condition: each phi1(D_i) 2-divisible in J2(Q).  On
   the quintic y^2 = h(x) := x(x-1)(x-a)(x-b)(x-c) this is the repo-standard
   tangent condition: h(x) - q_i(x) (Mx+N)^2 = (cubic)^2-type square system
   with q_i the pair quadratic — one cover layer, historically often
   parametrizable (cf. the m18_m14 first cover; Nicholls' own conic-chain
   trick (5.123)-(5.128) is the model).  Deliverable: branch equations and,
   if possible, a rational parametrization of the (8,4)-candidate subfamily.
   (iii) Lift layer: on that subfamily impose E_i rational (a29-square
   condition via Nicholls Sec. 5.9.2 one level up / Kummer duplication
   delta(xi) = xi(D_i)), then the second generator's stage-1 + lift.
   Dimension budget: 3 - 1 - 1 = 1, so if the covers stay rational the (8,8)
   locus is a curve — vs codim-2-in-2 in the old chart.
   EARLY WARNING STEP: before any search, compute mod-p function-field genera
   of the stage-1 cover fibers (heta_genus.m methodology) — if genus explodes
   like H_eta (>= 21), abandon the branch immediately.
2. **Other Lambda_ijk branches (Tables 5.2/5.3).**  The thesis parametrizes
   several further pointwise-(4,4) families (rows 311-344, 211-213 etc.; some
   need u^2 = quartic, i.e. elliptic-fibered bases).  The halving-cover geometry
   varies by branch; if Lambda_334's stage-1 cover is general type, rerun
   route 1 on the next parametrizable branch before giving up on the substrate.
3. **Fallback (in-repo): special curves on the old X_i surface.**  Per
   claude_tier2_88_richelot_sweep.md: fiber-genus-drop s0 values (genus 21 at
   p=97, s0=5), the involution swapping the conjugate deg-8 z-components,
   boundary strata.  A rational point must live on such a curve if
   Bombieri-Lang holds.  Low prior; only after routes 1-2 hit walls.

## Exact commands

```
cd /tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/t88
nice -n 15 magma -b calibrate.m   # order-30 calibration        (~40 s)
nice -n 15 magma -b family.m      # 12-tuple spot check          (~1 min)
nice -n 15 magma -b cert.m        # simplicity cert on (2,3,1)   (~1 min)
nice -n 15 magma -b scan.m        # h<=3 box + (8,4) filter      (~2 min)  -> scan_out.txt
nice -n 15 magma -b scan4.m       # h<=4 box, 10-prime filter    (~1 min)  -> scan4_out.txt
nice -n 15 magma -b last.m        # exact torsion (4,4,1/2)      (~30 s)
nice -n 15 magma -b scan_j2.m     # quintic-side stage-1 filter  (~1 min)
nice -n 15 magma -b j2exact.m     # exact check of 4 survivors   (~2 min)
```

## Files

All in the t88 scratchpad directory listed above: calibrate.m, family.m,
cert.m, scan.m + scan_out.txt, scan4.m + scan4_out.txt, last.m, scan_j2.m,
j2exact.m.  Thesis text: scratchpad/nicholls.txt (Prop 5.9.6 at lines
5885-5902, Example 5.9.4 at 5747, Tables 5.1-5.3 at 5993-6060).

## Unverified / caveats

- "Generically simple" rests on one certified member ((2,3,1), p = 31); other
  members were not certified (the (4,4)-kernel makes splitness a real
  possibility on subloci — check per hit).
- The identification of the prescribed classes phi1(D_i) was argued, not
  computed; step (i) of route 1 must pin them down before the symbolic stage-1
  derivation (the push-forward argument itself — 4E_i = T_i in ker phi1, so
  phi1(E_i) has exact order 4 — is elementary and solid).
- The (8,4) local filter checks group-structure embeddability, not divisibility
  of the specific generators D_i; it is a NECESSARY-condition filter only
  (fine for the no-hits conclusion, too coarse to certify a hit).
- Torsion computations used the squarefree-part twist of d2; same Q-isomorphism
  class as y^2 = d2 f1, but coefficient sizes (hence minimization behavior)
  differ from Nicholls' printed models.
- Box symmetries: s -> -s and t -> -t fix the curve (only s^2, u^2 enter after
  simplification; t -> -t flips u's sign only), so scans used s,t > 0; v has no
  such symmetry and ran over both signs.  v = 0 and s = 0 are degenerate
  (b = c resp. b = a), excluded structurally.
