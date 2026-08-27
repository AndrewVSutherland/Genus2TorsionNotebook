# Session log 2026-08-13: executing the [2,2,6,6] / [12,12] plan

Execution session for `notes/claude_split_2266_1212_plan_2026_08_13.md`.
Both targets have order 144 (would beat the current split-torsion order
record 128 = [2,2,4,8]).  This note records what was PROVED, what was
RULED OUT, and what is still running/open.  All scripts in `product/code/`,
logs in `product/logs/`.

## 1. Machinery built and validated

### 1.1 Class-(a) Prop-12 slot machinery (`lane_2266_derive.m`, log `derive2266.log`)

Symbolic iota tables (2-descent images of 2-torsion) for the HLP universal
curves E26 = E^t_{2,6} and E28 = E^u_{2,8} over Q(t), mod squares.  E26 pair
classes (ordered i<j; transposes flip sign):

    a12 = (t+3)(t-5),   a13 = 2(t-3),   a23 = -(t-1)(t-9)

The full two-gain system (rational halves of TWO independent graph classes
(T_i, T'_sigma(i)), which forces the third, i.e. FULL rational J[2]) for a
gluing with root-matching sigma in S3 collapses to exactly THREE
unordered-pair class-matching conditions

    a_{ik}(t) * b_{sigma(i)sigma(k)}(u) = square,  {i,k} in {12,13,23}.

**Anchors, all PASS**: (1) iota(T1) = HLP's printed value; (2) E28 has a
trivial iota row (T1 = 2*(rational 4-torsion)); (3) sigma=id T1-subsystem =
HLP system (6) verbatim; (4) E28xE28 sigma=id collapses to the single HLP
condition (5) (two conditions auto because of the trivial row);
(5) the [2,2,24]-witness (t,u)=(241/81,1/3) satisfies exactly a full
T1-subsystem, at sigma=[2,3,1].

### 1.2 Class-(b) halving (delta-match) machinery (`lane_1212.m`, log `lane1212_H120.log`)

For a class-(b) gluing of two X1(N)-curves (N in {8,10,12}; one rational
2-torsion point T_rat = (N/2)g1; conjugate pair over K = Q(sqrt(d))), every
order-144 (resp. 100) enlargement of the [6,12] (resp. [5,10]) image comes
from a rational half of im(3g1,3g2)-type classes, and the descent condition
collapses (delta(kg) = delta(g) for odd k; delta(T_rat) trivial for 4|N) to

    beta_t == beta_u  or  conj(beta_u)   in K*/(K*)^2,
    beta := x(g1) - e_conj  on a y^2 = cubic model,

with the rational-slot condition alpha_t == alpha_u (alpha := x(g1)-e_rat)
implied (norm checksum) and used as a hash-join prefilter.  A [2,6,12]-type
2-rank-3 gain is IMPOSSIBLE for class-(b) gluings (A[2](Q) has rank <= 2:
coset analysis; only gamma in {0, (T_rat,T'_rat)} can carry rational points).

**Validation, all PASS**:
- checksums (delta(T_rat) trivial for N=8,12; delta(T_rat) == delta(g1) for
  N=10) on ~35k pooled curves — zero violations;
- explicit negative controls: the four exact-[6,12] x1fam pairs
  (2/17,-3/7), (1/6,-1/6), (1/7,-2), (-3,-3/5) fail both variants;
- in-sweep negative funnels: unmatched d-matched class-(b) [8]x[8] pairs
  funnel to exactly [4,8] (the no-gain image), e.g. (1/8,1/24), (1/8,-1/112);
- **in-sweep positive controls: all 27 diagonal N=8 survivors funnel to
  exactly [8,8]** (see 2.1).

A pitfall found on the way: the x1fam "[8,8] hit pairs" (t=1/7 etc on X1(8))
are CLASS-(a) specializations ([2,8]-torsion, 3 rational 2-torsion points) —
not usable as class-(b) controls.  (Mazur: [2,8] allowed, so X1(8) has a
class-(a) locus; [2,10]/[2,12] don't exist, consistent with the rank-0
computations below.)

## 2. New results

### 2.1 Infinite family of [8,8] split Jacobians (self-gluing diagonal)

A class-(b) curve can be glued to ITSELF along the conjugation psi
(T_c <-> T_c-bar; valid because Aut(E) = +-1 acts trivially on E[2]; the
returned model is Genus2Elliptic2(E,E), whose nondegeneracy test kills
sigma=id automatically).  The halving condition becomes beta == conj(beta),
i.e. N(beta) = alpha a square in K:  **alpha == 1 or d mod rational
squares**.  Symbolically (`lane_1212_sym.m`, log `lane1212_sym.log`):

- N=8:  alpha == t(1-t) mod squares  ==> square iff t = a^2/(a^2+b^2):
  a CONIC — an infinite family.  All 27 such t of height <= 120 funnel to
  **exact [8,8]** (log `lane1212_H120.log`, tags `l1212diag|N=8`).
  [8,8] was previously known from isolated 2026-08-12 witnesses; this gives
  a clean 1-parameter supply t = a^2/(a^2+b^2) on the Kubert X1(8) chart.

### 2.2 Impossibility: [2,2,2,24] from [2,6] x [2,8] (2,2)-gluings

[2,2,2,24] requires full rational J[2] (2-rank 4), hence ALL THREE pair
conditions; E28's trivial iota row makes two of them t-only for every sigma:

    {2(t-3) or -2(t-3) square} AND {(t+3)(t-5) or -(t+3)(t-5)... }

mechanically (six sigmas) exactly two inequivalent t-side systems arise, and
both parametrize to rank-0 genus-1 quartics ((m^2-4)(m^2+12) resp.
-3m^4+2m^2+1; both with Jacobian [0,-1,0,1,0] = 24a, MW = Z/4).  All their
rational points map to degenerate t in {3,-3,1,5,9}.  **Hence no
(2,2)-gluing of a [2,6]-curve with a [2,8]-curve has torsion containing
[2,2,2,24].**  (`lane_2266_derive.m` + scratch rank computation.)

### 2.3 Impossibility: [2,2,6,6] from SELF-gluings of E26

On the diagonal t=u the only candidate psi are the two 3-cycles (sigma=id is
the degenerate product; transpositions carry a -1 obstruction), giving the
two-condition system C1: -2(t-3)(t+3)(t-5) = sq, C2: -(t+3)(t-5)(t-1)(t-9)
= sq (C3 = C1*C2).  All three genus-1 subcovers have Jacobian [0,1,0,-4,-4]
with MW = (Z/2)^2 — **rank 0** — and their finitely many points are all
degenerate.  So no E26-self-gluing works (`lane_2266_self.m`, log
`lane2266_self.log`; also confirms an H=200 direct scan finds nothing).

### 2.4 Diagonal ([N]-self-gluing) for N=10,12: untwisted case impossible

- N=12: alpha == t(1-t)(3t^2-3t+1) mod squares; the condition curve
  w^2 = t(1-t)(3t^2-3t+1) has Jacobian 24a ([0,-1,0,1,0]), MW = Z/4,
  rank 0: only degenerate t in {0,1/2,1}.  **No [12,12] from the
  conjugation self-gluing with alpha a rational square.**
- N=10: w^2 = t(t-1)(t^2-3t+1): Jacobian 20a ([0,1,0,-1,0]), MW = Z/6,
  rank 0: same conclusion for [10,10].
- The TWISTED variants (alpha*Delta = square) are genus-3 curves with
  **RankBound(J) = 0** (lane_misc2.m; 2-descent, GRH class-group bounds) and
  all rational points to height 2*10^4 degenerate — so the conjugation
  self-gluing diagonal is CLOSED for both N=10 and N=12 (all four condition
  curves have rank 0; only the routine rank-0 torsion-point enumeration
  remains to make it a fully formal proof).
- sigma-surface census (lane_misc2.m): each of the five non-identity
  sigma-surfaces for [2,2,6,6] has 231-388 two-of-three-condition points at
  H=40 — fibration attacks on those surfaces are possible future work
  (the sigma=id fibration itself was dry).
- Cute byproduct: the class-(a) loci of the X1(12)/X1(10) Kubert charts are
  w^2 = 12t^4-24t^3+20t^2-8t+1 (Jacobian 24a) resp. w^2 = 8t^3-8t^2+1
  (Jacobian 20a), both rank 0 with only degenerate points — a chart-level
  re-proof that [2,12] and [2,10] torsion don't occur over Q (Mazur).

### 2.5 One-sided beta-square carriers: empty to H=150

For MIXED gluings ([10]x[12] -> [2,60] order 120; [8]x[12]-class-(b) ->
[4,24]-type order 96) the halving condition is one-sided: beta itself a
K-square on one factor.  Scan (`lane_goodbeta.m`): **0 carriers among ~27k
class-(b) curves per family (N=8,10,12) up to height 150.**

## 3. Sweeps (negative so far)

- `lane_2266.m` (H=60): the 3-condition hash-join found 1561 matched (t,u)
  pairs — ALL j-equal (family deck transformations; a twist analysis shows
  matched j-equal pairs are necessarily isomorphic, and their gluings are
  the self-gluings ruled out in 2.3).  Zero genuine cross pairs.
- `lane_1212.m` (H=120): cross-pair beta-match survivors: 0 for N=8 (6298
  d-matched pairs), 0 for N=10 (273), 0 for N=12 (38).

## 4. Running / next

- `lane_2266_fiber.m`: MW-enumeration of the HLP (6)-surface fibers over
  each small u, testing the second-gain condition C3 (this reaches t of
  large height on the C1&C2-locus).
- `lane_1212_fiber.m`: same trick for Target B — fibers w^2 = D(t0)D(u) of
  the Delta-matching surface, MW-enumerated, beta-match tested.
- `lane_congtors.m`: NEW ROUTE — (N,N)-gluings along Frengley's LMFDB
  (N,r)-congruent pair lists with an anti-isometric scaling available
  ((5,1),(7,3),(8,7),(10,1),(11,2),(12,11),(13,1),(14,3),(16,3)): for
  gcd(N, #T1*#T2)=1 the product T1xT2 INJECTS into the glued Jacobian's
  rational torsion — no descent conditions at all.  Any pair with
  T1 x T2 not in KNOWN is a candidate new group; a [2,8]x[2,8] pair would
  give order 256.  Construction of the genus-2 model for a treasure pair
  would be analytic (period matrix of ExF extended by the graph lattice;
  Rosenhain -> Igusa -> rationalize -> twist-match by L-factors) — to build.

## 4.5 Afternoon developments

- **Section-arithmetic fiber attacks** (`lane_2266_fiber2.m`, `lane_1212_fiber2.m`):
  using only small combinations of the KNOWN sections (no descent) is fast
  (~5 s) but SHALLOW — section multiples' heights explode, so effectively
  |n| <= 3 gets tested.  A: 793 fibers, 2380 points on the C1&C2-locus,
  0 pass C3.  B: 796 fibers per N, ~85 Delta-matched partners tested,
  0 beta-match.  The deep versions (`lane_2266_fiber.m`, `lane_1212_fiber.m`,
  MordellWeilGroup per fiber + GRH class-group bounds) run in background.
- **Analytic (N,N)-gluing pipeline** (`analytic_glue.m`): periods of E,F;
  Lambda' = Lambda + graph_M/N for M in GL2(Z/N) with det = -1 mod N
  (Lagrangian precheck); form N*E_prod (integrality+unimodularity is the
  anti-isometry test); Frobenius symplectic basis -> tau -> Rosenhain ->
  Igusa-Clebsch -> rationality detection (Galois-equivariant M's only) ->
  Mestre -> twist fixed by #C(F_p) = p+1-a_E-a_F.  Pitfalls fixed: Magma's
  Sqrt(rational) returns a real (use IsSquare for exact roots);
  NEVER call ReducedMinimalWeierstrassModel on huge-coefficient Mestre output
  (it factors the discriminant and hangs) — use PARI hyperellred (Stoll
  reduction, no factoring) and reduce the IC tuple by weighted trial-division
  content first.  Validation at N=2 vs Genus2Elliptic2 in progress.
- **(N,N)-congruence treasure hunt** (`lane_congtors.m` over Frengley's
  LMFDB pair lists; usable anti-isometry classes (5,1),(7,3),(8,7),(10,1),
  (11,2),(12,11),(13,1),(14,3),(16,3)): for non-isogenous N-congruent (E,F)
  with gcd(N,#T1*#T2)=1, T1 x T2 INJECTS into the glued Jacobian's torsion —
  no descent conditions at all.  Also `lane_cong_sieve.m`: trace-congruence
  sieve mod 3/5/7 ACROSS the high-torsion family pools (E26, E28, X1(9),
  X1(10), X1(12)) — reaches conductors far beyond the 500k list limit; mod-3
  hits glue ALGEBRAICALLY via Genus2Elliptic3.  Best cases: [2,8]x[2,8] ->
  [2,2,8,8] (order 256), [2,8]x[2,6] -> [2,2,2,24] (192), [12]x[2,6] ->
  [2,6,12] (144), [2,6]x[2,6] -> [2,2,6,6] (144), [12]x[12] -> [12,12] (144).
  Theory note: for m-isogenous pairs the only Galois-equivariant
  anti-isometries are +-phi|E[N] (End_Gal(E[N]) = scalars generically), both
  degenerate — so ISOGENOUS sieve hits are useless; flagged and skipped.
- Local sanity: [2,8]-curves automatically have bad reduction at 2,3,5,7
  (16 | #E(F_p) beats the Weil bound), so no extra local obstruction to
  [2,2,8,8] beyond what [2,8]-curves already satisfy.

## 4.6 Congruence/isogeny gluing: pipeline validated, harvest running

- **Control PASS (algebraic)**: the H=60 trace sieve's only genuine
  non-isogenous congruent pair, x7(t=2/11) x x10(t=-1/2) (3-congruent, 44
  agreeing primes), glued via Genus2Elliptic3 to
  y^2 = 4x^6+36x^5-35x^4-390x^3+1237x^2-924x+4356 with exact torsion [70]
  = the full product, as the injection theory predicts.  (Known group -- but
  the chain trace-sieve -> congruence -> gluing -> product-torsion is now
  validated end to end.  The second partner x10(1/3) gives 0 curves: right
  congruence class, wrong twist-alignment; also as expected.)
- **Isogeny-class harvest** (`lane_isoglue.m`): IsogenousCurves-based scan of
  1873 classes from the family pools found 3280 candidates with product
  torsion of order >= 100 not in KNOWN: [2,6]x[12] -> [2,6,12] (144),
  [8]x[2,8] -> [2,8,8] (128, would tie the record), [10]x[10] -> [10,10]
  (100), at conductors up to 10^17.  The first (7,7)-gluing attempt on an
  isogenous [2,6]x[12] pair FOUND a Galois-equivariant anti-isometry with
  rational Igusa-Clebsch invariants (~90-digit heights) — isogenous-pair
  gluings are not automatically degenerate.  Certification (L_p-split +
  exact torsion) in progress.
- Sieve H=60 (16514 curves, 6 families incl. X1(7)): only the one x7-x10
  congruence above among non-isogenous pairs; everything else is the
  isogenous background (E26 <-> X1(12) partners are degree-2/4-connected;
  X1(10) has its own partner structure).
- **LMFDB congruent-pair lists exhausted** (`lane_congtors.m`, all 218k
  pairs in the nine anti-isometry-compatible (N,r) classes): the largest
  coprime torsion products are order 16-28 ([2,8], [2,2,4], [18], [2,14]),
  ALL known groups.  No new split-torsion group is reachable from
  conductor <= 5*10^5 congruent pairs by free product injection; the
  frontier lives at the large conductors only the family/isogeny routes
  reach.
- Kani degeneracy worked out explicitly: psi = c*phi|E[N] for an m-isogeny
  phi is degenerate iff N(a^2+m b^2) is a perfect square for some integers
  (a,b) with b = c a mod N.  For (N,m) = (3,2) both scaling classes ARE
  degenerate ((a,b) = (1,+-1) gives 9) — no algebraic (3,3)-shortcut for
  2-isogenous pairs.  For (N,m) = (5,6): a^2+6b^2 = 5u^2 has no nonzero
  solutions (3-adic descent) — such gluings are NEVER degenerate.

## 4.7 Final resolution of the (N,N)-gluing pivot (end of session)

- The "8 anti-isometries with rational invariants per isogenous pair" were an
  artifact of a too-weak rationality test: a precision-stability check
  (prec 300 vs 700) shows the recognized values TRACK the precision — the
  invariants are REAL but IRRATIONAL.  The 8/12-hit selection detects
  definability over R (M's commuting with complex conjugation on the period
  lattices), not over Q.  Lesson recorded: BestApproximation-based rational
  recognition MUST verify err < 10^-(2*height(q)+margin), which requires
  prec >= 2*height + margin; an |err| < eps test with eps ~ 1/hbound passes
  essentially every real number.
- Q-rationality needs Galois equivariance; for our classes (2-power isogeny
  graphs, End_Gal(E[5]) = scalars) the only Q-rational anti-isometries are
  the isogeny-induced psi = c*phi|E[N], available iff c^2 m = -1 mod N.
  Kani windows (rigorous, via d^2 = N(a^2+mb^2) and p-adic descents):
  (5,6), (5,24), (7,10) never degenerate; (7,6),(7,3),(7,12),(5,4),(5,9),
  (5,11),(3,2) degenerate.
- **Degree sweep over all Cremona classes to conductor 3*10^5** (degsweep):
  every same-class pair with torsion product >= 96 connects by degree
  EXACTLY 2 — never 6, 24, or 10.  Structural reason: the odd-prime isogeny
  step a golden window requires kills the torsion the product needs.
  => The (N,N)-gluing route to order >= 96 split torsion is obstructed at
  every accessible instance by the three-way tension {anti-isometry
  existence, Kani nondegeneracy, torsion preservation}.
- Session verdict on the plan's targets: no new group realized; the
  machineries, controls, impossibility proofs, the [8,8] infinite family,
  and the [70] control curve stand.  The (2,2)/(3,3)-routes (with the
  complete collapsed condition systems derived here) remain the viable
  frontier for [2,2,6,6]/[12,12]; the five non-id sigma-surfaces (each with
  rational points) are the concrete next attack.

## 5. Dead ends confirmed this session (do not repeat)

- [2,2,2,24] via [2,6]x[2,8] (2,2)-gluings: impossible (2.2).
- [2,2,6,6] via E26-self-gluings: impossible (2.3).
- [12,12]/[10,10] via conjugation self-gluing, untwisted: impossible (2.4).
- One-sided beta-square mixed routes: no carriers to H=150 (2.5).
- (Recall from plan: blind volume on any of these is hopeless.)

## 6. Post-close deep sweeps (left running / completed at session end)

- `lane_2266.m` H=200 (log `lane2266_H200.log`): 48921 values, 17711 matched
  pairs — ALL isomorphic deck-pairs again; zero genuine cross pairs for the
  [2,2,6,6] three-condition system up to height 200.
- `lane_1212.m` H=250 (log `lane1212_H250.log`, completed overnight):
  pools of 47454/57070/47558 curves (N=8/10/12); 59 more N=8 diagonal
  survivors (the a^2/(a^2+b^2) family at larger height, all exact [8,8]);
  N=10/12 diagonals empty (consistent with the rank-0 proofs); cross-pair
  beta-match survivors 0/0/0 among 35107/642/74 d-matched pairs.
  [12,12] and [10,10] cross pairs are now empty to height 250.
