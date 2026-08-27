# prod_09: (8,8) via Nicholls Lambda_334 — prescribed classes computed, stage-1 covers RATIONAL/ELLIPTIC (genus gate PASSED), production locus built

Date: 2026-07-18.  Target #9: J(Q)_tors = Z/8 x Z/8, geometrically simple.
Continuation of notes/claude_top10_08_88.md (substrate: Nicholls Prop 5.9.6, 3-parameter
(s,t,v) family, pointwise-rational (4,4) = Sigma on J1; J1 --phi1--> J2 = Jac(quintic)).
Working dir (all scripts + raw logs):
/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/t88prod/

## Strategy recap (3 lines)

(8,8) on J1 forces, via pushforward, both prescribed 2-torsion classes phi1(D_i) to be
2-divisible in J2(Q) (quintic side, full rational 2-torsion, exact x-T criterion).  We
computed the classes explicitly (task A), derived the divisibility square-systems in
Q(s,t,v) (task B), and ran the genus gate (task C): everything is genus 0/1 — GO.

## Task A: the prescribed classes (COMPUTED, not guessed)

Implemented the Richelot correspondence Gamma on C2 x C1 explicitly (corr.m):
eq1: G1(x)M1(z)+G2(x)M2(z) = 0, eq2: y*w = rho*G1(x)M1(z)(x-z), M_i = cross-brackets
[G_j,G_k]; verified the identity sum_i G_i(x)M_i(z) = k(x-z)^2 symbolically; determined
rho^2 empirically (constant, square in F_p at every good prime tried — the correspondence
is F_p-rational on Nicholls' d2-twist, confirming d2 = Richelot twist mod squares).
Transfer J1(K) -> J2(K), K = GF(p^4), via Mumford-style images (no root-finding):
u_im = eq1-quadratic/lc, v_im = rho*M1(z)/w * x*(x-z) mod u_im.
Checks all passed: homomorphism on random points, ker phi1 = {[L1],[L2],[L3]} -> 0,
mixed 2-torsion class -> nonzero, images of the global [4,4] generators are nonzero
2-torsion.  Applied to reductions of the global Sigma generators (members (2,3,1) at
p=29 and (5,2,1) at p=23,31):

  ker phi2 = phi1(Sigma) = { 0, [{0,a}], [{c,oo}], [{1,b}] }     (partition 0a | c,oo | 1b)

NOT the partition of Example 5.9.4 (0c|1,oo|ab).  (8,8) needs BOTH generators' images
divisible; since [{1,b}] = [{0,a}]+[{c,oo}], "two divisible" = "all three divisible".

## Task B: stage-1 square systems (x-T criterion, exact for divisibility)

Criterion (Schaefer, odd degree; validated 195 negative Q-tests in muval.m and 180
F_p brute-force tests incl. 12 positives in muval2.m, 0 mismatches): for a pair class T
of the quintic y^2 = h, T in 2J2(Q) iff all 5 x-T coordinates are squares.  Computing
the coordinates in Q(s,t,v) with (a,b,c) substituted (stage1.m; family squares are then
automatic) gives, with A := s^2-t^4+t^2:

- class [{0,a}]:   t^2-1 = square  AND  (t^2-1)(s^2-(t^2-1)^2) = square       (v-FREE!)
- class [{c,oo}]:  t^2-1 = square  AND  Psi := v(s^2 v+1)(A v+1) = square
- class [{1,b}]:   product conditions (automatic given the other two)

Both {0,a}-conditions are conics with rational parametrization (gate.m verified):
  t = (m^2+1)/(2m)  (alpha := (m^2-1)/(2m), t^2-1 = alpha^2)
  s = (n^2+alpha^4)/(2n)  (beta := (n^2-alpha^4)/(2n), s^2-alpha^4 = beta^2),  v free.
The {c,oo}-extra condition, after z := 1/v, is the ELLIPTIC fibration
  E_{s,t}:  w^2 = z(z+s^2)(z+A)          (full rational 2-torsion)
whose rational points z (away from 2-torsion) give exactly the valid v = 1/z.
Halves of (0,0) on E (exist iff A square) land precisely on the degenerate locus
s^2 A v^2 = 1 of the family — a cute consistency check.

## Task C: genus gate — PASSED (gate.m)

Stage-1 cover for {0,a}: iterated conics, fiber genus 0, rational 2-parameter base
(m,n) x free v = rational 3-fold.  Stage-1 cover for {c,oo} over that base: genus-1
v-fibers (cubic disc nonzero generically; mod-p checks at p=97,101,103 all genus 1).
Old chart's second-halving H_eta had fiber genus 21-31; this substrate is genus 0/1.
Verdict: GO (no Lambda_ijk fallback needed).

## Production results (exact Magma verifications)

1. STAGE-1 LOCUS VALIDATED (validate.m): 10/10 members of the (m,n,v) parametrization:
   TorsionSubgroup(J2) = [2,2,2,4] exactly, with the order-4 point doubling to the
   prescribed class {0,a} every time.  J1 stays [4,4] on these (lift layer is the
   remaining obstruction).  Data: data/claude_prod_09_88_stage1_members.txt.
2. DOUBLE-STAGE-1 IS ABUNDANT (edia.m, partial — killed for load management after 16
   bases): E_{s,t} has MW rank >= 1 on 16/16 sampled bases (ranks 1..3, e.g. (m,n)=(3,5)
   rank 3, (2,1),(3,1/3) rank 2); torsion of E always just (2,2).  So over EVERY sampled
   base there are infinitely many v with ALL THREE prescribed classes 2-divisible in
   J2(Q).  Naive v-boxes (scan2.m: 120 bases x ~1700 v, 0 hits) miss them because MW
   generator heights are large (min v-height 14575 at base (3,1/3)).
   EXACT CONFIRMATION (rank2c.m, base (3,1/3), rank 2, minimal model
   [0,-1,0,-2201424,-1046684736]): v = -729/17500, 26244/7975, 729/38425 all give
   TorsionSubgroup(J2) = [2,2,4,4] EXACTLY, with order-4 points doubling to BOTH
   prescribed classes — data/claude_prod_09_88_double_stage1.txt.  J1 = [4,4] at these
   three points (both lifts fail there; lift layer is the whole remaining (8,8) gap).
3. (8,4)-lift box scan (scan3.m): ABORTED by design review — the 10-prime local filter
   is WEAK on this locus (2-part of #J1(F_p) >= 32 automatic since J1 is isogenous to
   J2 whose [2,2,2,4] is global): survivors accumulated at ~30-50%/member (raw list in
   scan3_out.txt, m=2 slice only), so the planned exact phase was hopeless; killed by
   PID.  Correct lift testing = exact TorsionSubgroup(J1) on structured members, or the
   symbolic lift locus (next steps).
4. LIFT-LAYER DIAGNOSTIC (liftdiag.m) — key structural finding: on every one of the 10
   stage-1 members, the even-degree x-T NORM presieve (N_j = Res(monic l_j, u_D) square,
   j=1,2,3) is passed IN FULL by exactly 4 of the 12 order-4 elements of Sigma (pattern
   111; all others 001 — the l3-norm condition is identically satisfied family-wide).
   Since exact J1 torsion is still [4,4], the true lift obstruction is only the finer
   lambda-compatibility (u_D(theta_j) square in the quadratic field K_j after a common
   rational twist) — a 2-cover/conic-type condition, NOT an H_eta-type high-genus wall.
   Strongly suggests the lift locus is again low genus once derived symbolically.
5. Simplicity certificates (cert2.m / cert2_out2.txt): L-poly deg-4 irreducible AND
   12th-power transform irreducible, two primes per curve, for J2 and J1 of
   representative stage-1 members (m,n,v) = (2,1,1), (3,2,-1).
6. BANKED WIN (cert3.m): the double-stage-1 member (3,1/3,-729/17500) is CERTIFIED
   geometrically simple (p=67 and p=83) with J2 torsion [2,2,4,4] (order 64) and a
   compact minimal model y^2 = 2374249539600 x^5 + ... (see
   data/claude_prod_09_88_double_stage1.txt) — a geometrically simple (2,2,4,4)
   Jacobian living in a provably infinite structured locus (positive-rank elliptic
   fibration over every sampled base of a rational surface), unlike the isolated
   tor2244 box-search examples.

## Files (all in t88prod/ unless noted)

corr.m (+corr_out.txt) task A; muval.m, muval2.m x-T validation; stage1.m (+_out) task B
(all 15 classes' square classes); gate.m (+_out) task C; validate.m (+_out) stage-1
exact validation; scan2.m, edia.m, rank2.m/rank2b.m (aborted, non-minimal models),
rank2c.m (+_out, the one that worked) double-stage-1; scan3.m lift scan (aborted);
liftdiag.m (+_out) lift norm-presieve; cert2.m (+cert2_out2.txt) simplicity.
Repo data: data/claude_prod_09_88_stage1_members.txt (10 validated stage-1 members),
data/claude_prod_09_88_stage1_conditions.txt (all 15 classes' symbolic conditions),
data/claude_prod_09_88_double_stage1.txt ([2,2,4,4] members),
data/claude_prod_09_88_certified_curves.txt (2 two-prime simplicity certificates).

## Resume state / next steps (ranked)

0. Resume commands: everything reruns headless from the t88prod dir with
   `nice -n 10 magma -b <script>.m`.  rank2c.m is the working double-stage-1 extractor
   (MinimalModel BEFORE Generators — the non-minimal rank2.m/rank2b.m hang); at session
   end it had finished base (3,1/3) (3 exact [2,2,4,4] hits) and was still inside
   Generators for base (2,1) under a 780 s timeout; edit the `for base in [...]` list to
   add bases (rank>=1 guaranteed at all 16 sampled: see edia_out.txt).
1. Extend rank2c to more bases/more MW combinations: each new base gives fresh
   [2,2,4,4] members; run cert3-style simplicity per keeper.  J1 jackpot test is
   included in rank2c (prints "J1 UPGRADE" loudly).
2. LIFT LAYER (the remaining wall for (8,4)/(8,8)): condition is D_i in 2 J1(Q).
   Concrete symbolic route: the Kummer coordinates of D_i are computable in Q(s,t,v)
   (Nicholls Sec 5.9.2: eigenvector alpha1 v1 + alpha2 v2 of W_{T_i} with alpha from the
   split quadratic u_ijk — for Lambda_334 the relevant discriminants are squares by
   construction), and Kummer coords give the Mumford u_{D_i}(x) = x^2 - (xi1/xi0)x + xi2/xi0
   directly.  Then the even-degree x-T presieve on J1 is explicit: for each quadratic
   factor l_j of f1, N_j := Res(l_j^monic, u_{D_i}) must be a square in Q (norm-form
   necessary conditions; then the common-lambda compatibility in L*/(L*)^2 Q*).  This
   yields the lift locus as explicit square conditions over the stage-1 rational 3-fold;
   run the same genus-gate methodology on them.  If again conic/elliptic, (8,4) is in
   reach and (8,8) is a finite search on the double-stage-1 locus.
3. If lift-layer symbolic derivation stalls: large (m,n,v) scans testing J1-halving
   locally are pointless (filter weak); instead test exact TorsionSubgroup(J1) directly
   on double-stage-1 members (they are the highest-probability spots: both 4-chains
   already rational on J2).

## Unverified / caveats

- The class identification is a mod-p computation (2 members x {29},{23,31}) — but the
  answer is a discrete label constant in the family; three consistent prime/member runs
  and the passed kernel/homomorphism checks make mislabeling essentially impossible.
- The x-T criterion used is exact over Q (Schaefer kernel = 2J(Q), odd degree); over
  F_p it was validated empirically as well.
- edia.m rank bounds computed WITHOUT GRH class-group bounds for the first 16 bases
  (unconditional), except where Magma printed GRH warnings — treat printed ranks as
  lower bound 1 certificates (a rational non-2-torsion point exists) rather than exact
  ranks where the warning applied.
- scan2/scan3 grids are small; their negative results are height statements only.
