# [2,24] top-10 test run: Mechanism B on G1 EXECUTED — reduces to a genus-3 plane quartic, rank >= 1, no small points, all 7 known points degenerate

*(Claude agent, 2026-07-17. Scripts/outputs in session scratchpad
`top10_02_224/`: `validate_g1_t8.m`, `g1_mechB.m/.out`, `g1_mechB_quot.m/.out`,
`fcq_model.m/.out`, `fcq_quartic.m/.out`, `fcq_arith.m/.out`,
`fcq_final.m/.out` under
/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/.
Predecessors: notes/claude_tier2_224_m24_family.md, claude_next_224_secondlocus.md,
claude_next_224c_G1_family.md; prior-session artifacts still live in
/tmp/claude-1000/-home-claude-torsion-jac/1ef8f36a-f397-479b-90a0-3e27b1b94a84/scratchpad/
next_224b/ (resolventA.m) and next_224c/ (g1_mechA.m, per-fiber tables).)*

## Dossier recap

Target (2,24), order 48 — conspicuous gap: (2,22), (2,26), (2,28) all
simple-realized; [24] itself is richly realized (17 curves + TWO rational
one-parameter families on the M(12) halving surface).  Model
`y^2 + (x-r)(T+1)y = a x^2 T(T+1)`, `T = a x^2 - x + r`, marked D of order 12,
halving D gives order 24.  Main family (r,z,a) = (-2/(1+t^2), (1+t^2)/(2t),
(1-t^4)/(16t^2)): PROVEN [2,24]-free (both extra-2 mechanisms closed by rank-0
theorems: genus-2 quotient with RankBounds (0,0), and the conductor-54 elliptic
curve y^2+xy+y = x^3-x^2-14x+29 with #E(Q)=9).  New G1 family

    z = (1+t^2)/(2t),  r = (t^2-1)^2/(8(1+t^2)),
    a = -(t^2-1)^2(1+t^2)/(2t^2(t^2+3)^2)

all sampled members torsion exactly [24].  Mechanism A on G1 (rational
Weierstrass root of the residual quartic Q4 = W/(T+1)): genus-6 cover,
genus-3 quotient FQ(u,s), t-height <= 400 scan clean, Chabauty spec'd
(g1_mechA.m), OPEN but Faltings-finite.  Mechanism B on G1 (rational QUADRATIC
factor of Q4): NEVER run before this session — the single most concrete
unfinished [2,24] step in the repo.  Also: 25 of the 51 off-main halving points
(data/m24_offmain_51points.txt) fit NEITHER family — a third component exists,
unfitted.  p=5 forces boundary behavior on all off-main fibers.

## Strategy (ranked)

1. **Third-component hunt (the only live REALIZATION route in this world).**
   Both known components are sterile (main: proven; G1: this session shows
   Faltings-finite with empty small-height locus).  A [2,24] hit, if it exists
   here, lives on the unfitted component(s) through the z = 25/7 fiber (10
   non-G1 non-main r-values there alone).  Plan: scale `code/m24_fibersieve.c`
   + `code/m24_gentab.gp` per-fiber presieves to r-height 1000-2000 on fibers
   t = 2..9 (respecting the mod-5 boundary constraint z = 0,+-1 mod 5 or
   5 | den(z)), accumulate >= 60 non-G1 points, exact nullspace fit in
   (Z=z^2, r) at bidegree (4,4)+ (the G1 derivation pattern: the fit succeeded
   because Z^2 - Z = z^2(z^2-1) is a square on every fiber), verify
   symbolically, then run BOTH mechanism covers on the new component hoping for
   a low-genus or rational-point-bearing cover.
2. **Close G1 as a theorem (this session's route, now fully mapped).**
   Mechanism B reduces to an explicit genus-3 plane quartic X_B (below) with
   rank(Jac) >= 1, torsion | 3, trivial automorphisms, likely simple Jacobian:
   needs plane-quartic 2-descent (Bruin-Poonen-Schaefer-Stoll style) for a rank
   upper bound, then Chabauty-Coleman (Tuitman's code) if rank <= 2.  Mechanism
   A closure = the already-spec'd rank/Chabauty on FQ genus 3 (g1_mechA.m).
   Both are theorem-grade tasks, not realization hopes.
3. **Orthogonal fallbacks.**  (a) [2,12] line a=(1-r)/4: halve D or D+Tind at
   r-height 300 -> 2000 with the C presieve, mod-5-aware (5|den(r) or r = 0,1,2
   mod 5) — the line was only dry to height 300.  (b) Richelot-glue route from
   the abundant [24] curves, mirroring notes/m244_to_248_route.md.

## Test run (validation + Mechanism B on G1, ~20 CPU-min, all single-threaded nice -15)

### (a) Validation — reproduces the dossier exactly (validate_g1_t8.m)

    cd <scratchpad>/top10_02_224 && nice -n 15 magma -b validate_g1_t8.m

t=8: r = 3969/520, z = 65/16, off-main confirmed; ord(D) = 12,
IsDivisibleBy(D,2) = true, TorsionSubgroup = [24]; simplicity certificate at
p=31: L_p = 961x^4 - 310x^3 + 54x^2 - 10x + 1 irreducible AND charpoly of
alpha^12 irreducible of degree 4.  Pipeline calibrated.

### (b) Mechanism B cover for G1 (g1_mechB.m — resolventA.m with G1 substitution)

Q4 = product of two rational quadratics (x^2+ux+v)(x^2-ux+w)  <=>  rational
point (u,t), u != 0, on the sextic resolvent cover

    FB:  u^6 + 2 p(t) u^4 + (p(t)^2 - 4 r0(t)) u^2 - q(t)^2 = 0

(x^4 + p x^2 + q x + r0 = monic depressed Q4).  Results:

- **Branch accounting**: u=0 branch needs q(t)=0; low-degree factors of
  numerator(q) are exactly t and t^2+3 (both degenerate: t=0, and t^2=-3 is
  den(a)=0); remaining factor irreducible of degree 28.  lc(Q4) vanishes
  rationally only at t = +-1 (degenerate).  So FB captures Mechanism B
  completely.
- Resolvent cubic irreducible over Q(t): no generic factorization.
- **FB(u,t) is irreducible over Q, u-deg 6, t-deg 80, geometric genus 13,
  even in BOTH u and t** (Q4 is even in t on G1 because r, a are even and z
  enters only via z^2 — same structural luck as the main component).

### (c) Quotient tower (g1_mechB_quot.m, fcq_final.m)

    FB(u,t)  g=13
      |  theta=u^2          |  s=t^2
    FC(theta,t)  g=7      FBq(u,s)  g=6
            \               /
          FCq(theta,s)  g=3   (irreducible; theta-deg 3, s-deg 40)

(Main-component analog of FCq had genus 1 and rank 0 — G1 is strictly harder.)

### (d) FCq analysis (fcq_model.m, fcq_quartic.m, fcq_arith.m, fcq_final.m)

- FCq is NON-hyperelliptic; canonical model is the smooth plane quartic

      X_B:  X^4 - 30X^3Y + 1344X^2Y^2 - 6384XY^3 + 34496Y^4 - 1614X^3Z
            + 13904X^2YZ - 399504XY^2Z + 807296Y^3Z + 486032X^2Z^2
            - 1183696XYZ^2 + 23657984Y^2Z^2 - 55771184XZ^3
            + 24888960YZ^3 + 2216854336Z^4 = 0.

- PointSearch bound 5000: **7 rational points** (400/3:17/3:1), (96:-13:1),
  (376:-13:1), (1144/9:43/9:1), (136:7:1), (131:9/2:1), (1912/17:59/17:1).
  Every one pulls back only into the base/singular locus of the plane model,
  (theta,s) in {(0,0), (0,-3), infinity} — i.e. **all 7 lie over the
  DEGENERATE fibers t in {0, +-sqrt(-3), inf}** (t^2+3 | num(q), den(a)).
  None is a Mechanism-B solution.
- Direct affine scan (cubic in theta at every rational s of height <= 40):
  **zero nonzero-theta rational points**; only coefficient poles at s = +-1
  (t = +-1 degenerate).  Consistent with the prior t-height <= 400 factoring
  scan (0 hits).
- #J(F_p) for the quartic: p=11: 1224, 13: 2682, 17: 6018, 19: 9210,
  23: 12132, 29: 31155, 31: 28431, 37: 49686; **gcd = 3, so
  #Jac(X_B)(Q)_tors | 3**.  Since X_B has >= 7 rational points and C(Q)
  injects into J(Q): **rank(Jac(X_B)) >= 1**.  NO rank-0 shortcut (unlike the
  main component).
- Automorphism group of X_B over Q: **trivial** (order 1) — no quotient maps.
- L-polynomial of X_B irreducible of degree 6 at all 10 primes 11..43 —
  Jac(X_B) is (almost certainly) geometrically simple of dimension 3: no
  elliptic factor to descend to.

## Verdict after test

- (2,24) remains unrealized.  **Mechanism B on G1 — the last untried concrete
  step — is now mapped: it is Faltings-finite (genus-3 quartic X_B), empty at
  all searched heights (s <= 40 here, t <= 400 prior), and every known
  rational point of X_B is degenerate.**  G1 will almost certainly never give
  [2,24], but closing it as a theorem needs plane-quartic descent + Chabauty
  (rank >= 1, tors | 3, trivial Aut, simple 3-dim Jacobian) — genuinely harder
  than the main component's rank-0 collapse.
- Realization hopes inside the M(12) halving world therefore rest ENTIRELY on
  the unfitted third component (>= 10 non-G1 points on the z = 25/7 fiber
  alone).  That is where production compute should go.

## Next steps

1. **Production**: third-component fit — scale the per-fiber presieve
   (code/m24_fibersieve.c + code/m24_gentab.gp, tables to p=199, ~89 ms/fiber)
   over fibers t = 2..9, r-height 1000-2000, exact-filter main+G1, nullspace
   fit in (Z, r) bidegree >= (4,4); then run both mechanism covers on it.
2. **Theorem**: rank bound on Jac(X_B) via plane-quartic descent (Magma:
   TwoCoverDescent-style / BPS-S; or reduce along the degree-3 map to P^1_s
   with cubic-resolvent descent); if rank <= 2, Chabauty-Coleman with
   Tuitman's QCMod/Coleman code; enumerate X_B(Q), pull back with the
   theta=square, s=square tests; combine with g1_mechA.m closure of Mechanism
   A to prove "G1 never gives [2,24]".
3. Cheap insurance: mod-25 / Q_5 obstruction analysis of the off-main locus
   (all fibers are 5-adically boundary-forced; a clean local obstruction to
   extra 2-torsion would kill [2,24] in this whole world at once).
