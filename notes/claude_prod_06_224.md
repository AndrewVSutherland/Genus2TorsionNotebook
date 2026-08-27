# PROD #6 (2,24): third component is a dense 2-parameter [2,12] SURFACE; mech-A/B x halvable empty at 6,426 points, no local obstruction; (2,24) still unrealized

*(Claude agent, production campaign 2026-07-18.  All scripts/outputs in session
scratchpad `t224_prod/` under
/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/.
Predecessor: notes/claude_top10_02_224.md.  New repo data:
data/claude_prod_06_224_{new37_tor24,residual62,mechA_212family,mechA_212family_deep,mechB_points}.txt.)*

## Strategy recap (3 lines)

Ran the planned per-fiber Cassels presieve at scale (179 fibers t=tn/td, tn<=24,
r-height 20000; 11 fibers to 60000-100000), exact-verified everything, and fit the
residual points — the fit is rank-full at ALL bidegrees, which led to the real
discovery: the "third component" is not a curve but a dense 2-parameter surface,
and the correct [2,24] hunt is Q4-factorization (mech A/B) sieves on that surface.

## Pipeline (all reusable, in t224_prod/)

1. `gentab.gp` + `fibersieve2.c` + `prod_exact.m` — halvability presieve
   (tables p<=199) + Magma exact stage (order-12 D, IsDivisibleBy(D,2),
   TorsionSubgroup, jackpot flag w/ multi-prime simplicity certs).
2. NEW `gentabA.gp` / `mechA_exact.gp` / `dumpA.gp` — mechanism-A sieve:
   N(k,r,z) := numerator of Q4(k) (bidegree (4,4) in (k,r), file q4num.txt;
   Q4 = f5/(X - a/z)); tables of forbidden k-residues ("quartic in r has no
   root mod p, incl. r=infinity"), same fibersieve2 binary, then exact k->r.
3. NEW `gentabB.gp` / `mechB_exact.gp` — mechanism-B sieve via the even sextic
   resolvent FB(u) = u^6+2p u^4+(p^2-4s)u^2-q^2 of depressed Q4 (deg_r 12).
4. NEW `obstructA.gp` — per-fiber, per-prime local compatibility test of
   {mech-A} x {halvable} (joint r-residue analysis incl. infinity conventions).
5. NEW `mechA_tors.m` — direct TorsionSubgroup scan (no halvability gate),
   `cert3.m` — FAST simplicity certs by direct point counting + resultant
   transform of Frobenius^12.  (CAVEAT: the EulerFactor(BaseChange(C,GF(p)))
   route used in prod_exact.m/mechA_tors.m jackpot branches HANGS on these
   models — if a jackpot ever fires, re-certify with the cert3.m method.)

## Results

### (a) Production halvability sweep (reused aborted-run state, completed here)
- 179 fibers, H=20000 (+11 fibers deep to 60000/100000): 533 unique survivors.
- Exact stage: main=358, G1=83, known=25, **NEW=37** — every NEW point has
  TorsionSubgroup exactly **[24]** (data/claude_prod_06_224_new37_tor24.txt).
  Zero [2,24].
- Residual (non-main, non-G1) set now 62 points
  (data/claude_prod_06_224_residual62.txt).

### (b) The "third component" is NOT a curve
- Nullspace fits of the 62 residual points in (Z=z^2, r): **rank-full at every
  bidegree (1..6)x(1..6)** — no polynomial component.
- 12 residual points on fiber z=25/7 alone: impossible for any r-degree<=6 curve.
- Structure found: f5 ALWAYS has the rational root rho0 = a/z (because
  disc(T+1) = z^2 and w, w' = 2(r+1)/(1 pm z) are roots of W with
  x_P - rho0 = -(1+z)^2/(4z(r+1))); halvability forces lc(f5)(x_P-rho0) in Q*^2
  — square-class conditions, not algebraic ones.  The residual halvable locus
  is (empirically) Zariski-dense in the (z,r) surface: sporadic points of a
  2-descent-cut cover, NOT a fittable curve.  Fits can never succeed; this
  closes the "third component curve fit" line permanently.

### (c) MAIN STRUCTURAL FIND: dense 2-parameter [2,12] family (mechanism A)
- Mech-A locus {N(k,r,z)=0}: a second rational Weierstrass root k gives a
  second rational 2-torsion class T1, so torsion contains <D> x <T1> = [2,12].
- The per-fiber mech-A curve is rich: 4,887 rational points found at k-height
  20000 across 179 fibers (~27/fiber), +1,533 more in the 20000..60000 band on
  17 rich fibers.  **All 6,420 points: TorsionSubgroup EXACTLY [2,12]**
  (Magma, zero exceptions, zero degenerates).
- Geometric simplicity certified on samples (t=2 r=629/100 at p=13,29,31;
  t=2 r=247/28 at p=59,73,101; t=7 r=-16/7 at p=37,47,59): L-poly irreducible
  deg 4 AND charpoly(Frob^12) irreducible deg 4.
- Data: claude_prod_06_224_mechA_212family(.deep).txt.  This is a **banked
  win**: a dense 2-parameter family of geometrically simple Jacobians with
  torsion Z/2 x Z/12 (previous repo knowledge: only the 1-param a=(1-r)/4 line).

### (d) [2,24] = (halvable) x (mech A/B): empty, and why that is meaningful
- [2,24] needs BOTH: Q4 factors (extra 2-torsion) AND D or D+T1 halvable.
  TorsionSubgroup on all 6,420 mech-A + 6 mech-B points covers BOTH branches
  (it would show [2,24]/[4,12]/[2,48] if any combination were halvable): all
  are exactly [2,12].  Conversely all 62+37 halvable points have Q4
  irreducible ([1,4] type, checked).
- Local analysis (`obstructA.gp`): **all 179 fibers are LIVE** — at every
  table prime there exist compatible r-residues; NO single-prime local
  obstruction explains the emptiness.
- First global square condition S = 2z(1+z+2k(r+1)) in Q*^2 (necessary for
  halvability at a mech-A point): passes at ~8% of mech-A points (11/142 on
  fiber 7/1) — the obstruction sits in the residual-cubic-field condition
  c(x_P - theta_3) in K3*^2, i.e. a genuinely global 2-descent correlation.
- Mech-B (rational quadratic factor, sextic resolvent sieve, 17 rich fibers,
  u-height 20000): locus is thin — 6 points total, all [2,12], none halvable.
- Triple-rational-root scan (torsion >= [2,2,12]): 0 hits among 6,420.

## Verdict

(2,24) remains unrealized.  The M(12)-halving world now looks like this:
halvable points (torsion [24]) and extra-2-torsion points (torsion [2,12])
are BOTH plentiful — two dense families on the same surface — but their
intersection is empty through 6,426 exact checks with no local obstruction.
This smells like a theorem (a Weil-pairing / Cassels-kernel argument that
mu(D) can never land in mu(J[2](Q)) when Q4 factors), and proving it would
retire the M(12) route for (2,24) entirely; failing that, the surface gives
unlimited cheap [2,12] material for glue/twist routes.

## Resume state / next steps

- All state in scratchpad t224_prod/ (tables tab_*/tabA_*/tabB_*, survivors
  surv*/, allkA.txt, mapts*.txt, logs tors_chunk*.log torsdeep_chunk*.log).
  No processes left running.
- Resume commands (from t224_prod/):
  - extend mech-A sweep: `./fibersieve2 tn td H tabA_tn_td.txt [Hlo] > out`
    then `gp -q runAdeep.gp` (edit file name), dump via dumpAd.gp, scan via
    `magma -b fn:=<pts> lo:=1 hi:=N mechA_tors.m`.
  - new fibers: add to fiblist.txt, `gp genA_all.gp`-style calls of
    genfiberA(tn,td,199,...), plus gentab.gp genfiber(...) for halv tables.
- Ranked next steps:
  1. THEOREM attempt: prove mech-A x halvable = empty on M(12) (Cassels
     pairing computation at a generic point of the mech-A surface; the
     cubic-component condition is the obstruction to formalize).
  2. If instead a hit is still hoped for: sieve the D+T1 branch DIRECTLY at
     scale (new gentab variant testing divisibility of D+T1 in J(F_p)), then
     deep k-sweeps (H up to 10^6) on fibers where the S-condition passes
     often; also widen fibers tn<=60 (tables are 0.02s each).
  3. Exploit the [2,12] family for OTHER targets: (2,2,12)/(4,12)/(2,36)-type
     hunts via dedicated codim-1 sieves on the family (triple-root sieve =
     self-intersections of the mech-A curve; halving T1 instead of D; etc.).
  4. The 37 new [24] points also feed the (2,2,8)/glue programs as source
     curves.
