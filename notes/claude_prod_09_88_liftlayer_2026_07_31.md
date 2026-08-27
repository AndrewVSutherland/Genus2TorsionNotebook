# [8,8] lift layer: symbolic sections, closed-form lambda classes, genus gate = UNIFORM GENUS 3 (2026-07-31)

Primary lane of the 2026-07-31 campaign; continuation of
`notes/claude_prod_09_88.md` (whose /tmp scratchpad is lost — the family was
first REBUILT in-repo and re-validated on all 15 recorded ground-truth rows:
`code/claude_prod_09_88_defs.m` + `code/claude_prod_09_88_family.m`,
log `results/claude_prod_09_88_family_rebuild.log`).

## 1. Numeric lift survey (part A) — liftdiag reproduced + new fact

`code/claude_prod_09_88_liftlocus.m` (log `..._liftlocus_partA.log`), at the
10 stage-1 + 3 double-stage-1 members (J1 = [4,4] exactly at all 13):

- liftdiag reproduced exactly: at every stage-1 member, exactly 4 of the 12
  order-4 elements pass the full norm presieve (pattern 111), the l3-norm
  never fails, and **all four passers double to T3 = [l3] = [x^2-a]**.
- NEW: at the three double-stage-1 members, **all 12 order-4 elements pass
  the norm presieve** — over the double-stage-1 locus the lambda-compatibility
  is the entire remaining (8,4)/(8,8) obstruction.
- 0 anomalies: no element has a common lambda (consistent with exact [4,4],
  and with the x-T criterion being effectively exact here).

## 2. Closed-form lambda classes (new this session)

For val_j = u_D(theta_j) = Pt + Qt*theta' with theta'^2 = e_j (complete the
square in l_j), a rational twist lambda with lambda*val_j square in K_j
exists iff w^2 := Pt^2 - e_j Qt^2 = Norm(val_j) is a rational square, and
then the only two candidates mod squares are
    lambda_j^{+-}  =  2*(Pt +- w),
whose product is e_j mod squares.  So "D in 2 J1(Q)" (up to the even-degree
kernel caveat, which is harmless for a search prefilter) becomes: the three
2-element sets {lambda_j^+, lambda_j^-} share a common square class — i.e.
two explicit product-square conditions.  No Kummer/W_T machinery needed.

## 3. Symbolic sections (part B)

The order-4 sections over T_j solve beta^2 - g1 = kappa u^2 l_j; evaluating
at the roots of l_j forces monic(l_j) | beta, so beta = monic(l_j)(mu x+nu)
and the identity divides down to "explicit quartic in x = kappa*(monic
quadratic)^2": two closure equations F1 (deg 6), F2 (deg 8) in (mu,nu).
Implemented per stage-1 base over Q(v) in `code/claude_prod_09_88_liftgate.m`:
resultant, factor, rational roots => sections u_D = x^2 + p(v)x + q(v).
VALIDATED: at bases (2,1) and (3,1/3) the specialized sections at v=1 match
the exact torsion computation's u-polynomials verbatim (assert passes).
(The full trivariate Q(s,t,v) version `code/claude_prod_09_88_liftsym.m`
was killed at ~1h in the resultant — superseded by the per-base pipeline.)

Structure found: over each tested base there are exactly 2 rational
sections over T3 (the +- pairs seen numerically), and for BOTH sections all
three norms Norm(val_j) are IDENTICALLY SQUARE in Q(v) — w_j is a rational
function, so lambda_j^{+-} are explicit rational functions of v.

## 4. GENUS GATE: the lift cover has v-fiber genus 3, uniformly

The lift locus over a base is the set of rational points of
  X: y12^2 = lambda_1^{e1} lambda_2^{e2},  y13^2 = lambda_1^{e1} lambda_3^{e3}
(2 sections x 8 sign combos = 16 components per base; squarefree models have
deg12 = 6, deg13 = 4).  Exact genus over Q, computed via function fields:

| base | components | genus |
|---|---|---|
| (m,n) = (2,1)   | 16/16 | **3** |
| (m,n) = (3,1/3) | 16/16 | **3** |
| (m,n) = (2,3)   | 16/16 | **3** |

Logs: `results/claude_prod_09_88_liftgate_m2n1_l3.log`, `..._m3n13_l3.log`,
`..._m2n3_l3.log`.

**Decision (gate branch genus >= 2): the free-v lift route on Lambda_334 is
a wall** — but a well-shaped one: for every rational base (m,n) the (8,4)
candidates are the rational points of an EXPLICIT genus-3 curve (Faltings-
finite per base).  Contrast the old M_1(8,4) chart's H_eta wall (fiber genus
21–31): this substrate got the wall down to genus 3, but not below the
Faltings line.  Directions l1/l2 on the free v-line fail already at the norm
layer (not identically square — extra cover layers), consistent with part A.

## 4b. Base-grid scan results (same day)

The lambda-compatibility scan ran twice over the full 96-base grid
(m in {2,3,4,5,1/2,1/3,2/3,3/2,5/2,2/5,4/3,3/4} x n in
{1,2,3,5,1/2,1/3,2/5,5/2}), all clean (no errors, degenerate bases, or
timeouts):

| v-height bound | bases | LIFTHITs |
|---|---|---|
| 60  | 96/96 | 0 |
| 250 | 96/96 | 0 |

So the T3-direction lift locus has NO rational point with v-height <= 250
over any of these bases — a height statement across ~1500 explicit genus-3
cover components (master log
`results/claude_prod_09_88_liftscan_master.log`).  The j1exact probe
(independent route: exact J1 at MW-generated double-stage-1 members) ran in
parallel: through 5 completed bases + (3,1/3) previously, ~90 exact
TorsionSubgroup(J1) computations, all exactly [4,4]; base (4,1) timed out
in MordellWeilGroup (revisit with RankBounds-only if desired).

FINAL (17:08): the j1exact sweep completed — 12 bases attempted, 9 finished
(bases (4,1), (3,4), (2,1) timed out in MordellWeilGroup at 90 min each),
**164 exact TorsionSubgroup(J1) computations at MW-generated double-stage-1
members, every one exactly [4,4], zero upgrades**
(results/claude_prod_09_88_j1exact_master.log + per-base logs).  Both
independent routes (lambda-scan to v-height 250, exact-J1 probe) are now
closed negatives at their stated heights; the lane's frontier is deeper MW
combinations, the timed-out bases, and the Tables-5.2/5.3 Lambda branches.

## 5. What continues to run / next moves (ranked)

1. **Base-grid lambda-scan** (`code/claude_prod_09_88_liftscan_run.sh`,
   liftgate.m `Scan:=1`): sweep ~96 bases, testing rational v of height<=60
   against the exact lambda-compatibility (microseconds per v).  Any LIFTHIT
   is an (8,4) candidate for exact `TorsionSubgroup(J1)`.  This is the
   sporadic-point hunt on the genus-3 family — the same shape as the ABC
   search that found [2,2,2,12] curve #3.
2. **j1exact probe** (`code/claude_prod_09_88_j1exact_run.sh`, running):
   exact J1 torsion at MW-generated double-stage-1 members across 12 bases —
   the highest-probability jackpot spots; prints `J1 UPGRADE` loudly.
3. **Other Lambda_ijk branches**: Nicholls Tables 5.2/5.3 (C. Nicholls,
   *Descent methods and torsion on Jacobians of higher genus curves*,
   DPhil thesis, University of Oxford, 2018; consult the thesis directly).  The halving/lambda machinery here
   applies verbatim to any parametrizable branch; the lift-cover geometry
   (genus) can differ per branch.  Also note thesis Remark 5.9.7 ("we
   haven't been able to solve Delta_1 = v_1^2 in general") — the unsolved
   discriminant condition is closely related to this lift layer.
4. If a LIFTHIT verifies as exact [8,4]: certify + bank, then the (8,8)
   question at that member is the SECOND direction's lambda condition at the
   same point — evaluate immediately.

## Caveats

- The lambda-compatibility is derived from the x-T (Cassels) map on an
  even-degree model; its kernel can exceed 2J(Q) by index up to 2, so the
  condition is used as a NECESSARY prefilter only — every hit goes to exact
  `TorsionSubgroup(J1)`.
- Genus-3 verdicts are per-base computations at 3 bases; a base with
  degenerate lambda functions (component collapse) could in principle have
  lower genus — the base-grid scan implicitly probes this (degenerate bases
  surface as errors/TIMEOUTs in the master log; inspect any such base).
