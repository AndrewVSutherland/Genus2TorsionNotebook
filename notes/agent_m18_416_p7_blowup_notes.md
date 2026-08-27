# FRONT 1B (continued): p=7 blowup / Hensel tower for `[4,16]`

Continuation of `agent_m18_416_search_notes.md`.  The affine and denominator
p=7 CRT/height searches were negative through height 30.  Their own conclusion
was: *"the next `[4,16]` step should not be another blind height increase; it
should be a deeper blowup/Hensel analysis of the surviving p=7 boundary
strata."*  This note does that first blowup level.

## Script

```text
code/agent_m18_416_p7_blowup.m       (log: data/agent_m18_416_p7_blowup.log)
```

It builds the cleared integer second-halving equations `E416` directly (the
`[4,16]` condition `f - ell^2 = c4*(x+R)*q^2`, `q=x^2+ax+b`, `ell=cx^2+dx+e`,
from `code/m18_m14_second_halving_equations.m`) and runs a `p`-adic blowup:

- confirm the good-open emptiness over `F_p`;
- at every boundary residue `(R0,w0)` collect all mod-`p` aux solutions
  `s=(a,b,c,d,e)`;
- for each `s` linearise: lift `R=R0+p*rho`, `w=w0+p*omega`, and require the
  aux next digit `delta` to solve `J_aux*delta = -carry(rho,omega)`.  With
  `P` = left kernel (cokernel) of the 5x5 aux-Jacobian `J_aux`, solvability is
  the linear condition `P*(c0 + gR*rho + gw*omega) = 0` on `(rho,omega)`.  The
  surviving `(rho,omega)` are the live tangent directions, i.e. the mod-`p^2`
  refinement of the boundary closure.
- if `J_aux` is invertible for some `s`, that `s` Hensel-lifts to a **smooth
  `Q_p` point** for all `(rho,omega)` — a `p`-adic `[4,16]` with bad curve
  reduction at `p`, the prime candidate for a rational point.

## Results at p=7

```text
E416 solvable residues: good_open = 0,  boundary = 21
```

So the second halving (`P_R` divisible by 2) is exactly what is empty in the
good-open chart: confirmed independently.  `21` boundary residues carry a
mod-7 `E416` solution (superset of the 19 `target416` residues; the extra 2
lack the first `[4,8]` cover).

**First blowup, mod 7 -> mod 49.**  Of the naive `21*49 = 1029` mod-49 lifts,
`449` survive:

```text
residue <3,0>          : 0/49   <-- KILLED (first-order local obstruction)
residues <0,2>,<0,5>   : 2/49   (almost killed)
<1,2>,<1,3>,<1,4>,<1,5>: 7/49
<0,3>,<0,4>,<5,2>,<5,5>: 24/49
<3,3>,<3,4>            : 23/49
<1,1>,<1,6>,<6,1>,<6,6>: 25/49
<4,0>                  : 28/49
<0,0>,<1,0>,<5,0>      : 49/49  (blowup flat in (rho,omega); all w0=0)
TOTAL surviving (R,w) mod 49 = 449
```

**Key structural finding:**

```text
Hensel-smooth boundary residues = 0
```

At *every* boundary residue, *every* mod-7 aux solution has a **singular**
aux-Jacobian.  Hence there is **no smooth `Q_7` point** of the `[4,16]` cover
reducing to the `p=7` boundary.  Any rational `[4,16]` in `M_1(8,4)` must
reduce to a *singular* point of the `E416` cover mod 7 and can only be reached
through an **iterated** blowup — a single Hensel step never suffices.  This
explains cleanly why the height-30 affine/denominator CRT searches found
nothing: the relevant lifts, if they exist, live deep in the blowup tower, not
at bounded height with generic reduction.

## What this buys us

1. A genuine partial obstruction: stratum `<3,0>` is dead at first order, and
   `<0,2>,<0,5>` are all but dead (`2/49`).  The live closure shrinks
   `1029 -> 449`.
2. The `449` surviving `(R,w) mod 49` are a strictly stronger CRT gate than the
   21 mod-7 residues for any resumed height search
   (`agent_m18_416_search_crt.m`); use them to prune before exact tests.
3. The three flat strata `w0=0` (`<0,0>,<1,0>,<5,0>`) are where the first-order
   blowup is degenerate; they are the natural targets for the **second** blowup
   level (mod 343).

## Second blowup level (mod 49 -> mod 343)

```text
code/agent_m18_416_p7_blowup_level2.m
```

Resolves the three flat `w0=0` strata by explicit Hensel lifting of the exact
`E416` system through mod 343 (per (R,w) mod 343, lift each mod-7 aux solution
digit by digit, enumerating the aux nullspace coset).  Aux degeneracy measured
first: most boundary strata have tiny aux fibers (`#aux = 2..6`, `rank J_aux =
3..4`); the exception is `<1,0>`, where `J_aux == 0` identically (nullity 5).

Results (of the `2401` mod-343 lifts reducing to each residue):

```text
<0,0>: survive 2401 / 2401   (fully live -- blowup stays flat)
<5,0>: survive 1715, KILLED 686   (partial second-order obstruction)
<1,0>: all 2401 CAPPED (J_aux==0, nullspace too large to enumerate; unresolved)
```

**Interpretation.**  The blowup does NOT produce a full 7-adic obstruction.
Instead it concentrates the 7-adic `[4,16]` locus:

- `(R,w) == (0,0) mod 7` is fully live to mod 343 -- i.e. `7 | R` and `7 | w`.
  This is a genuine live stratum (curve with bad reduction at 7), and it is
  exactly the region the affine height-20/30 CRT scans barely sampled (only
  `R = 7,14,...`, `w = 7,14,...`).  It is the prime refined search target.
- `(5,0)` is `~71%` live; its 1715 surviving mod-343 residues gate a search.
- `(1,0)` is deeply degenerate (`J_aux==0`); the aux variables enter `E416`
  only quadratically mod 7 there, so it needs a dedicated Hessian / weighted
  blowup, not the linear Newton step.

## Hessian / weighted blowup of `<1,0>`

```text
code/agent_m18_416_p7_hessian_10.m   (log: data/agent_m18_416_p7_hessian_10.log)
```

Resolves the `J_aux == 0` stratum using the exact quadratic structure:
`E416` is exactly quadratic in `aux=(a,b,c,d,e)`, so
`E(aux0+7*sigma) = E(aux0) + 7*Grad.sigma + 49*Qf(sigma)` exactly, with the
pure quadratic part `Qf` constant on the stratum (mod 7).  With
`Grad(aux0) == 0 mod 7` (verified for all 49 aux solutions):

- mod-49 solvability is a condition on `(R,w) mod 49` alone: **all 49**
  sub-residues survive, with 385 `(pair, aux0)` branches;
- mod-343 needs `sigma` with `E/49 + G'.sigma + Qf(sigma) == 0 (mod 7)`;
  the third-digit dependence is linear, so one sigma-scan per branch covers
  all 49 third digits.

Result:

```text
STAGE2 mod-343: live (R,w) residues = 343 / 2401, depth-2 Hensel-smooth = 0
```

So `<1,0>` is cut down by a factor 7 (2058 of 2401 residues KILLED -- exactly
one linear condition on the third digits), and it still carries **no depth-2
Hensel-smooth branch**: even after the weighted blowup there is no smooth
`Q_7` point over `<1,0>`.  The everywhere-singular character of the `[4,16]`
cover over `Z_7` persists through depth 2 on this stratum.

## Refined live-stratum search (mod-49 CRT gate)

```text
code/agent_m18_416_live_stratum_search.m
data/agent_m18_416_live_stratum_h70_part{0..5}.log
```

Replaces the 21-residue mod-7 CRT gate of `agent_m18_416_search_crt.m` by the
449-residue mod-49 survivor gate (recomputed inline), keeps the `11,13`
good-open aux filters and the exact chain
(`FirstCoverPossible -> TangentCandidates -> exact Tx half -> exact P_R half
-> torsion`), and pushes the rational height to 70 (previous passes: 20/30),
partitioned across 6 jobs.  Height-20 validation run reproduces Filip's
counts under the stricter gate.

Height-70 results (6 parts combined):

```text
parameters 5975, mod49-gated pair budget 5,109,378
candidates_after_aux 3,092,096   aux_killed 2,012,610
family_smooth 3,079,864
first_possible 233
tangent_bases 23
exact_tests 23  first_verified 23  (all Tx halves verified)
PR_halved 0    hits 0
```

Negative through height 70 on the full mod-49-gated live closure: rational
first halves (`[4,8]` tangent points) keep appearing, but the second halving
of `P_R` never has a rational point.  Together with the smoothness results
below, this is consistent with a structural `Q_7`-level obstruction rather
than a height artifact.

## Smoothness scan: implicit-function reduction (MAJOR UPDATE)

```text
code/agent_m18_416_p7_smooth_scan.m
data/agent_m18_416_p7_smooth_scan_all_l2.log
data/agent_m18_416_p7_smooth_scan_00_l4.log
```

For aux solutions with `r = rank(J_aux mod 7)` in `1..4`, Hensel-solve `r`
pivot aux variables (invertible `r x r` minor) as exact `Z_7`-analytic
functions, leaving a reduced system `G` of `5-r` equations in
`(R, w, u_1..u_{5-r})`.  A node where `grad(G) mod 7` has full rank `5-r`
carries genuine smooth `Q_7` points of the `[4,16]` cover (implicit function
theorem).  Cross-check: the free-aux columns of `grad(G)` vanish at level 1
(Schur complement of an exact-rank block), and the `(R,w)` columns reproduce
the level-1 cokernel conditions.

**Result: the `[4,16]` cover has genuine smooth `Q_7` points.**  The earlier
"no Hensel-smooth point" statements referred to the naive fixed-`(R,w)` test
(rank-5 `J_aux`), which is too strong; letting `(R,w)` move, six strata are
transverse already at level 1:

```text
SMOOTH (Q_7 points exist): <3,3>, <3,4>, <4,0>, <5,0>, <5,2>, <5,5>
  (via rank-4 aux branches, and some rank-3 branches with rank-2 (R,w)-block)
```

**There is NO 7-adic obstruction to `[4,16]` in `M_1(8,4)`.**  Since the
good-open target416 residues at `p=11,13` are nonempty (4 each), the cover is
plausibly everywhere locally soluble.  The absence of rational points through
height 70 is therefore a height phenomenon or a GLOBAL obstruction, not a
local one at 7.

## Final p=7 stratum map for [4,16]

```text
SMOOTH (Q_7 points):      <3,3> <3,4> <4,0> <5,0> <5,2> <5,5>
dead (level 1):           <3,0>
nearly dead (2..24/49):   <0,2> <0,5> <0,3> <0,4>   (rank 3; no smooth
                          node through level 2, exhaustive)
persistent degenerate:    <0,0>  rank 3 but grad(G) == 0 mod 7 at EVERY node
                          through level 3; survival cut ~1/7 per level in the
                          u-directions while the (R,w)-projection stays FULL
                          (2401/2401 mod 343, matching the level-2 script).
                          This signature (never-smooth + persistent survival +
                          clean factor-7 cuts) suggests the residual system
                          has a MULTIPLE (square) structure along the branch;
                          if G = unit * H^2 + O(7^k), the true equation is H.
rank-0 everywhere (J_aux == 0, need Hessian treatment):
                          <1,0> (done: 343/2401 live mod 343, no depth-2
                          smooth), <1,1> <1,2> <1,3> <1,4> <1,5> <1,6>
                          <6,1> <6,6> (all aux branches rank 0; the corners
                          <1,1>,<1,6>,<6,1>,<6,6> have E416 == 0 identically
                          in aux mod 7, i.e. 16807 solutions)
```

Note the +-ell symmetry: aux solutions come in pairs `(a,b,c,d,e)` and
`(a,b,-c,-d,-e)` (replace `ell` by `-ell`); the aux fibers above reflect it.

## Symbolic reduction: the cover as a surface in A^4

```text
code/agent_m18_416_e1e0_reduction.m
code/agent_m18_416_branch_discriminant.m
```

Eliminating `ell = (c,d,e)` exactly: `P_R` halves iff the quartic
`D := f - c4*(x+R)*q^2` is a perfect square `ell^2`.  For `deg D = 4` this
is the two classical conditions

```text
E1 = 8*d4^2*d1 - d3*(4*d4*d2 - d3^2) = 0,
E0 = 64*d4^3*d0 - (4*d4*d2 - d3^2)^2 = 0,
```

plus the squareclass condition `d4 = square`.  So the `[4,16]` cover is a
SURFACE `Sigma' : E1 = E0 = 0` in `A^4_{R,w,a,b}` (E1: deg 22, E0: deg 30;
both with only `(R-1)`-power boundary factors).  At the origin (`<0,0>`
branch) both cores have lowest form `b^2` (times `-8R`): the surface is a
DOUBLE BRANCH in `b` -- the explicit square structure behind the 7-adic
never-smooth tower.

`Res_b(E1core, E0core)` (computed mod 101/103/107; stable) factors as

```text
(w^2-1)^3 * (R+1)^8 * F5(R,w,a)^8 * F52(R,w,a)
```

with `F5` of degree 5, LINEAR in `a`.

## The degenerate-ell branch and its plane curve

```text
code/agent_m18_416_collision_family.m
code/agent_m18_416_degenerate_branch.m
code/agent_m18_416_sigma_hunt.m
```

Rational reconstruction + symbolic verification shows `F5` is exactly the
numerator of `d4` divided by `(R-1)`: the multiplicity-8 factor is the locus
where `D` DROPS DEGREE.  On it `E0 = -d3^4` forces `d3 = 0` as well, so the
branch is
`d4 = 0 (solves a = a(R,w)), d3 = 0 (solves b = b(R,w))`, leaving
`D = d2*x^2 + d1*x + d0` and the degenerate halving condition `D = (linear)^2`:

```text
Sigma(R,w) := d1^2 - 4*d2*d0 = 0   and   d2 = square.
```

`Sigma` factors as `(R-1)^2` times an irreducible plane curve of degree 40
(degw 32, even in w).  Verified: on `Sigma = 0` with `d2` square the
identity `f - ell^2 = c4*(x+R)*q^2` holds with `ell` linear -- each rational
point gives a rational half of `P_R` (order-16 point, torsion >= [2,16]).

Hunt results (`agent_m18_416_sigma_hunt.m`):

```text
mod-p: core points exist at every p in 11..23, nearly all with d2 square
       (no local obstruction on the branch);
rational: R-sweep height 30 -- ZERO rational points (not even Sigma-points).
```

The degree-40 curve is high genus, so its rational points are finite and
special; a deeper sweep (height 150) is running.  Note this branch is NOT
the `<0,0>` tower (the tower has `a == 0 mod 7`; the branch forces
`a == 2 mod 7` at the origin), so the tower's square structure lives on the
main sheet `F52` and remains to be extracted.

## Main-sheet analysis via exact 2-descent (the right formulation)

```text
code/agent_m18_416_descent_conditions.m    (+ log h20)
code/agent_m18_416_descent_param_search.m  (+ logs h50v2 part0..5)
code/agent_m18_416_sstage_obstruction.m    (+ log h15v2)
```

The fiber of `Sigma'` over `(R,w)` has degree 16 = #J[2]: it is the halving
torsor of `P_R`.  For the odd-degree model the halving condition is EXACT
2-descent (Schaefer, `x - T` map, odd degree => ker delta = 2J(Q)):

```text
P_R in 2J(Q)  <=>  u = X_R - T  is a square in Q[T]/F,
F = T*At*Bt (monicized),  X_R = -c4*R,  c4 = R+2+4t = 2*(R^2-1)/(w^2-1).
```

Componentwise:

```text
C1:  -c4*R      = sq   (rational component)
C2:  A(-R)      = sq   (norm at the A-quadratic; = -2*R*(R-1)^2*Qfac/(w^2-1))
C3:  c4*B(-R)   = sq   (implied: C1*C2*C3 = c4^4*Y_R^2 * even factors)
S_A, S_B: u actually a square in each quadratic component
          (2*(alpha +- n) = sq for one sign; the two sign-values multiply
           to d * square).
```

VALIDATED exactly: 30/30 rational samples agree with `IsDivisibleBy(P_R,2)`,
and 84/84 finite-field samples over `F_11, F_13` including 12 positive
halvable controls.

**Beware Magma `Factorization` dropping unit constants**: `c4` carries a
factor 2 and `C2` a factor `-2`; both bit us once (fixed; the h50 v1 logs
are invalid, v2 logs are correct).

**Pell trick.**  With `W = w^2`, `K = -2R(R^2-1)`:
`C1 & (w rational) <=> w^2 - K*s^2 = 1`, a conic with the trivial point,
giving the FREE parametrization `w = (m^2+K)/(m^2-K)`.  On it C2 collapses
to one square condition

```text
G(R,m) = 2*(R^2-1)*(R*(2R+1) - W*(R+2)) = sq,   W = w(R,m)^2.
```

**Results (h50 v2, 6 parts):**

```text
9.57M (R,m) pairs -> 4.78M distinct C1-points (effective w-heights ~10^3)
-> 5,632 C1&C2 points -> SECOND STAGE KILLS ALL 5,632.
```

The second stage (S_A, S_B) is a wall: 0 passes ever, matching every prior
negative but now at vastly larger effective height and with the failure
localized precisely.  Diagnostics (`agent_m18_416_sstage_obstruction.m`):
the ramification of the failing quaternion classes `(2(alpha+n), d)` VARIES
from point to point -- there is no single constant ramified place, so the
obstruction (if structural) is a relation among the classes, not a constant
Brauer class.  Useful exact identities for the follow-up:

```text
2(alpha+n) * 2(alpha-n) = 4*d*beta^2  == d  mod squares;
N_A(R,m) = G(R,m)*h(R,m)^2 for an explicit rational h
  => on the double cover g^2 = G, n_A = +-g*h is an explicit function,
     and S_A becomes: 2*(alpha_A(R,m) +- g*h(R,m)) = sq
  -- an explicit iterated 2-cover of the conic bundle {g^2 = G}.
```

## Symbolic second-stage tower (completed)

```text
code/agent_m18_416_sstage_symbolic.m       (log: ...sstage_symbolic.log)
code/agent_m18_416_sstage_crossproducts.m  (log: ...sstage_crossproducts.log)
```

On the Pell-parametrized C1 locus everything is explicit over `Q(R,m)`
(`P3 := R^3 - R + m^2/2 = (m^2-K)/2`):

```text
N_A/G and N_B/G are BOTH perfect squares:
  h_A = P3^3/(m^3*R*(R+1)),   h_B = P3^2/(m^2*R*(R+1)),
so on the cover g^2 = G both second-stage norms trivialize at once and
  n_A = g*h_A,  n_B = g*h_B.
beta = +-1/2 exactly, giving the EXACT identities
  V_A+ * V_A- = d_A,   V_B+ * V_B- = d_B,   V_*pm = 2*alpha_* +- 2*g*h_*.
S_A <=> V_A+ in {1, d_A} mod squares;  same for B.
```

**No function-level obstruction exists.**  Squareness of `X + Y*g` in
`Q(R,m)(g)` forces its norm `X^2 - G*Y^2` to be a square in `Q(R,m)`;
`N(V_A) = d_A`, `N(V_B) = d_B`, `N(V_A V_B) = d_A d_B`, and base twists
change norms by squares.  All three classes are NONtrivial, so no twist
of any `V` or product is a function-field square: the 5,632/5,632
second-stage kill is genuinely arithmetic (point-dependent), not a
constant Brauer class.  Consistently, the failing quaternion ramification
varies from point to point.

**Reduction to per-fiber 2-covers.**  Eliminating `g`, the S-conditions
become explicit plane curves over the `m`-line for each fixed `R`:

```text
S_A-cover:  y^4 - 4*alpha_A(m)*y^2 + d_A(m) = 0   (deg_m alpha_A = 8/4, d_A = 16/8)
S_B-cover:  y^4 - 4*alpha_B(m)*y^2 + d_B(m) = 0   (deg_m alpha_B = 4/2, d_B = 4/2)
```

`P_R` halvable over a rational point of the fiber `E_R : g^2 = G(m)`
(genus 1)  <=>  a rational point on the fiber product of the two covers.
The B-cover is small (bidegree (4,4)-ish); the A-cover is the big one.

## Per-R local solvability scan and ELS fiber attack (completed)

```text
code/agent_m18_416_fiber_local_scan.m   (log: ..._fiber_local_scan_h30.log)
code/agent_m18_416_els_fiber_attack.m   (log: ..._els_fiber_attack.log)
code/agent_m18_416_els_mw_deep.m        (logs: ..._els_mw_deep_R25_4/R8/R29_8.log)
```

Scan of 60 C1&C2 fibers at p <= 23 and infinity (representatives mod p^2
plus valuation shifts):

```text
57 fibers locally obstructed -- at VARYING places (2:36, 3:28, 7:23,
   5:13, 23:11, 11:9, 19:7, 13:4, 17:4); no uniform place, consistent
   with the no-constant-class theorem.
3 fibers everywhere locally solvable:
   R = -8    (E_R rank 1, torsion [2,2])
   R = -25/4 (E_R rank 3, torsion [2,2])
   R = -29/8 (E_R rank 1, torsion [2,2])
```

Deep Mordell-Weil enumeration on the ELS fibers (full multi-generator
combos, exact second-stage descent on every accessible fiber point):

```text
R = -25/4: |n| <= 4 over 3 generators x [2,2]: 1,658 distinct m -- 0 passes
R = -8:    |n| <= 25: 64 distinct m (height-capped) -- 0 passes
R = -29/8: |n| <= 25: 24 distinct m (height-capped) -- 0 passes
```

**Interpretation.**  Everywhere-locally-solvable halving covers on
positive-rank fibers whose accessible global points ALL fail is the
classic signature of a nontrivial Sha[2]-type class: the S-condition
2-covers of E_R appear to be ELS torsors without rational points.
Combined with the no-function-level-identity theorem, the picture is:
the [4,16]/P_R-halving obstruction in M_1(8,4) is a genuinely global,
fiberwise Sha-flavored phenomenon -- invisible to local analysis at any
single place, and to any constant Brauer class.

## Next steps

1. **Certify the fiber obstruction.**  On R = -25/4 (rank 3, richest),
   run a genuine descent on the S_B 2-cover (the small quartic
   y^4 - 4*alpha_B(m)*y^2 + d_B(m) = 0 pulled back to E_R): e.g.
   TwoCoverDescent / Cassels-Tate pairing to prove the torsor class is a
   nontrivial element of Sha(E_R)[2] (or find a large point).  A proof on
   one fiber demonstrates the mechanism; uniformity across fibers would
   be the route to full nonexistence in this chart.
2. ELS scan over MANY more fibers (758 known at h30; scan was capped at
   60): more ELS fibers = more chances one has an accessible point.
3. Same descent formulation for the `[2,2,16]` `m3222` chart (applies
   verbatim).
