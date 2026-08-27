# [4,16] production campaign on M_1(8,4): full per-R funnel b<=100, 627 fiber-closure theorems, global conic disjunction

*(Claude production agent, 2026-07-18.  Target #3 of the top-10.  Session
scratchpad `t416/prod/`: `killsets.m`+`kill_p*.txt`, `tier1.gp`,
`tier1_b*.txt`, `postproc.gp`, `rank0.gp`+`rank0_c[1-4].txt`, `closure.gp`,
`pass2.gp`, `pass2post.gp`, `cert_bank.m`, `sweep1e6.sh`+`tier2_b*.txt`.
Repo deliverables: `data/claude_prod_416_killsets_p47_199.txt`,
`data/claude_prod_416_rank0_fibers.txt`, `data/claude_prod_416_new_points.txt`.
No repo file modified.  An earlier aborted launch of this campaign had already
completed the kill tables and the tier-1 sweep; this session reused them.)*

## Strategy recap (3 lines)

Per-R elliptic solve on the `(R,w)` chart of M_1(8,4) (notes/claude_top10_03_416.md):
membership in the family is rational-point-finding on two genus-1 quartics
`C_R^+/-` per fiber; candidates are filtered by sound all-twist kill tables
(now p<=199) and finished by the exact Magma stage; independently, `ellrank`
rank-0 fibers + torsion bookkeeping yield unconditional per-R nonexistence
theorems, and a new descent computation collapses all exponent-16 routes into
two explicit conic conditions.

## What ran (exact commands)

```text
cd <scratchpad>/t416/prod
./runkills.sh                 # (A) killsets.m for p=47..199, 32 primes, ~2.5h
                              #     -> kill_p47.txt .. kill_p199.txt, 411392 KILL residues
./runtier1.sh                 # (B) tier-1: all 6086 fibers R=a/b, 0<|a|<b<=100,
                              #     gcd=1 (R->1/R involution halves work),
                              #     hyperellratpoints on C_R^+/- at w-height 1e5
                              #     -> 22898 candidate points (tier1_b*.txt)
KMAX=199 OUTPRE=post199 gp -q postproc.gp   # kill filter, 42 primes
magma -b infile:=exact_in_prod.txt ../../code/m18_m14_416_exact.m
INFILE=fibers_all.txt I0=.. I1=.. OUT=rank0_cN.txt gp -q rank0.gp   # (C) x4 chunks
gp -q closure.gp              # closure bookkeeping on both-rank-0 fibers
gp -q pass2.gp                # conic-stratum octic solve, all fibers, Hm=2000
gp -q pass2post.gp            # kill filter on pass-2 points
magma -b cert_bank.m          # simplicity certificates for new points
RANGES="2:40 41:60" TAG=j1 ./sweep1e6.sh    # (B cont.) tier-2 w-height 1e6
RANGES="61:75 76:84" TAG=j2 ./sweep1e6.sh   # (j3 = b85:100 stopped for load)
```

## Results

1. **No [4,16] / [2,4,16].  The funnel is clean at b<=100, w-height 1e5
   (tier 1) and [partially, see resume] 1e6 (tier 2).**  22898 tier-1
   candidates -> 14 survive the p<=199 kill tables -> all 14 are DEGENERATE
   boundary points on the line `Lplus: w=(3R+1)/(R+3)` (discB=0, f not
   squarefree; they appear with y!=0 on `C_R^-` and can never be killed
   because boundary residues are never emitted as KILL).  The exact stage
   rejects all 14 (`read 14 smooth 0`).  Zero genuine candidates survive.
   Kill histogram of the 22884 killed: 11:3155, 13:5251, 17:5587, 19:2886,
   23:1692, 29:1920, 31:1034, 37:661, 41:355, 43:204, 47:82, 53:42, 59:7,
   61:5, 67:1, 71:1, 73:1 — no kill needed a prime > 73.  Only 16 of the
   22898 candidates lie on the coset-1 conic stratum (all killed by
   p<=31): the conic is razor-thin in the w-sweep, which is why the
   m-parameterized pass-2 (result 5) is the right way to search it.

2. **(A) Kill tables extended to p<=199** (32 new primes, 411392 residues,
   `data/claude_prod_416_killsets_p47_199.txt`, generator `killsets.m`
   audited criterion; p=11 spot-run byte-identical to the repo base file).
   Potency check: the old height-5093 in-family point `(1/11,-1843/5093)`
   (exact: [4,8], div false) is now killed by p=47 — the extended tables
   independently re-derive the exact-stage negatives.

3. **(C) 627 unconditional per-R nonexistence theorems**
   (`data/claude_prod_416_rank0_fibers.txt`).  All 6086 fibers were run
   through `ellrank` on both elliptic models of `C_R^+/-` (PARI r2 upper
   bound = Selmer - torsion - Cassels-pairing rank: unconditional; zero
   failures/timeouts).  627 fibers have r2=0 on BOTH branches; for these
   the full rational point set of both quartics is the (computed) torsion,
   and closure bookkeeping (`closure.gp`) shows every torsion-image w is
   either degenerate or killed by a p<=199 residue.  Hence: **for these 627
   values of R (and their 1/R mirrors: 1254 fibers), no rational w
   whatsoever gives exponent-16 torsion** — a statement no height-bounded
   sweep can make.  Rank profile of all 6086 fibers in `rankagg.txt`
   (scratchpad): 627 both-0, rest have a positive-rank branch.

4. **(D) NEW STRUCTURE: the exponent-16 locus lies on two explicit conic
   strata; all four J[2]-twists share ONE conic.**  Machine-checked
   identities on the chart: `lc(B) = c4 = 2(R^2-1)/(w^2-1)`, `B(0) = w^2 c4`,
   `A(0) = R^4`.  Consequences for the 2-descent (x-T) map's Q-component
   (first component of `L = Q x Q[T]/(g)`), writing eps1 for its value mod
   squares:

   - `eps1(H_x) = c4 * g(0) = c4^2 R^4 w^2 == 1`  (auto-trivial),
   - `eps1(s_A) = A(0) = R^4 == 1`, `eps1(s_B) = B(0) lc(B) = w^2 c4^2 == 1`:
     **every rational 2-torsion class has trivial first component**, so all
     four order-8 classes `P_R + kappa, kappa in J[2](Q)` have the SAME
     first-component condition
     `eps1(P_R) = -c4 R in Q^2  <=>  -2R(R^2-1)(w^2-1) in Q^2` — the known
     conic.  (In fact `kappa in 2*J(Q)_tors` always — `J[2] = {0, 2e1,
     4P_R, 2e1+4P_R}` — so `div(P_R+kappa) <=> div(P_R)` outright,
     explaining why div_PR and div_PRHx have never disagreed.)
   - The only other coset of order-8 elements mod 2J(Q) is `P_R + h` with
     `2h = H_x` (h an order-4 class; h is forced rational if an order-16
     point exists via this coset, since h = 2x - P_R).  The halving identity
     for H_x (`g - x l1^2 = c4 q^2` with u_h = q) forces `q(0)^2 = R^4 w^2`,
     i.e. `q(0) = +-R^2 w` (matching the two tangent branches `V = +-R^2 w`
     of the membership machinery), so `eps1(P_R+h) = -c4 R q(0) =
     -+ c4 R w`.  Necessary condition for coset 2:
     `c4 R w in +-Q^2  <=>  2R(R^2-1)(w^2-1) w in +-Q^2`.

   **Disjunction theorem (first-component necessary condition): any (R,w)
   with exponent-16 torsion satisfies `-2R(R^2-1)(w^2-1) in Q^2` (coset 1,
   rationally parameterized by `w=(m^2+K)/(m^2-K), K=-2R(R^2-1)`) OR
   `2R(R^2-1)(w^2-1)w in +-Q^2` (coset 2, an elliptic condition
   `Y^2 = +-2R(R^2-1) w(w^2-1)` per fiber).**  Both flags are false on every
   known in-family point ((-8,6), (-16/11,14/11), (513/560,-663/700),
   (1/11,-1843/5093), (2,4), (5,25)) — consistent with all being exponent 8.

5. **Complete coset-1 sweep (pass 2).**  The conic parameterization makes
   the coset-1 stratum a rational curve per fiber; `scanRm` (octic
   `y^2 = q8_R(m)`) run on ALL 6086 fibers at m-height 2000 (6 s total!)
   reaches effective w-heights ~1e6*|K| and found exactly 7 points:
   the 2 known gold parasites, 1 known in-family point (-5/8,19/46), and
   4 new — all killed by p<=31; two of the new ones are **genuine new
   in-family points at record heights**:

   ```text
   (R,w) = (-32/65, 1816/4489)  h=4489  torsion [4,8]  div false  kill 17
   (R,w) = (-72/85, 5886/7459)  h=7459  torsion [4,8]  div false  kill 19
   ```

   Both **certified geometrically simple at multiple primes**
   (R=-32/65: p=53,83; R=-72/85: p=31,43,53; charpoly irreducible + 12th
   power transform irreducible; `data/claude_prod_416_new_points.txt` has
   integral models).  New family height record 7459 (old: 5093).  Banked.

6. **w=R^2 is an in-family section** (Magma-confirmed [4,8] at R=2,5):
   it appears as a torsion point on the `C^-` elliptic model of every
   fiber, which is why naive "no family points" fiber closure never fires;
   closure must (and now does) go through kill/exact bookkeeping of the
   torsion images.  Also: the old val_run "para" expectations for (2,4),
   (5,25) were wrong — `infam` agrees with Magma on all cross-checks.

7. **Coset-structure verification and coset-2 exclusions (2108 fibers).**
   Magma check on 3 fibers: `4*P_R = s_A` (never `H_x`), confirming the
   quarter-twist coset is the `H_x`-half exactly as in the kill-table
   criterion and in the (D) derivation.  The coset-2 stratum of fiber
   `R=a/b` is the elliptic curve `E2_R: Y^2 = c w(w^2-1)`,
   `c = 2ab(a^2-b^2)` (Weierstrass `X^3-c^2X`, congruent-number type;
   the `+-` branches fold via `w -> -w`).  `ellrank` on all 6086 fibers
   (`cos2.gp`, 72 s with caching on `core(c)`): **2108 fibers have r2=0 and
   torsion (Z/2)^2, whose w-images {0,+-1} are boundary — the coset-2
   stratum is empty there** (`data/claude_prod_416_coset2_closed.txt`).
   On these 2108 fibers any [4,16] must satisfy the coset-1 conic, whose
   complete octic sweep (result 5) is empty to effective height ~1e6*|K|.
   Octic-quotient Jacobian ranks on these fibers (`cos1jac.gp`,
   `ellfromeqn`+`ellrank`): 210 have both quotients r2=0 — all 210 already
   inside the 627 fully-closed set (the two rank conditions correlate), so
   no additional full closures from this angle yet; positive-rank
   coset-1 fibers need MW-lattice enumeration (next steps).

   Amusing structure: `E2_R` is the classical congruent-number curve with
   parameter `c = 2ab(a^2-b^2) = 4*Area` of the Pythagorean triangle
   generated by `(a,b)`; a fiber is coset-2-open iff `core(c)` is a
   congruent number.  Rank stats over the 6086 fibers: C+ r2 dist
   0:2073 1:2951 2:1006 3:56; C- r2 dist 0:1792 1:2822 2:1238 3:218 4:16;
   both-0: 627.  Coset-2-closed and fully-closed overlap in exactly the
   210 octic-rank-0 fibers; union coverage: 627 fully closed + 1898
   coset-2-only = 2525/6086 fibers (41.5%) with a rigorous stratum
   exclusion beyond the sweeps.

## Interpretation

The [4,16] obstruction on this chart is now structural on three
independent levels: (i) heights — nothing to w-height 1e5 for all 12172
R of height <=100 (1e6 partially, see resume); (ii) fiberwise theorems —
1254 fibers closed outright; (iii) descent — the entire exponent-16 locus
is confined to two thin conic/elliptic strata, and the rationally-sweepable
one (coset 1) is empty to effective height ~1e9 on every fiber.  Either
[4,16] lives at genuinely large height on a positive-rank fiber inside the
coset-2 stratum, or it does not exist on this chart and the proof route is
closing the coset-2 stratum fiberwise (see next steps).

## Resume state / next steps

- **Tier-2 (w-height 1e6) sweep**: jobs j1 (b=2..60) and j2 (b=61..84)
  running at session end (logs `tier2prog_j1.log`/`j2`, outputs
  `tier2_b*.txt` in `t416/prod/`); j3 (b=85..100) stopped for load-average
  compliance.  Resume:
  `cd <scratchpad>/t416/prod && RANGES="85:93 94:100" TAG=j3 ./sweep1e6.sh`
  then re-run `postproc.gp` with the tier2 chunk list (edit `chunks`) and
  exact-test survivors as above.  Expect the same Lplus parasites plus
  possibly new genuine points to bank.
- **Coset-2 closure (the sharp remaining weapon)**: per fiber the coset-2
  stratum is `Y^2 = +-2R(R^2-1) w(w^2-1)` (elliptic, full rational
  2-torsion) intersected with C_R^{+-}; a fiber is [4,16]-closed for
  coset 2 whenever this elliptic curve has rank 0 (torsion w-images are 0,
  +-1: all boundary!) OR H_x not in 2J(Q) (computable by the same descent).
  Wire `ellrank` on these cubics into `rank0.gp` — every fiber where either
  E^+ or the C-branch condition dies closes completely; combined with the
  existing 627 this could plausibly close most of the chart.
- Kill tables can be pushed to p<=499 with the same `killsets.m` (~4h);
  marginal (no kill above 73 was ever needed), lower priority than coset-2.
- The two new banked points and the height-record march suggest the family's
  positive-rank fibers keep producing [4,8] points indefinitely; the [4,16]
  wall is real.

## Addendum (2026-07-18 evening): tier-2 sweep COMPLETE — chart negative extended to w-height 1e6
All ranges redone cleanly after the aborted first launch (partial chunk files discarded —
no completion markers existed). Final: 6,086 fibers x 2 quartics to w-height 1e6, 30,571
cover candidates, 42-prime kill tables leave 14, and all 14 are Lplus boundary parasites
(w=(3R+1)/(R+3)) with singular models (poldisc(f)=0, verified 14/14). Zero smooth in-family
points in (1e5, 1e6]; zero exponent-16. Combined statement: no [4,16] or [2,4,16] on the
M_1(8,4) chart for any R of height <=100 and w of height <=1e6. Summary + survivors:
data/claude_prod_416_tier2_summary.txt. Next escalations per the main note: coset-2 closure
program on the 3,978 open fibers; R-height 200-300; Route 3 (deformation at the split
[4,16] moduli point).
