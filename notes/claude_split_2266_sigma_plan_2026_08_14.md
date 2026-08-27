# Plan: [2,2,6,6] via the five non-identity sigma-surfaces (next session)

*Written 2026-08-14, for execution in a fresh session.  Successor to
`notes/claude_split_2266_1212_plan_2026_08_13.md`; the 2026-08-13 execution
session (`notes/claude_split_2266_1212_session_2026_08_13.md`, rendered
report `reports/split-144-hunt/`) built the machinery, closed the self-gluing
and (N,N)-congruence routes, and left exactly one structured attack surface
untouched: the fibrations of the five non-identity sigma-surfaces.*

## 0. Background established 2026-08-13 (do NOT re-derive; re-verify only)

- **The complete criterion.**  A (2,2)-gluing of E26(t) x E26(u) with
  root-matching sigma in S3 has J(Q)_tors >= [2,2,6,6] (order 144, odd part
  [3,3] automatic, full rational J[2]) **iff** the three unordered-pair
  class-matching conditions hold:

      a_o(i,k)(t) * a_o(sigma(i),sigma(k))(u) = square,  {i,k} in {12,13,23},

  with ordered classes a(1,2) = (t+3)(t-5), a(1,3) = 2(t-3),
  a(2,3) = -(t-1)(t-9), and a(k,i) = -a(i,k).  This is anchored five ways in
  `product/logs/derive2266.log` (iota(T1) = HLP printed value; sigma=id =
  HLP (6); E28xE28 = HLP (5); trivial iota row of E28; the [2,2,24]-witness).
  **Step 0 of the new session: rerun `lane_2266_derive.m` (2 min) and use its
  MECHANICAL sigma-system table — never hand-transcribed conditions.**
- The six sigma-systems, for reference (variable-wise: LA(t)*LB(u) square,
  from the derive log; re-generate before use):

      id    : (t+3)(t-5)~(u+3)(u-5) | 2(t-3)~2(u-3)       | -(t-1)(t-9)~-(u-1)(u-9)
      (12)  : (t+3)(t-5)~-(u+3)(u-5)| 2(t-3)~-(u-1)(u-9)  | -(t-1)(t-9)~2(u-3)
      (13)  : (t+3)(t-5)~(u-1)(u-9) | 2(t-3)~-2(u-3)      | -(t-1)(t-9)~-(u+3)(u-5)
      (23)  : (t+3)(t-5)~2(u-3)     | 2(t-3)~(u+3)(u-5)   | -(t-1)(t-9)~(u-1)(u-9)
      (123) : (t+3)(t-5)~-(u-1)(u-9)| 2(t-3)~-(u+3)(u-5)  | -(t-1)(t-9)~-2(u-3)
      (132) : (t+3)(t-5)~-2(u-3)    | 2(t-3)~(u-1)(u-9)   | -(t-1)(t-9)~(u+3)(u-5)

  ("~" = product must be a square.)  Swap symmetry: (t,u) -> (u,t) maps the
  sigma-system to the sigma^{-1}-system, so (123) and (132) are ONE surface;
  transpositions are self-paired.  The E26 deck transformations may induce
  further identifications between sigma-surfaces — dedupe empirically via
  glued-model G2-invariants, as always.
- **What is already closed (do not repeat):** sigma=id cross pairs empty to
  height 200 (join) and its (6)-surface fibration dry (MW and
  section-arithmetic versions, plus all 21 symbolic section-multiple
  condition curves — genus 7+, no nondegenerate points); self-gluings
  impossible (rank 0); [2,2,2,24] via [2,6]x[2,8] impossible; the whole
  (N,N)-congruence/isogeny gluing route structurally obstructed; LMFDB
  congruent-pair lists exhausted.  See session note sections 2, 4.5-4.7, 5.
- **What is known positive:** `lane_misc2.log` census — each non-id
  sigma-surface has 231-388 rational points satisfying 2 of its 3 conditions
  at height 40 (every drop-choice).  These surfaces have never been fibered.
- Funnel infrastructure unchanged (`split_lab.m`; [2,2,6,6] NOT in KNOWN, so
  a hit prints `HIT ... invs=[ 2, 2, 6, 6 ]`).  GlueFunnel + OddInvs [3,3].
  House rules: <= 3 concurrent Magma jobs, `-b` buffers (progress
  side-channel via `System("echo ... >> log")` every N iterations — this
  saved the session twice), excluded parameter values {3,-3,1,5,9} both
  sides, j-equal pairs are deck pairs (skip).

## 1. Primary attack: join-driven fiber analysis of the four new surfaces

The sigma=id fibration had a universal section (HLP's (y,z)=(-1,1)); the
non-id surfaces do not (no diagonal: t=u fails their conditions), so their
fibers are genus-1 TORSORS.  The census shows many fibers are nonetheless
soluble — the join finds their points, and each populated fiber becomes an
elliptic curve whose Mordell-Weil group reaches t-heights the join cannot.

1. **Build `lane_2266_sigma.m`** (variant of `lane_2266.m`): for each
   sigma in {(12),(13),(23),(123)} and each drop-choice d in {1,2,3}, run
   the 2-of-3 hash-join at H = 150 (kernels via SFrat as before; join is
   O(N log N), minutes).  For every surface point (t,u) found, IMMEDIATELY
   test the dropped third condition (one IsSquare) — any pass is a full
   [2,2,6,6] candidate: funnel on the spot.  Expected volume: censuses
   ~300/sigma at H=40 suggest 2-5k points per (sigma, drop) at H=150.
   Output: per-sigma point lists `product/data/sigma_<s>_pts.txt`
   (t, u, drop) — the raw material for everything below.
2. **Fiber organization.**  For each sigma pick the drop-choice whose two
   REMAINING conditions contain a class linear in t (2(t-3) appears in
   every system on one side or the other): then for fixed u0 the fiber is
   `t = 3 + w^2/(2*c_B(u0))` substituted into the other condition — a single
   genus-1 quartic in w.  Group the step-1 points by u0.  Statistics to
   print: #populated fibers, #fibers with >= 2 points (rank >= 1 hints),
   heights.
3. **Per-fiber MW attack** (the discipline that failed and the fix, from
   2026-08-13: unbounded MordellWeilGroup per fiber stalled; GRH bounds
   alone were not enough):
   - fiber has a point P0 from step 1 -> EllipticCurve(C, P0);
   - `RankBounds(E)` first; if bounds are NOT equal, log `MWSKIP` and fall
     back to `Points(C : Bound := 10^5)` only;
   - if rank 0: test the torsion points only (finite);
   - if conclusive rank >= 1: MordellWeilGroup, enumerate the coefficient
     box |n_i| <= 10 + torsion, map back to t, test the third condition at
     every point (exact IsSquare at any height), funnel survivors;
   - progress side-channel every 10 fibers; hard per-sigma funnel cap.
   Budget ~2-3h across the four surfaces; shard per-sigma as separate
   background jobs (3 at a time).
4. **Section detection.**  On each per-sigma point list, look for structure:
   (a) fibers of u0 with points at several heights (rank >= 1 confirmed);
   (b) parametric families: fit t as a rational function of u of degree
   <= 2 through subsets (exact fitting over Q on 5+ points; a hit is a
   SECTION of the fibration and upgrades the surface to positive rank over
   Q(u) — then run the full section-multiple machinery as in
   `lane_2266_sections.m`, which is written and debugged).
5. Any HIT: standard checklist (fresh-session exact TorsionSubgroup on the
   integral model, L-poly split at good p < 200 against the two E26 factor
   traces, Richelot SetCart factor identification) and the recording
   pipeline of the 08-13 plan section 3 (witnesses JSON -> SRC -> KNOWN ->
   table regen + drift-guard verifier; the order-record line changes:
   144 beats 128).

## 2. Secondary lane: genuine 5/7-congruent [2,6]-pairs at family scale

The only [2,2,6,6]-route outside (2,2)-gluings that survives the 08-13
obstructions is a NON-isogenous (5,1)- or (7,3)-congruent pair of
[2,6]-curves (product injection; both moduli allow it — no small-conductor
instances exist, but the family reaches conductors ~10^15).

1. Extend `lane_cong_sieve.m` to the E26 family alone at H = 150 (~25k
   curves; the optimized pair loop does ~310M pairs in ~15 min).  Sieve
   mod 5 AND mod 7 (trace vectors, MINCMP 14).  Flag NON-isogenous
   candidates only (`IsIsogenous` on survivors).
2. Any genuine pair -> `GlueScan` (analytic_glue.m) at prec 200; a
   Q-rational anti-isometry now REQUIRES the strict rationality test
   (err < 10^-(2*height(q)+15), heights <= (prec-30)/2 — the 08-13 session's
   hard-won lesson; the weak test accepts every real number).  If verified
   rational and heights <= ~60 digits: Mestre via the I2=1-normalized
   ReduceIC tuple; expect Mestre to be the bottleneck (documented) — cap
   attempts, record invariants regardless (a certified Q-surface with
   [2,2,6,6]-injection is reportable even if the curve model lags).

## 3. If everything is dry: decide thin-vs-obstructed, and close cleanly

1. **Local analysis** (`local-obstructions` skill): for each sigma, test the
   full 3-condition variety for Q_p-points, p <= 50, and R-points (the
   (Z/2)^2-cover of the (t,u)-plane; solubility = existence of local
   (t,u,z1,z2,z3)).  Everywhere-locally-soluble + globally dry = genuine
   thinness (report as such); a local obstruction on EVERY sigma-component
   would be a strong structural theorem — check carefully before claiming.
2. **Geometry of the surfaces** (optional, half-day cap): Kodaira dimension
   of one sigma-surface S~ (triple (Z/2)^2-cover of P^1 x P^1 branched over
   three (2,2)-curves): general type would justify "finitely many points
   expected on each" (Bombieri-Lang heuristic) for the writeup; K3/rational
   would argue for pushing the fiber search harder.
3. **Cleanup theorem** (30 min): finish the twisted-diagonal impossibility
   formally — the two genus-3 curves have RankBound(J) = 0; enumerate
   J(Q)_tors explicitly (reduction injectivity at two good primes) and
   verify the known degenerate points exhaust C(Q).  Turns session-note
   item 2.4 into an unconditional statement.
4. Update `product/split_torsion_table.md` commentary and the session-note
   frontier accordingly; the [8,8] infinite family (t = a^2/(a^2+b^2))
   should also be offered to the paper as a remark regardless of the
   [2,2,6,6] outcome.

## 4. Logistics

- New lanes from `product/code/`: `lane_2266_sigma.m` (join + fiber file),
  `lane_2266_sigma_fib.m` (per-fiber MW, sharded per sigma),
  `lane_cong_sieve.m` (H:=150, E26-only switch).  Reuse GlueFunnel verbatim.
- Validation gates before trusting anything new: (i) rerun
  `lane_2266_derive.m` — all five anchors must PASS; (ii) the sigma=id
  2-of-3 join must reproduce the known (6)-surface point counts (positive
  control for the join code); (iii) one deliberately-wrong sigma-label run
  must produce ~zero matches (negative control).
- ETAs: joins minutes; fiber MW hours (shard, cap, side-channel); sieve
  ~20 min.  Keep <= 3 Magma jobs; MemGB 4-6.
- Known dead ends (session note section 5 + 4.5-4.7): do not re-run
  self-gluings, diagonals, (N,N)-congruence/isogeny gluings, LMFDB lists,
  sigma=id fibrations, or blind volume.
