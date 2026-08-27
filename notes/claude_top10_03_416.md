# [4,16] per-R elliptic solve on M_1(8,4): the O(H^4) wall is down

*(Claude agent, 2026-07-17.  Target #3 of the top-10: torsion `[4,16]`
(order 64) on a geometrically simple genus-2 Jacobian over Q.  Session
scratchpad `t416/`: `validate416_known.m`, `perR.gp`, `scan_all.gp`,
`cert_new.m`; outputs `validate416_known.out`, `pass1_all.txt`,
`survivors_new.txt`, `gold.txt`, `exact_out.txt`, `cert_new.out`.
No repo file modified; total ~8 CPU-min, single-threaded, nice 15.)*

## 0. Dossier context

Split `[4,16]` is known (Sutherland split-not-product); the paper's old
candidate `y^2 = 2025x^5+11484x^4+9846x^3+11484x^2+2025x` is palindromic
hence split (`NotesAndTodo.tex` (8,4) section, at `(s,u,t)=(-3,0,4)`,
i.e. the `u=0` symmetric slice; on the `(R,w)` chart this is the
boundary `R=-1`).  Route on the `M_1(8,4)` `[4,8]` family (chart
`(R,w)`, `n=1`): halve the order-8 class `P_R = [(-R,Y_R)] - inf`.
State of play before this session (`notes/claude_next_416_route_revived.md`,
`notes/claude_next_416b_search.md`): no local obstruction (218 ALIVE
7-adic reps on `R+-w = 0 mod 7`), corrected all-twist kill tables
p=11..43 on file (`data/m18_m14_416_killsets_corrected.txt`, 4368
residues), rigorous negative to `(R,w)`-height 800 via O(H^4) blind pair
enumeration (6e11 pairs -> 96 sound covers -> 52 in-family points, all
`[4,8]`/`[2,4,8]`), near-misses at height 700: `(513/560, +-663/700)`.
The spec'd-but-unbuilt unlock was the per-R solve.  **This session built
it, validated it end to end, and it works.**

## 1. Headline results

1. **Validation clean.**  The known one-split point `R=-16/11, w=14/11`
   reproduces exactly: `t=42/55`, `disc(A)` square / `disc(B)` not,
   `Order(P_R)=8`, `J(Q)_tors = [2,4,8]`, charpoly at 47
   `T^4-4T^3+30T^2-188T+2209` irreducible with irreducible 12th-power
   transform, digit-for-digit the `g12` in `notes/m18_m14_halving.md`.

2. **The per-R fibration is explicit and genus 1.**  For fixed `R=a/b`,
   family membership (the discriminant-cover condition) says `w` is a
   rational point on one of two genus-1 quartics with integer
   coefficients (`b^6` denominators and `4(a-b)^2` are squares, dropped):

   ```text
   C_R^+ : y^2 = -(a+b) (b^2 w^2 - a^2) (w+1) ((a+3b)w - (3a+b))
   C_R^- : y^2 =  (a+b) (w-1) ((a+3b)w + (3a+b))
                        (b^2(a^2+2ab-b^2) w^2 + a^2(a^2-2ab-b^2))
   ```

   Symbolic identity vs the notes' `Delta_plus/Delta_minus` machine-checked
   (`symcheck()` in `perR.gp`).  All 9 known in-family `(R,w)` points
   (5 base + 4 near-miss signs) are re-found on the expected curve
   (`knowncheck()`; base points and near-misses on `C^+` except the
   `-16/11, -11/16` pair on `C^-`).  Solving one fiber to w-height 1e5
   with `hyperellratpoints` costs **~0.1 s** (PARI 2.18).

3. **New aux-free necessary condition for halving `P_R` (the conic).**
   On the chart, `c4 := lc(f) = 2(R^2-1)/(w^2-1)` (machine-checked).
   The constant coefficient of the (416) identity
   `f - ell^2 = c4 (x+R) q^2` at `x=0` forces (with `q = x^2+ax+b_q`)
   `e^2 = -c4 R b_q^2`, and `b_q = 0` is impossible off-boundary
   (degree count forces `A(0)B(0)=0`).  Hence

   ```text
   P_R in 2J(Q)  ==>  -c4 R = -2R(R^2-1)/(w^2-1) is a rational square
                 <==>  y^2 = K(w^2-1) has a point with that w,
                       K := -2R(R^2-1).
   ```

   The conic has the rational point `(w,y)=(1,0)`, so the whole
   div(P_R)-compatible stratum is rationally parameterized:
   `w = (m^2+K)/(m^2-K)`.  Substituting into `C_R^+-` gives an even
   octic `y^2 = q8_R(m)` per fiber — a second, sharper per-R solve
   targeting the untwisted class specifically (`scanRm` in `perR.gp`).
   Consistency: conic false at `(-8,6)`, `(-16/11,14/11)` (both known
   non-divisible), and every conic-true candidate found below is
   kill-filtered and indeed exact-tests to `divPR=false`.

4. **Test scan (186 R values, w-height 1e5, 30 s):** all reduced
   `R=a/b` of height <= 12 plus the known/near-miss specials.  Found
   **1280** distinct cover-level `(R,|w|)` candidates: 760 in the old
   territory (`max height <= 800`), **520 NEW** beyond the old wall.
   Sound kill filter (p=11..43 tables) leaves **6 unkilled**:
   the known pair `(-11,4)/(-1/11,1/4)`, the dossier's near-miss pair
   `(513/560,-663/700)` (exact cross-check of the old funnel), and a
   **new pair at height 5093**: `(R,w) = (1/11, -1843/5093)` and its
   `R -> 1/R, w -> 1/w` mirror.  Kill-prime histogram of the rest:
   13:358, 17:342, 19:162, 11:122, 29:108, 23:94, 31:42, 37:24, 41:14,
   43:8.

5. **The height-5093 candidate is a genuine NEW in-family point** —
   the first beyond the height-800 frontier (a blind-pair sweep to 5100
   would have cost ~4000x the height-800 run; the per-R solve found it
   in 30 s).  Exact stage (`code/m18_m14_416_exact.m`, unmodified):

   ```text
   INFAMILY R=1/11 w=-1843/5093 torsion=[4,8] divPR=false divPRHx=false
   ```

   Not `[4,16]` — but it extends the family census 6.4x in height and
   is **geometrically simple**: certificate at p=43,
   `chi = T^4+4T^3+22T^2+172T+1849` irreducible, 12th-power transform
   irreducible (`cert_new.out`).  (Banked: new highest-height simple
   `[4,8]` member.  Note its conic flag is 0, so `P_R` could not halve;
   only the twisted class was live, and it fails too.)

6. **Rank profile of the fibration** (PARI `ellrank`, all bounds tight):

   ```text
   R      : C+ rank  C- rank      R      : C+ rank  C- rank
   1/11   : 1        2            -6/5   : 1        2
   -8     : 0        0            1/3    : 1        0
   2      : 0        0            9/16   : 2        0
   3/2    : 0        0            5      : 0        0
   ```

   Both-rank-0 fibers (`R = -8, 2, 3/2, 5` here) are **completely and
   unconditionally closed**: the only rational `w` are torsion points
   (all boundary), so no `[4,16]` exists over those `R` for ANY `w` —
   a per-fiber nonexistence theorem, strictly stronger than any height
   bound.  Positive-rank fibers (e.g. `1/11`, ranks 1+2 — explaining
   its 47 cover points to height 1e5) are where all candidates live and
   where MW-lattice enumeration reaches heights no sieve can.

## 2. Strategy (ranked)

**Route 1 — production per-R funnel (validated here, primary).**
Enumerate `R` by height (<= 100-1000; use the `R -> 1/R` involution to
halve work).  Per fiber: `hyperellratpoints` on `C_R^+-` at w-height
1e5-1e6 (~0.1-1 s) -> sound kill filter -> exact stage
(`m18_m14_416_exact.m`).  Two force multipliers: (i) extend the kill
tables from p<=43 to p<=199 with the existing corrected all-twist
criterion (each new prime kills ~60-70%; expect ~(1/3)^k thinning —
the test's 6/1280 survival at 10 primes extrapolates to ~1 survivor
per ~10^4 fibers at 20 primes, so exact-stage load stays trivial);
(ii) for any fiber with an unkilled candidate, run `ellrank`: rank 0
closes the fiber unconditionally, rank >= 1 gives MW generators for
exhaustive lattice enumeration far beyond sieve heights.  Do NOT
hard-filter on the heuristic mod-7 alive discs (the near-miss pair has
`7 | den`); use them for prioritization only.

**Route 2 — class-targeted conic solves (the sharp weapon).**  (a) For
the untwisted class: per-R even octic `y^2 = q8_R(m)` from the conic
parameterization (built, `scanRm`); being even it covers two elliptic
quotients in `M = m^2`, so the `ellrank`-closure trick applies to a
genus-3 curve — rank-0 quotients kill the whole div(P_R) stratum of a
fiber.  (b) Derive the analogous aux-free constant-coefficient
condition for the twisted class `P_R + H_x` (halving normal form
`f - ell^2 = c4 * u2 * q^2` with `u2` the degree-2 Mumford u-polynomial
of `P_R + H_x`, computable once symbolically on the first-cover chart:
condition `c4 * u2(0) in Q^2` — a second conic).  The test showed the
kill-filter survivors have conic=0, i.e. they live entirely on the
twisted branch — (b) is exactly the filter that addresses them.

**Route 3 — deform off the split point in M_1(8,2^w).**  The `[4,16]`
locus in the 3-fold `(s,u,t)` chart is nonempty: it contains the split
point `(-3,0,4)` on the palindromic slice `u=0` (`NotesAndTodo.tex`).
Compute the local equations of the `[4,16]` locus at that point
(halving cover of the order-8 point plus the order-4 structure) and
trace its component transverse to `u=0`; any rational point with
`u != 0` off the palindromic symmetry is a simple candidate.  This
bypasses the `(R,w)` chart's lost strata (`R in {0,+-1}`, `b_q=0`,
`n=0`).  Riskier (the component may be contained in `u=0`), but it is
the only route anchored at an actual `[4,16]` moduli point.

## 3. Exact test commands

```text
cd <scratchpad>/t416
nice -n 15 magma -b validate416_known.m          # validation, ~35 s
nice -n 15 gp -q scan_all.gp                     # passes 1+2, ~31 s
#   pass1: 186 fibers x 2 quartics, hyperellratpoints H=1e5
#   pass2: conic-targeted octics, Hm=2000
nice -n 15 magma -b infile:=exact_in.txt \
    /home/claude/torsion_jac/code/m18_m14_416_exact.m   # ~3 min
nice -n 15 magma -b cert_new.m                   # simplicity cert, ~30 s
# ellrank probe: probe() in perR.gp, 16 curves, ~80 s
```

`exact_in.txt` = 2 kill-survivors + 4 conic-true (gold) points; result:
4 in-family, all `[4,8]`, all divisibility false, 0 hits — every layer
(conic, kill tables, x-T exact, torsion) mutually consistent.

## 4. Corrections / cautions

- The conic condition is necessary only for the UNTWISTED class; do not
  use it as a global filter (kill-survivors were conic=0 and still
  legitimately live via `P_R + H_x`).
- `hyperellratpoints` output includes `y=0` and boundary `w`; filter
  `w in {0, +-1, +-R}` and `y != 0` (singular strata `Lpm, Q, Quartic`
  are factors of the quartics, hence appear as `y=0`).
- `(R,w)` and `(R,-w)` give the same curve (`f` depends on `w^2` only);
  dedupe on `|w|`, and `TangentCandidates` tries both `V = +-R^2 w`.

## 5. Next steps

1. Extend kill tables to p <= 199 (Magma, group-theoretic all-twist
   criterion as in the audit; ~1 h single-core, reusable data file).
2. Production Route 1: all `R` of height <= 100 (~ 12k fibers after the
   involution), w-height 1e6 on unkilled fibers; exact stage on
   survivors.  ~1-2 h on a few cores, all rigorous.
3. Symbolic derivation of the twisted-class conic (`c4 * u2(0)` square)
   and its per-R octic; wire into pass 2.
4. `ellrank` bookkeeping: log every both-rank-0 fiber as an
   unconditional per-R nonexistence certificate (these accumulate into
   a statement of the form "no `[4,16]` with `R` in explicit infinite
   list").
5. Bank by-products: every new unkilled in-family point gets exact
   torsion + simplicity certificate (the height-5093 simple `[4,8]`
   example above is the first).
