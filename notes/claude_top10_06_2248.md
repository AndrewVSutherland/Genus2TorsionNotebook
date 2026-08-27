# Top-10 target 6: `(2,2,4,8)` — strategy and first test run

*(Claude, 2026-07-17.  Scratchpad: `/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/` — scripts `validate_2248.py`, `calib_29.m`/`calib_29b.m`, `twisted2228.c`, `sweep_twisted.py`, `genus3_search.py`, `rk_e1e2.m`; outputs `tw200.txt`, `tw1000.txt`, `tw3000.txt`.  No repo file modified.)*

## 0. Dossier summary (what is known going in)

Target `(2,2,4,8)`, order 128: would tie the all-time genus-2 record (split HPL,
order 128) and beat Elkies' geometrically-simple record of 80.  Known points are
split HPL models only (`data/ten2248models_abcd.txt`, heights 10^19+).  Criterion
(root `main.tex` Lemma) on `y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)`: torsion
contains `(2,2,4,8)` iff BOTH the `(2,2,4,4)` conditions
(`(C-A)(C-B), (C-A)(D-A), (C-B)(D-B), (D-A)(D-B)` squares, `A=a^2` etc.) AND a
`(2,2,2,8)`-type condition set hold.  Exhaustive negatives
(`claude_tier1_item1_2248_descent.md`): all 30,387 primitive `(2,2,4,4)` tuples
`d <= 65535`, all 3 components x 16 twists — 0 hits; K3 route killed.  Live
routes: (i) the TWISTED `(2,2,2,8)` family (component I, `e(A,B)`-twist) —
square set `{abcd, a(a+b)(c-a)(d-a), b(a+b)(c-b)(d-b), c(c+d)(c-a)(c-b),
d(c+d)(d-a)(d-b)}` — 52 tuples `d <= 200`, never intersected with the
`(2,2,4,4)` side; (ii) the HPL `(Z/2)^4`-cover threefold `M(2,2,4,8)` (HPL is
NOT on the Faltings-controlled surface `S'`); (iii) the genus-3 split-locus
curve `{y^2=(u-3)(u+1)(u^2-6u+1), z^2=-(u-1)(u^2-6u+1)}` of `S'`, rank-1
quotients, Faltings-finite, undecided.

## 1. Strategy: ranked routes

**Route A (best odds per CPU-hour): enumerate the twisted `(2,2,2,8)` family and
intersect with the `(2,2,4,4)` side via 2-descent.**  Mathematical content: a
twisted-family tuple has a rational `G` with `2G = D0 + e(A,B)`, hence `4G = T_0`
— an order-8 point.  Torsion then contains `(2,2,4,8)` iff some OTHER 2-torsion
class `T' != T_0` is 2-divisible, i.e. `delta(T')` trivial (15 classes to test;
the Weierstrass-pair delta convention of `claude_twist_sweep.py`, Magma-validated
with 0/600 mismatches).  This is the **second component** of the `(2,2,4,8)`
locus over the twisted family and had never been enumerated at any height —
the 65535-exhaustive negative covers only tuples that already satisfy the
`(2,2,4,4)` square set, a *different* slice.  Enumeration is cheap: `abcd`
square forces `d = sqfree(abc)*m^2`, so the search is `O(B^3)` gcds, not
`O(B^4)`; production target `B = 10^5`-`10^6` with `__int128` square tests.
Bonus: the same list feeds twisted analogues of components II/III (halving of
`T_AB`/`T_CD` on twisted-family curves — condition sets not yet derived; the
`F_i(rho,sigma,tau)` of `paper/m2248_equations_gpt55.tex` give component II
after transport by the `e(A,B)` twist).

**Route B (structural, infinite-family prize): the HPL threefold
`M(2,2,4,8) -> A(2,2,4,4)`.**  HPL lies on the `(Z/2)^4`-cover threefold, not on
`S'`, so no fibration-finiteness applies.  Concrete step: fiber the threefold
over the `(rho, sigma)`-plane instead of `P^1_t` and hunt for a *horizontal*
rational curve through the HPL point (the HPL models are multiples `nP` on a
positive-rank elliptic surface — find that surface inside `M(2,2,4,8)` in the
`(rho,sigma,tau)` chart, `rho0 = 58466134224/53109477625`, and search for other
low-height sections).  Expensive, but it is where an infinite family would live.

**Route C (well-posed finiteness problem): the genus-3 split-locus curve.**
Both genus-1 quotients have rank exactly 1 (reconfirmed below), so classical
descent doesn't finish; the tools are elliptic-Chabauty over `Q(sqrt2)`
(the curve maps to rank-1 elliptic curves over the field where `u^2-6u+1`
splits) or two-cover descent.  Any nondegenerate point gives a NEW `(2,2,4,8)`
curve on `S'` itself; proving emptiness finishes the impossibility program on
`S'`.  Low prior of a point (search below: none to height 2000), high value
either way.

## 2. Validation (all dossier facts reproduced)

* `validate_2248.py` (pure integer): smallest HPL tuple
  `a=16336390342285800000, ...` satisfies all four `(2,2,4,4)` conditions AND
  all four `(2,2,2,8)` conditions — the Lemma's criterion checks out on the
  known record curve.  `(29,121,125,145)`, `(2,3,12,18)`, `(1,2,4,50)`,
  `(34,41,68,82)` satisfy all 5 twisted conditions and FAIL all 4 untwisted
  ones (`n_a(2,3,12,18) = 2800`, as stated).
* Magma (`calib_29b.m`): `(29,121,125,145)` has torsion invariants exactly
  `[2,2,2,8]` and passes the simplicity certificate at `p = 61`
  (charpoly `x^4+8x^3+70x^2+488x+3721` irreducible, 12th-power transform
  irreducible of degree 4).  Note `p = 13, ...` give reducible charpolys —
  certificate needs the prime scan.
* Enumerator cross-check: at `B = 200` the new enumerator returns **exactly the
  52 tuples** of the dossier (25 with `gcd(a,b,c,d)=1`), including all three
  quoted examples, in 0.03 s.
* Descent-code cross-check: `delta(T_0)` is trivial for **52/52** tuples at
  `B=200` (forced by `T_0 = 4G` on the family) — enumerator and delta code
  agree perfectly.
* `rk_e1e2.m`: `E1: y^2=(u-3)(u+1)(u^2-6u+1)` and `E2: z^2=-(u-1)(u^2-6u+1)`
  both have `RankBounds = (1,1)`, torsion `Z/2` — as in
  `claude_tier1_item2_genus_drop.md`.

## 3. Test run (route A first step) — commands and results

* `twisted2228.c` (single-threaded, sqfree-forced `d`, mod-64 + exact-sqrt
  square tests; compile `gcc -O2 -march=native`):
  * `B=200`: 52 tuples, 0.03 s.  `B=1000`: 489 tuples, 4.3 s.
  * `B=3000`: **1947 tuples (439 gcd-primitive)**, 129.7 s single-threaded
    nice-15 on the loaded box (list: scratchpad `tw3000.txt`; e.g. largest
    primitive `[648,1875,2000,3000]`).  The family keeps growing roughly
    linearly in `B` (52 / 489 / 1947 at 200 / 1000 / 3000).
* `sweep_twisted.py <list>` (imports `delta_of_U` from
  `code/claude_twist_sweep.py`): for each tuple tests all 15 nonzero 2-torsion
  classes for 2-divisibility.  Any trivial `delta(T')` with `T' != T_0` is a
  `(2,2,4,8)` hit.
  * Result on all 1947 tuples (~40 s): `delta(T_0)` trivial **1947/1947**
    (the forced consistency check — enumerator and descent code agree on the
    whole list); best-square-component histogram over the other 14 classes
    `{1: 1519, 2: 212, 3: 216}`; **FULL HITS: 0**.  4/5 never occurs
    (product-relation parity, same as component I), and 3/5 is common (11%),
    mirroring the structure of the 65535-exhaustive first-component sweep.
* Route (iii) bonus, `genus3_search.py`: exhaustive `u = p/q`,
  `|p|, q <= 2000`, on the genus-3 curve: **only the 3 degenerate points
  `u = -1, 1, 3`** (roots of `A*B`); no nondegenerate point to height 2000.

## 4. Verdict and next steps

**Verdict: pipeline validated end-to-end; second component NEGATIVE to
`d <= 3000`; realizability still open but the odds shift slightly down.**  The
second `(2,2,4,8)` component behaves statistically like the first (best 3/5,
never 4/5): the descent conditions look independent enough that a hit at
height `H` has probability roughly `H^{-1}`-ish per candidate, so a
positive-density family of twisted tuples (linear growth observed) gives a
log-divergent expected hit count — the search is worth pushing 2-3 orders of
magnitude higher before drawing structural conclusions, and the machinery to
do so now exists and is calibrated.  Meanwhile route C stays the only finite,
decidable sub-problem (empty to height 2000), and route B the only credible
source of an infinite family.

1. Production run of `twisted2228.c` (parallel, `__int128`) to `B = 10^5`,
   then `sweep_twisted.py` (ported to C for the mod-square tests) over the
   full list — the honest first exhaustive statement for the second component.
2. Derive the twisted analogues of components II/III (halving `T_AB`, `T_CD`
   over the twisted family) from the `F_i` of `m2248_equations_gpt55.tex` and
   add them to the sweep — 3x the coverage per tuple.
3. Elliptic-Chabauty over `Q(sqrt(2))` on the genus-3 split-locus curve
   (rank-1 quotients make it a textbook case) — would close `S'` entirely.
4. Route B: locate the positive-rank elliptic surface through the HPL point
   inside `M(2,2,4,8)` in the `(rho,sigma,tau)` chart and search for
   independent sections of lower height.
