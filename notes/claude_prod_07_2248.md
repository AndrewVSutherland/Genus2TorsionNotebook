# PROD 07: (2,2,4,8) — production run, split-locus closure, and the sign-reduction theorem

*(Claude, 2026-07-18.  Scratchpad: `/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/t2248prod/` — all scripts and raw outputs there; found tuples in `data/claude_prod_07_*`.)*

## Strategy recap (3 lines)

(A) Production-enumerate the twisted (2,2,2,8) family (component I, e(A,B)-twist) far past the prior d<=3000, and sweep every tuple for a second 2-divisible 2-torsion class (= (2,2,4,8) jackpot).
(B) Determine the rational points of the genus-3 split-locus curve `{y^2=(u-3)(u+1)(u^2-6u+1), z^2=-(u-1)(u^2-6u+1)}` of `S'` completely.
(C) Derive twisted analogues of components II/III and add to the sweep.

## Headline results

1. **Task B SOLVED, no Chabauty needed.**  The genus-3 curve has a *third*
   elliptic quotient, missed by the tier-1 analysis (which only used E1, E2 of
   rank 1): `w := yz/(u^2-6u+1)` satisfies `w^2 = -(u-3)(u+1)(u-1)` (exact
   identity, Magma-verified).  That curve `E3` (`y^2=(x-1)(x+1)(x+3)` after
   `x=-u`) has **RankBounds (0,0)** (unconditional 2-descent: 2-Selmer dim 2 =
   torsion dim) and torsion `(Z/2)^2`, so `E3(Q) = {O,(1,0),(-1,0),(-3,0)}`
   and every rational point of the genus-3 curve has `u in {3,-1,1,infinity}` —
   exactly the 8 known degenerate points (`(3,0,±4), (-1,0,±4), (1,±4,0)`, two
   points over `u=infinity`), all mapping to boundary `t in {0,infinity}`.
   **The split locus of `S'` is now completely closed**: family 1 empty (sign
   obstruction, tier-1), family 2 empty of nondegenerate points (rank-0
   descent, this session).  Script/output: `e3_close.m` / `e3_close.out`.

2. **Sign-reduction theorem (new; collapses the second-component sweep).**  On
   `y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)` with `0<a<b<c<d`, computing the
   Magma-validated delta vectors symbolically:
   * `delta(e(0,inf))` is identically trivial (the order-4 class `D0` is always rational);
   * the 8 classes `e(0,X), e(X,inf)` have `delta` component at `x=0` equal to
     `-X = -x^2` mod squares — identically negative, never trivial;
   * `delta(e(A,C)) = delta(e(B,D))` has component `(A-B)(C-B) < 0`,
     `delta(e(A,D)) = delta(e(B,C))` has component `(A-B)(D-B) < 0` — never trivial;
   * `delta(e(A,B)) = delta(e(C,D)) = [1, (C-A)(D-A), (C-B)(D-B), (A-C)(B-C), (A-D)(B-D)]`.
   Hence the ONLY possible second 2-divisible class is `e(A,B) ~ e(C,D)`, and
   the 15-class sweep reduces to **4 integer square conditions** — which are
   precisely the partition-`(2,2,4,4)` conditions of Lemma 2244(c).
   Numerically confirmed: 295 sampled tuples run through the full 15-class
   delta machinery, 0 mismatches with the theorem (`fastsweep.py`).

   **Corollary (subsumption).**  A (2,2,4,8) hit on ANY tuple of this model
   must satisfy the partition-(2,2,4,4) squares, i.e. be a `tor2244` tuple, and
   the order-8 chain is then one of the 3 components x 16 twists tested in the
   tier-1 exhaustive audit (`claude_tier1_item1_2248_descent.md`: 30,387
   primitives `d<=65535`, 0 hits).  Together with delta-injectivity this gives:
   **no curve `y^2 = x prod(x+x_i^2)` with `0<a<b<c<d<=65535` has torsion
   containing (2,2,4,8)** — through ANY route, twisted families included.  The
   twisted-family sweep below 65535 was therefore mathematically redundant (and
   indeed returned 0); new territory starts at `d > 65535`.

3. **Task A production enumeration** (`twisted2248_prod.c`, forced-squarefree-kernel
   `O(~B^2 polylog)` algorithm, `__int128` exact square tests, pthreads —
   validated: exact match with prior lists at B=200/1000/3000, fresh-run match
   at B=10000, and the K3 section tuple `[98,144,147,294]` is found):
   * `B=30000`: **29,426 tuples** (275 s, 3 threads), file `v30000.txt`.
     Fast 4-condition sweep: **0 hits**; #satisfied-of-4 histogram
     `{0: 28565, 1: 409, 2: 452}` (nothing even reaches 3/4).
   * `B=100000`: 3 sequential b-chunks `[2,63300),[63300,84500),[84500,100000)`,
     3 threads (`v100k_c{1,2,3}.txt`) — results in the RESULTS section below /
     resume state if incomplete.
   * Fresh-tuple end-to-end check: `[6,68,1881,21318]` and `[80,81,150,28830]`
     have torsion exactly `[2,2,2,8]` in Magma (`fresh_check.m`) — the deep
     enumeration produces genuine order-8 curves.

4. **Task C resolved structurally.**
   * Components II/III (`4G = T_AB / T_CD`) require `T_AB` (resp `T_CD`)
     2-divisible, hence the partition squares, hence a `tor2244` tuple — they
     have **no tuples outside the tor2244 list** and are fully covered to
     65535 by tier-1.  Nothing new to enumerate below that height.
   * The 16 twists of component I collapse by the sign analysis to exactly
     **two sign-viable families**: the e(A,B)-twist (our enumerated family) and
     the identity twist `{abcd, a(a+b)(a+c)(a+d), ..., d(a+d)(b+d)(c+d)}` —
     and the latter **is precisely the known untwisted tor2228 family**
     (spot-checks: `[1,55,99,125], [2,4,23,46], [4,11,16,44]` all in
     `data/tor2228.txt`).  A fast production enumerator for it now exists
     (`identity2248_prod.c`, validated against brute force at B=400,
     40/40 exact), ready to extend `tor2228` beyond its old d<=16384.

## Exact commands

```
# enumerate (in scratchpad/t2248prod)
gcc -O2 -march=native -o twisted2248_prod twisted2248_prod.c -lpthread -lm
./twisted2248_prod 30000 2 30000 3 v30000            # 275 s
for spec in "2 63300 c1" "63300 84500 c2" "84500 100000 c3"; do
  set -- $spec; ./twisted2248_prod 100000 $1 $2 3 v100k_$3; done
# sweep (theorem-reduced 4-condition test + full-delta cross-check on sample)
python3 fastsweep.py v30000.txt 100
# task B
magma -b e3_close.m
# jackpot template (if a hit appears)
magma -b a:=A b:=B c:=C d:=D jackpot_verify.m
```

## Results table

| run | tuples | (2,2,4,8) hits | note |
|---|---|---|---|
| twisted family B=3000 (prior) | 1,947 | 0 | prior session |
| twisted family B=30000 | 29,426 | 0 | redundant below 65535 by the theorem, confirms it |
| twisted family B=100000 | see below | see below | fresh territory is d in (65535, 100000] |
| genus-3 split curve | C(Q) = 8 degenerate pts | — | S' split locus CLOSED |

## Resume state

* Scratchpad `t2248prod/`: `twisted2248_prod.c` (+ binary), `identity2248_prod.c`
  (+ binary), `fastsweep.py`, `identity_brute.py`, `e3_close.m/.out`,
  `jackpot_verify.m`, `fresh_check.m`, lists `v200/v1000/v3000/v10000/v30000.txt`,
  `v100k_c*.txt`.
* To go deeper: rerun chunks at larger B (cost ~ B^2.2-2.4; B=200000 est. 5-8 h
  on 3 threads); only tuples with `d > 65535` need the 4-condition test
  (`fastsweep.py`), and any hit goes through `jackpot_verify.m` (prime-scan
  simplicity certificate).
* To extend tor2228: `./identity2248_prod B 2 B nthreads out`.
* The tor2244-side deep search (`d > 65535`) remains the complementary route:
  by the theorem it is THE only route on this model; sign-viability analysis of
  the 48 component x twist tests (which twists can ever be all-positive on
  sorted tuples) would cut that work further — method demonstrated here for
  component I, not yet done for II/III (H_AB delta signs are tuple-dependent).

## Addendum (2026-07-18, post-handoff): B=100000 enumeration COMPLETE
The resumed background run finished: 114,425 twisted-family tuples to d<=100000
(43,158 with d>65535; 13,646 gcd-primitive), banked in
data/claude_prod_07_2248_twisted_B100000.txt. Subset consistency vs B=30000: MATCH.
15-class delta sweep (sign-reduction-theorem form, 229-tuple full-delta cross-check,
0 mismatches): 0 hits; satisfied-conditions histogram {0: 111379, 1: 1416, 2: 1630} —
nothing reaches 3/4. The twisted second component of (2,2,4,8) is now empty to
d<=100000 (the theorem covers d<=65535 for ALL components; this extends the twisted
slice 1.5x beyond it).
