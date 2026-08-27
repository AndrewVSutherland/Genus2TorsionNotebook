# [2,2,2,12] REALIZED (verified hit on the T5 pencil) + T3 route decided by a proven square-class identity

**JACKPOT: a geometrically simple genus-2 Jacobian over Q with torsion
EXACTLY (Z/2)^3 x Z/12 (order 96) was found and fully verified — see
Section 0a and data/claude_prod_02_22212_hit.txt.**

Date: 2026-07-18.  Production agent, target (2,2,2,12) (order 96).
Scripts/outputs: scratchpad `prod_02_22212/`.  Inputs:
notes/claude_top10_04_22212.md, notes/claude_tier1_item3_22212_diagnostic.md,
code/tor22212.c.  Companion data file:
data/claude_prod_02_22212_structures.txt (identity, normal forms, member
tables, search tallies).  No repo file modified.

## 0. Strategy recap (3 lines)

The T3 near-miss locus of M(2,2,2,6) fibers over the pencil
lambda·B3+mu·B5 (rho = B3/B5); Task A was to find where the 4th
condition V5 degenerates.  Answer: it degenerates NOWHERE special —
instead a global identity makes V5 redundant-or-fatal on every member,
splitting the pencil into a dead part and a codim-1 "w-locus" where
near-miss = hit; both parts were then attacked directly (symbolic scan,
exact integer sweeps, Mordell-Weil lifts, rank certificates).

## 0a. THE HIT

`t5sweep.c` (generic (u,rho') box N=150, R=2000, exact integer square
tests of the four T5 normal forms; validated by re-finding all five
known T5 near-misses) found u = -97/48 on member rho' = B4/B5 = -49/240
with ALL FOUR conditions square.  Magma verification (`hitverify.m`):

* (s,m,n) = (336396, -689185, -166464) on M(2,2,2,6);
  curve y^2 = prod(A_i + B_i x), A=[1,1,1,2,2],
  B = [282322361376, -8243383980, -64241207724, -114724491840, 561915878400];
* TorsionSubgroup(Jac) invariants = **[2, 2, 2, 12]** (exact, order 96);
* geometric simplicity: certificates at p = 37, 73, 113 (each chi
  irreducible deg 4 AND chi^12 irreducible deg 4 — multi-prime protocol);
* reduced minimal model: y^2 + (x^2+1)y = 756900x^6 + 737595570x^5
  + 150572203590x^4 - 15854483576121x^3 - 530648977741620x^2
  + 32014154874551031x + 830742747091037849.

A second representation of the SAME curve (identical G2Invariants) was
found by the targeted rn=+-a^2 sweep at u = 133/145, rho' = 289/240
(s,m,n) = (134995,-263568,147175) — the T5 pencil carries an involution
swapping the two members (analog of the proven T3 w<->1/w symmetry).
Empirical laws: all 8 known T5 near/hit points have rho'-numerator
+-odd^2; both hit members have denominator 240.  Post-hit scans:
member -49/240 has no second point to u-height 4000; targeted box
(a <= 61, rd <= 8000, u-height 150) contains exactly the two
representations.  This was the first (2,2,2,12) realization; whether
more curves exist (e.g. rd != 240 members) is now a focused question.

## 1. Headline results

1. **PROVEN IDENTITY** (exact `issquare` in Q[rho,u], `pencil_id.gp`):
   on the T3 pencil, **V1·V2·V4·V5·(1-2rho) is a perfect square**.
   Hence: (a) on any member with 1-2rho not a rational square,
   (2,2,2,12) T3-hits are **impossible** — this retroactively explains
   all 59 near-misses/0 hits (their cores 46,61,109,241,249... are
   exactly sqfree(1-2rho); numeric pin at (25,-26,-15): rho=117/722,
   1-2rho=4·61/361, V5 class 61); (b) on the **w-locus**
   rho=(1-w^2)/2, near-miss <=> hit.
2. **rho-scan (Task A as specified) completed** (`rhoscan.gp`, log in
   scratchpad): all degeneration loci of the V5 pullback S5 (deg-7,
   contains q(u) and one linear factor of each of V1,V2,V4 — the source
   of the identity) have rational roots only at rho in
   {0, 1/2} (lc/disc: the two degenerate members B3=0, V5=0) and
   {1, -1/2, 1/4} (S5(u=0), i.e. the excluded base point s=0).  **No
   special member exists**; confirmed independently by the elliptic
   surface discriminant Delta ~ rho^2(4rho^2-2rho+1)^2 (`jinv.gp`).
3. **The w-locus is uniform, unkilled, and empty as far as searched:**
   * Exact normal forms (constants validated on 117 random points):
     V1 == q·(1-(4u-2)rho), V2 == (2u-1)q·((2u-1)-(4u-4)rho),
     V4 == (u-1)q·(rho+u-1) mod squares, q = 4u^2-6u+3.
   * On EVERY member tested (29 members w=a/b, a,b<=7 in `members.m`;
     +22 more in `liftdrv`), the three condition covers are the SAME
     elliptic curve (E1=E2=E4 as minimal models), rank 1-2, never 0,
     torsion Z/2; w<->1/w gives the same curve;
     j(w) = 64(w^2-2)^3(2w^2-1)^3/((w-1)^2(w+1)^2(w^4-w^2+1)^2).
     **No rank-0 kill anywhere.**
   * Universal degenerate section u* = (1+2rho)/(4rho): the unique
     rational 2-torsion point of E1, lies on B2=0 (curve degenerates,
     `secchk.gp`) — it is the ONLY point the MW boxes ever produced.
   * **Searches (all 0 genuine hits):** blind exact-integer sweep
     `t3sweep.c` (validated vs PARI `xval.gp`): N=1000/W=100 in 19 s
     and N=4000/W=150 in 727 s (u=p/q up to 4000, w=a/b up to 150;
     51078 single-condition passes, 0 triple);  MW lifts `lift.m`/
     `liftdrv.sh`+`lift_one.m`: |n_i|<=12 (21 members), |n_i|<=16
     (all coprime a,b<=8), NMAX=200 deep runs on the two rank-1
     members (w=2,1/2) — zero points passing a second condition
     beyond the degenerate section.
3b. **Genus-3 quotient probe (`h7.m`, `h7b.m`)**: for ALL 22 members
   w = a/b with coprime a,b <= 6, the quotient H7_w: y^2 =
   q(u)·(five linear factors) (deg 7, genus 3) — into which every hit
   injects with y != 0 — has NO rational point with y != 0 up to
   height 1e5 (only the 6 Weierstrass points).  These members are
   hit-free to u-height 1e5; each is a clean Chabauty target (finish
   proof next session).
3c. **Global class obstructions (`identity2.gp`)**: unrestricted on
   M(2,2,2,6), prod of the 4 class-T_i conditions == c_i·C_i(s,m,n)
   mod squares with C_3 cubic-in-s (c=-1), C_5 cubic (c=8) — the
   "thin" classes — and deg-7/8 forms for T1/T2/T4 (stiff; explains
   why only T3/T5 near-misses occur).  Master obstruction surfaces:
   t^2 = -C3, t^2 = 2·C5 (double covers of P^2); the pencil identity
   is the restriction of the first to the conic pencil.
4. **T5 pencil**: V'1V'2V'3V'4 == (q·rho'-2u+1)(q·rho'-1) mod squares
   (proven; NOT a constant class — T5 has no analogous constant
   obstruction; the failing-class observed at (80,125,8) matches).
   Task B rank checks (`taskB.m`, `taskB5.m`): rho'=-1/143: EW3 rank
   [1,1] tors Z/4, EW4 rank [2,2]; rho'=-25/551: EW3 [2,2], EW4 [2,2];
   the two conic conditions W1,W2 are rational (hundreds of small
   points).  No kill; route open but 2-dimensional (see data file).
5. **Task B (T3 third 6-point member rho=-63/242)**: E1=E2=E4 =
   [0,1,0,-4775,563754] (cond 4407144), **rank exactly 3**, tors Z/2;
   pairwise products all [0,1,0,-199934824,992407557956]
   (cond 2230014864), rank 3, tors (2,2); IsIsomorphic confirmed.
   With the identity, the "correlated triple" mystery of the dossier is
   now explained structurally (one elliptic curve per member + constant
   twist class 1-2rho).
6. **Task C enumerator (chart M(2,2,2,4)+3, y^2=x(x+a^2)...(x+d^2))**:
   H=2000 chunked over d, cached -P127 tables (`t127.bin`), plus NEW
   `postfilter.c` (direct 96 | #J(F_p) for all good p in [131,397]).
   Results (ALL FOUR CHUNKS COMPLETE): d in [5,1000]: 1877 survivors
   -> 0 pass; [1001,1400]: 6396 -> 0; [1401,1700]: 11146 -> 0 (1560 s);
   [1701,2000]: 20369 -> 0 (2261 s).  Total 5.6e11 tuples checked,
   39788 first-stage survivors, ZERO pass the p<=397 postfilter.  Since
   96 | #J(F_p) at good p is NECESSARY for torsion >= (2,2,2,12), this
   is an **unconditional kill for ALL gcd-1 tuples a<b<c<d <= 2000**
   in the chart y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2) (vs the
   ~1-at-H=2500 heuristic: true count 0; the -P127-only survivor
   density is ~50x the p<=199 prediction, everything killed by
   p in [131,167] in practice).

## 2. Exact commands

```
cd /tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/prod_02_22212
gp -q pencil_id.gp        # identity proof + numeric pins (~10 s)
gp -q constants.gp        # exact normal forms + T5 analog (~10 s)
gp -q xval.gp             # 117-point random validation of C formulas
gp -q rhoscan.gp          # Task A degeneration scan (log: rhoscan.log)
gp -q jinv.gp             # j(rho), j(w), Delta of the elliptic surface
gp -q secchk.gp           # section u*: B2 == 0 proof
gcc -O3 -march=native -fopenmp -o t3sweep t3sweep.c -lm
./t3sweep 4000 150 2      # decisive blind sweep (727 s, 2 threads)
magma -b members.m        # per-member E1=E2=E4 + ranks (29 members)
magma -b wa:=A wb:=B nmax:=N lift_one.m   # MW lift of one member
bash liftdrv.sh           # all coprime a,b<=8, timeout 150 s each
magma -b taskB.m ; magma -b taskB5.m      # Task B rank certificates
bash runchunks.sh         # sieve chunks (skips done_*; 3 threads)
gcc -O3 -march=native -o postfilter postfilter.c -lm
./postfilter < surv_dLO_HI.txt > pf_dLO_HI.txt 2> pf_dLO_HI.log
```

## 3. Resume state

* Sieve: **ALL chunks d <= 2000 done and postfiltered (0 pass)** —
  markers done_d5_1000 ... done_d1701_2000 + ALLDONE set.  Next
  extension: H=3000+ (edit runchunks.sh; keep -P127 + postfilter
  pipeline, strictly stronger and cheaper than -P199 tables).
* T5 sweeps DONE: main box N=150/R=2000 (2341 s, 1.03e7 X1-passes,
  2 full hits = 1 curve); `--sq` box a<=61/rd<=8000/N=150; `--member`
  scan of -49/240 to u-height 4000.  Next boxes: bigger N (u-height)
  in main mode; rd sweep beyond 8000 in --sq mode.
* `lift5.m` (MW enumeration on the hit member via quartic V'3): both a
  400 s and a 1500 s attempt died in Rank/Generators (the E3 quartic
  Jacobian has large coefficients) — use TwoDescent/points-search with
  bounded effort, or work on the reduced minimal model, to hunt sibling
  hits at exponential heights next session.
* Deep lifts `liftdeep_2_1.log`, `liftdeep_1_2.log` (NMAX=200, w=2,1/2):
  802 pts each, 0 passes.
* Candidate validator: `candfull.gp` (T3), `hitverify.m` (both classes,
  arg-driven: `magma -b cls:=5 pp:=.. qq:=.. rn:=.. rd:=.. hitverify.m`).

## 4. Next steps (ranked)

1. **Close the w-locus**: the residual question is rational points on
   the (Z/2)^2-fiber-products X_w of three isomorphic-but-distinct
   double covers of E(w).  Since E(w) is now explicit with j(w) known,
   compute the genus-3 quotient H7_w: y^2 = q(2u-1)(u-1)·L1L2L4 and try
   TwoCoverDescent/Chabauty (Jac(H7) decomposes; rank data suggests
   rank-vs-genus favorable for some w).  A single w with a point kills
   the target; a uniform obstruction proof retires the T3 route.
2. T5 route: hits live over the rational surface S: z^2=(qr-1)(qr-2u+1)
   (Pythagorean: z^2+(u-1)^2 = (qr-1-(u-1))^2, parameter g), and ON S
   the identity gives V'4 == V'1·V'2·V'3 mod squares — i.e. only THREE
   independent conditions on a 2-parameter (u,g) family (same shape as
   the pre-identity T3 near-miss situation, one more degree of
   freedom).  Build the (u,g) sweep with r = r(u,g); the generic-(u,r)
   box sweep `t5sweep.c` (validated: re-finds all three known T5
   near-miss classes) is running N=150/R=2000 — see results below.
3. Enumerator: finish chunks to d=2000, then H=3000+ with -P199 tables
   on a quiet box (postfilter makes -P199 unnecessary: keep -P127 +
   postfilter pipeline, it is strictly stronger and cheaper).
4. Other classes T1/T2/T4 on M(2,2,2,6): run the same pencil/identity
   machinery (B-pencils lambda·B1+mu·B5 etc.) — 30 min of gp each now
   that the method is scripted.
