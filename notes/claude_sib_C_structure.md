# Lane C (structure): T5 pencil fully factored, involutions proven, LAW proven —
# and the law-targeted sweep produced the SECOND (2,2,2,12) curve

Date: 2026-07-18.  Sibling lane C of the (2,2,2,12) campaign (after the first
hit, notes/claude_prod_02_22212.md).  Scripts: scratchpad `sib_C_structure/`
(verify.gp, step2.gp, step3.gp, chartmatch.gp, sigmob.gp, t3fix.gp, lawsweep.c,
fixsearch.py, batch5.m, minmodel.m).  Data out:
**data/claude_sib_C_curve2.txt** (the new curve, all reps, certificates).

## 0. HEADLINES

1. **SECOND geometrically simple (2,2,2,12) curve /Q found and fully verified**
   (torsion exactly [2,2,2,12], 3 simplicity certificates at p=71,103,127):
   representative (u,rho') = (-23/75, -9025/3519), (s,m,n) = (2208,-8303,-7200),
   minimal model y^2 + (x^2+x)y = 36x^6 + 36750x^5 - 462983772x^4
   - 301623595823x^3 + 1518598238654317x^2 + 397058962729817115x
   - 1282993930035013443975.  Six pencil representations known (data file).
   Found in 71 s by the LAW sweep (below) — vs ~40 min for the blind box that
   found curve 1.
2. **Complete multiplicative structure of the T5 conditions** (all proven by
   exact square-class computation, verify.gp/step2.gp): with v = 2u-1,
   Q = q(u) = v^2-v+1, X = Q*rho':
   V'1 == Q·X·(X-v),  V'2 == Q·(X-Q)·(X-v),
   V'3 == v·Q·X·(X-1)·(X-Q)·(vX-Q),  V'4 == v·Q·(X-v)·(vX-Q)   (mod squares).
   Product == (X-1)(X-v) — recovers the prod_02 identity.
3. **THE LAW (proven): V'1·V'2 == rho'(rho'-1) mod squares.**  Hence any point
   passing conditions 1,2 (in particular every near-miss and every hit) has
   rn(rn-rd) = a nonzero square: **rn = ±a^2 and rd = ∓(c^2-a^2)** (signs
   matched, gcd(a,c)=1).  Corollaries: no hit has rho' in (0,1); the empirical
   "denominator 240" is just (a,c)=(7,17); the empirical "odd squares" law is
   FALSE — curve 2 has a member with (a,c)=(182,209), a even.
4. **Involutions (proven, explicit):**
   * tau = relabeling (1 2) of the A=1 slots (identity Moebius, exact same
     polynomial): **tau: (u,rho') -> ((4rho'u-3rho'-1)/(4rho'(u-1)), rho')** —
     member-preserving.  Explains all same-member companion pairs, and
     predicted the 3rd representation of curve 1 at u = 3637/7105 (verified
     all 4 conditions; height 7105 explains why the u-height-4000 member scan
     missed it).
   * sigma = relabeling (1 2)(5 inf) (slot 5 <-> Weierstrass point at
     infinity; found numerically by Moebius-matching, sigmob.gp):
     **sigma: (u,rho') -> ((4u-3)/(4u-4), 1-rho')**, i.e. v -> v/(v-1),
     rho' -> 1-rho'.  Maps the hit pair members -49/240 <-> 289/240.
   * <sigma,tau> = (Z/2)^2; sigma·tau = (5 inf): (u,rho') ->
     ((rho'u-1)/(rho'-1), 1-rho').
   * Search-halving: every hit is sigma-equivalent to one with rho' < 0, i.e.
     member rho' = -a^2/(c^2-a^2).
5. **A THIRD moduli identification exists (empirical, not yet explicit):**
   G2-dedupe shows hit members come in TRIPLES per curve:
   curve 1: (a,c) = (7,17), (13,77), (13,187); curve 2: (95,112), (65,88),
   (182,209).  Each member also carries its tau-pair of u-values.  The 6
   affine S3-relabelings and the 6 (45)-composed ones were computed
   (chartmatch.gp, step3.gp): only (12) is an affine chart symmetry, and
   NO (45)-composed relabeling is — so the triple-identification is NOT a
   slot relabeling with affine Moebius; deriving it explicitly (probably an
   order-6-point translation P -> P+T, or an inf-moving relabeling) is the
   top open item — it would cut future sweeps by another factor 3.
6. **On the surface S** (z^2 = (X-1)(X-v), parameter t: X = (t^2-v)/(t^2-1);
   hits all lie on S): conditions reduce to THREE conics over the v-line:
   C1 = Q(1-v)(t^2-v),  C2 = Q(v t^2-v+1),  D4 = vQ((v-1)t^2+1)
   (V'3 == C1·C2·D4 on S is then automatic).  sigma acts as t -> 1/t and
   swaps C1<->C2, fixes D4.  The known hit: v=-121/24, t=437/1013.
7. **Fixed-locus routes are provably DEAD (both classes):**
   * fix(tau_T5): rho' = -1/(v(v-2)).  Conditions collapse to v+1 = A^2,
     v-1 = B^2, plus one more; parameterizing v = (w^4+4)/(4w^2) leaves the
     genus-5 curve Y^2 = 2(w^4-2w^3+4w+4)(w^4-2w^2+4)(w^4+2w^3-4w+4), which
     has **no Q_3-points (unsolvable mod 27; fixsearch.py: no rational points
     to height 2500)**.
   * T3 pencil analog (t3fix.gp): tau_T3: u -> ((4r-1)u-3r)/((4r-2)u+1-4r)
     fixing rho; fix locus rho = (v^2-1)/(2v(v-2)).  Conditions collapse to
     2-v = A^2, v(v-2) = B^2, 1-2v = E^2  =>  v = -D^2, A^2-D^2 = 2,
     E^2-2D^2 = 1  =>  Y^2 = 2(h^4-2h^2+4), which is **unsolvable in Q_3
     (exact valuation-1 mod 9 / 8 mod 9 argument)**.  So the T3 tau-fixed
     locus contains no hits — consistent with (and structurally explaining
     part of) the total emptiness of the T3 route.

## 1. The law-targeted sweep (lawsweep.c)

Enumerate ONLY members rho' = -a^2/(c^2-a^2), coprime a<c (both parities!),
u = p/q with |p|,q <= N; X1-prefilter then all four exact i128 square tests
(code adapted from t5sweep.c; validated: --test refinds the hit; c<=61/N=500
box refinds all 5 known near-misses + hit in 0.8 s).

Results (3 threads):
* c<=301, N=1000: 71 s, 9 hit-representations, 10 near.
* c<=601, N=1000: ~310 s, same 9 hits (no new curve with larger members
  at u-height <= 1000).
* All 9 = exactly TWO curves after G2-dedupe (batch5.m):
  curve 1 (the known hit) as members (7,17),(13,187),(13,77);
  curve 2 (NEW) as members (95,112),(65,88),(182,209).
* Every representation independently re-verified: torsion exactly [2,2,2,12]
  + 3 simplicity certificates each (chi irred deg 4 AND chi^12 irred deg 4).

New near-misses found (3/4 conditions, X4 fails), all obeying the law:
(-383/4, -361/69335), (-9/8, -8836/14573), (-787/13, -121/14279),
(-527/58, -3136/34889), (125/242, -121/14279) + the 5 previously known.

## 2. Exact commands

```
cd .../scratchpad/sib_C_structure
gp -q verify.gp     # factored forms, LAW, product identity (~1 min)
gp -q step2.gp      # sigma transport, S-surface conditions, hit coords
gp -q chartmatch.gp # affine relabeling symmetries: only (12) => tau formula
gp -q step3.gp      # <sigma,tau> group, fix(tau) collapse, (45)-tests
gp -q sigmob.gp     # sigma = (12)(5 inf) via numeric Moebius matching
gp -q t3fix.gp      # T3 tau + fix locus (then hand/py 3-adic kill)
gcc -O3 -march=native -fopenmp -o lawsweep lawsweep.c -lm
./lawsweep 301 1000 3            # 71 s: 9 hit reps = 2 curves
magma -b batch5.m                # torsion + certs + G2 dedupe of all reps
magma -b minmodel.m              # curve 2 minimal model/Igusa/disc
python3 fixsearch.py 2500        # fix(tau_T5) point search: 0 points
```

## 3. Resume state

* lawsweep boxes DONE: (c<=61,N=500), (c<=301,N=1000), (c<=601,N=1000; 306 s),
  (c<=101,N=3000 deep-u; 78 s, 4 hit reps all known) — outputs law_*.txt in
  the scratchpad, consolidated in data/claude_sib_C_lawsweep.txt.
  No processes left running at session end.
* Next boxes: c<=1501/N=1000 (bigger members), c<=61/N=10000 (deep u; i128
  bounds checked to N~2000 for c~500 — recheck overflow before N>4000 with
  large c), and the (v,t)-surface direct sweep (t-heights ~1000+, since the
  known hit has t = 437/1013).
* OPEN (top): derive the third identification (member triples) explicitly;
  candidates: inf-moving relabelings ((5 inf) analogs on other slots) or
  order-6-generator changes.  With it, sweep only 1 member per 6.
* OPEN: is the number of (2,2,2,12) curves infinite?  The law reduces the
  question to squares of the three explicit conic classes on the (v,t)-plane;
  no rational curve on the triple-cover surface found yet (fix loci are dead;
  try curves where C1 or C2 degenerates, or conic-bundle sections over
  sub-families v = phi(s)).
* OPEN (Task 3 spec): T1/T2/T4 classes — (45)-affine relabelings are NOT
  chart symmetries, so T4 is not trivially T5-equivalent; run the same
  factored-normal-form derivation on the T4 pencil rho'' = B5/B4 and the
  deg-7/8 obstruction forms of identity2.gp; expect an analogous X-line
  factorization since the (v,X)-frame is class-agnostic.
```

## Addendum (Codex review, PR #4): checked-in scripts
The exact-commands sections above reference the session scratchpad. The key scripts are
now checked in under code/claude_sib_lanes/C/ (same filenames); sweep binaries rebuild
with the gcc lines given in the commands. Scratchpad paths remain valid only on the
discovery machine.
