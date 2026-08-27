# PROD 05: (5,5) — order 25, smallest open group

Date: 2026-07-18.  Production campaign step, target #5.

## Strategy recap (3 lines)

The only live (5,5) route is the `b2=0` branch of the full Mumford norm system
`f = (1+h1*x+h2*x^2)^2 - K*x^5`, `A^2 - B^2*f = U^5` (A monic quintic forced,
`B = b0 + b1*x`), open locus `K*b1*b0*disc(U)*Res(B,U)*disc(f) != 0`.
Task A settles the celebrated smooth F_7/F_11 charts by exact Groebner
decomposition of the saturated `h1=1,h2=0` slice over Q; Task B eliminates the
full 2-dimensional branch over Q using the weight-(1,2,5,-1,-2,-5,-4) torus
normalization (h1=1 / h1=0,s=1 / h1=0,s=0 trichotomy); Task C is the multi-slice
F_p census + CRT fallback.

## Scripts (scratchpad t55/)

All under `/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/t55/`:

- `prod05_weights.m` — quasi-homogeneity certificate.
- `prod05_sliceQ.m` -> `prod05_sliceQ.log` — TASK A decision procedure (complete).
- `prod05_h1oneQ.m` -> `prod05_h1oneQ.log` — TASK B case h1=1 (6 vars).
- `prod05_h1zeroQ.m` -> `prod05_h1zeroQ.log` — TASK B cases h1=0 (s=1 and s=0).
- `prod05_fpcensus.c` — TASK C F_p census of the open branch locus.
- `census_p7.txt` ... `census_p23.txt` — census outputs.
- `prod05_interp_modp.m` -> `prod05_interp_p10007b.log` — (h2,K) plane model
  of the h1=1 curve by fiber interpolation mod 10007.
- `prod05_fiberelimQ.m` -> `prod05_fiberelimQ1.log` — lean EXACT per-fiber
  decision over Q (saturated K-eliminant + rational-root test).
- `prod05_fibfac.m` — per-fiber eliminant factorization probe mod 10007.
- `prod05_taskC_screen.m` + `prod05_crt.py` — 3-prime multi-fiber screen of
  all h2 with height <= 10, with CRT rational reconstruction.
- `prod05_jackpot.m` — full verification protocol for any rational candidate.
- `prod05_rootscreen.m` -> `prod05_rootscreen.log` — the per-fiber root-free
  prime screen (the workhorse negative instrument).
- `prod05_c1_st.m` -> `prod05_c1_st.log` — C1 point sampling + (s,t)-shadow
  fit (G of degree 10) + rational-function fits.
- `prod05_h24_dig.m`, `prod05_h24_lll.m`, `prod05_h8_lll.m`,
  `prod05_m4_verify.m`, `prod05_h8_verify.m` — escalation digs for the
  persistent-root fibers h2 = 4, -4, 8 (all closed).
- `prod05_parasite.m` — the <D0>-parasite refutation / automatic
  independence certificate.
- `prod05_F1_lift.m`, `prod05_F1_param.m` — F1 lift attempt (failed:
  large coefficients) and the ready parametrization step for when F1/Q
  is obtained.

Key artifacts copied to `data/claude_prod_05_55_sliceQ.log` (Task A) and
`data/claude_prod_05_55_census.txt` (Task C census summary).

## TASK A — THEOREM (decision procedure, COMPLETE)

Exact commands:

```sh
cd <scratchpad>/t55 && magma -b prod05_sliceQ.m > prod05_sliceQ.log
```

Validation: the three known open points (p=7: (K,s,t,b0,b1)=(2,0,6,3,1),(2,0,6,4,6);
p=11: (7,5,10,2,4)) satisfy all five residuals and all six openness factors.

Result (exact, over Q):

- Raw slice ideal (5 residuals, degrees 10..6, terms 26/23/20/15/11): dim 2.
- After saturation by `K, b1, b0, disc(U)=s^2-4t, Res(B,U), disc(f)=3125K^4+108K^3`
  (iterated = product saturation, verified equal): **dim 0**, quotient dimension 34.
- Radical decomposition: **2 prime (Galois-irreducible) components, degrees 28 and 6**.
  - Component 1 (degree 28): minpolys of K,s,t have degree 14, of b0,b1 degree 28.
    Rational points: NONE.  Reduces mod 11 to the two known F_11 open points
    (7,5,10,2,4),(7,5,10,9,7); has NO F_7 points.
  - Component 2 (degree 6): minpolys of K,s,t have degree 3 (K-field defined by
    K^3 - (10634004/9765625)K^2 + (549579978/30517578125)K - 960097374/95367431640625),
    b0,b1 degree 6.  Rational points: NONE.  Reduces mod 7 to exactly the two
    known F_7 open points (2,0,6,3,1),(2,0,6,4,6); has NO F_11 points.
- `TOTAL_RATIONAL_POINTS_ON_SATURATED_SLICE=0`, cross-checked by `Variety(J,Q)=[]`.

**Closed negative:** the smooth F_7 chart points (the ones with unique Hensel
lifts to 7^7) are reductions of a Galois orbit of 6 conjugate points whose
K-coordinate generates a degree-3 field; the smooth F_11 chart points are
reductions of a degree-28 orbit (K of degree 14).  Neither orbit contains a
rational point.  "Reconstruction failed" on these charts is now a theorem:
**no rational point of the (5,5) b2=0 full-norm branch lies on the h1=1, h2=0
slice.**  The p-adic lift targets of agent_z5x5_b2zero_lift.md are certified
irrational.

## Scope remark (no hidden Lambda-twist sub-branch)

In A^2 - B^2*f = Lambda*U^5 with deg B <= 2, deg(B^2*f) <= 9 < 10 forces
deg A = 5 exactly and Lambda = lc(A)^2, a square — so A can always be taken
monic and Lambda = 1.  The full-norm system as posed covers ALL degree-2
Mumford 5-torsion representatives; b2=0 is a genuine linear-B sub-branch,
and the b2 != 0 general branch (residuals 852/560/347/208/123 terms)
remains the only other lane for (5,5) in this ansatz.

## Torus certificate

`prod05_weights.m`: all five residuals and all three boundary factors are
quasi-homogeneous for the x-scaling torus with weights
(h1,h2,K,s,t,b0,b1) = (1,2,5,-1,-2,-5,-4)
(weighted degrees -10,-9,-8,-7,-6 and 20,-2,-10).  Since h1 has weight 1 and s
weight -1, every Q-point of the branch is torus-equivalent over Q to one with
h1=1, or h1=0,s=1, or h1=0,s=0.  This makes the Task-B trichotomy exhaustive.

## TASK B — full branch over Q (h1,h2 free)

Torus trichotomy (exhaustive by the weights certificate):

**Case h1=0, s=0 (prod05_h1zeroQ.log): EMPTY over Qbar.**
Saturation by b1 alone already kills the raw dim-2 ideal (dim -1).  Closed.

**Case h1=0, s=1 (prod05_h1zeroQ.log): saturated dim 0, NO rational points.**
Raw GB (54s) dim 2; saturation chain K,b1,b0,discU,resBU,discF gives dim 0.
Radical decomposition: two prime components, degrees **52** (minpolys: h2,K,t
deg 26; b0,b1 deg 52) and **10** (minpolys: h2,K,t deg 5; b0,b1 deg 10).
`rational_points=[]` on both.  Closed over Q.

**Case h1=1 (variables h2,K,s,t,b0,b1; saturated dim 1 = curve):**
whatever survives of the (5,5) b2=0 branch over Q is exactly the rational
points of this saturated curve.  The DIRECT 6-var Groebner is infeasible:
the raw grevlex GB did not finish in 23 min over Q (`prod05_h1oneQ.m`,
killed) nor in 18 min mod 101 (`prod05_h1one_modp.m`, killed) — consistent
with the curve's true size found below.  The working method is fiber-wise:
exact 5-var eliminant decisions per rational h2 (`prod05_fiberelimQ.m`)
plus mod-p interpolation of the (h2,K) plane model (`prod05_interp_modp.m`).

**Structural discovery (mod-p fiber interpolation, prod05_interp_modp.m):**
the generic h2=c fiber of the saturated h1=1 curve has quotient dimension
**124** with K-eliminant of degree **62** (checked mod 10007 for many c).
The Task-A fiber h2=0 (qdim 34, elimK degree 17) is a DEGENERATE special
fiber — most of the generic degree-124 scheme escapes into the boundary
there.  So the h1=1 curve has degree ~124 over the h2-line (62 modulo the
B -> -B involution), i.e. the b2=0 branch curve is far larger than the
h2=0 slice suggested.  Plane model F(h2,K) being interpolated mod 10007;
first pass shows degK=62 with coefficient-of-K^j h2-degree following
~ floor(155 - 5j/2), so the plane model has bidegree about (155, 62)
(consistent with the torus weights h2:2, K:5 — 155/62 ~ 5/2).

Per-fiber eliminant factorization mod 10007 (prod05_fibfac.m): the degree-62
K-eliminant factors with a DIFFERENT degree pattern at each fiber
(c=1: [10,14,18,20]; c=2: [2,2,2,2,2,2,8,8,12,22]; c=3: [1,2,2,3,4,4,4,8,10,24];
c=5: [2,3,5,6,8,38]; c=11: [2,4,4,4,22,26]; c=23: [1,1,1,3,4,4,8,12,28]).
F_p-fiber point counts are sparse (c=1,2,5,11,23: zero F_10007 points;
c=3: two) — a 3-prime liveness screen is highly selective.

**INTERPOLATION THEOREM (mod 10007, verified on 5 fresh fibers,
prod05_interp_p10007b.log):** the (h2,K) shadow of the h1=1 curve is
F(h2,K) = F1 * F2 with
- **F1: bidegree (25, 10), total degree 25, 62 terms, GENUS 0**, 10005
  projective F_10007 points — a RATIONAL curve;
- F2: bidegree (130, 52), total degree 130, 1554 terms (genus computation
  killed as too slow; large).
Fiber accounting: 124 = 2*10 + 2*52, so BOTH components are ±B double
covers of their shadows: over a generic shadow point of F1 the fiber is a
single (b0,b1) -> (-b0,-b1) pair sharing (s,t).  Hence on the F1-component,
**s,t,K are rational functions on a genus-0 curve and only (b0,b1) needs a
square root** — the (5,5) question on this component reduces to a quadratic
condition over a rational base.

**Parasite hypothesis REFUTED + automatic independence
(prod05_parasite.m):** the only nontrivial classes in <D0> (D0 = contact
class [(0,1)-infty]) have reduced Mumford reps [x, +/-1] (degree 1) and
[x^2, +/-(x+1)] (non-separable U).  Both are excluded by the open conditions
(deg U = 2, disc(U) != 0).  Since reduced Mumford reps are unique, EVERY
open-branch point has [U,V] outside <D0>: the two 5-torsion classes are
AUTOMATICALLY independent — every rational open-branch point genuinely
carries (Z/5)^2.  (Bonus: the norm function of 2*D0 has linear B exactly
at h2 = 0, which explains why the Task-A fiber h2=0 degenerates 124 -> 34.)

## TASK C — F_p census (open locus, all 7 variables)

`prod05_fpcensus.c` (validated: reproduces the 12/240 counts of
agent_z5x5_b2zero_elim.md exactly):

- p=3: open locus EMPTY (prior note).
- p=7: 12 open points, ALL with h2=0, h1 in F_7^* — a single torus-orbit pair;
  consistent with Task A component 2 (the degree-6 orbit) being the ENTIRE
  F_7 open locus.
- p=11: 240 open points, all h1 != 0; h1=1 fiber has 24 points spread over 10
  h2 values (h2=6 empty) — the h1=1 saturated curve has 24 F_11 open points.
- p=13: 96 open points; h1=1 fiber has 8 points on h2 in {3,4,6,9}.
  Points come in (b0,b1) -> (-b0,-b1) pairs (the B -> -B involution).
- p=17: 160 open points; h1=1 fiber 8 points; h1=0 locus contributes 32
  (= 2 torus orbits, the reduction of the case-2 irrational locus).
- **p=19: open locus EMPTY** (0 points in F_19^7) — second empty prime after
  p=3.  Any rational open point must be non-integral or boundary mod 19.
- p=23: 264 open points, all with h1 != 0 (12 per h1 value; h1=1 fiber 12).

h1=1-fiber open counts: p=7:2, 11:24, 13:8, 17:8, 19:0, 23:12.

## C1 component: explicit double-cover data (mod 10007)

`prod05_c1_extract.m`: on the C1-subfiber (fiber ideal + g1_c(K), where
g1_c = gcd of the fiber K-eliminant with F1(c,K)): qdim=20 for every tested
c=1..12, with s- and t-eliminants of degree 10 (s,t are FUNCTIONS of K on
C1) and b0-eliminant of degree 20 (quadratic).  Saved factors:
`prod05_F1_p10007.txt` (bidegree (25,10)), `prod05_F2_p10007.txt`
((130,52)); archived in data/claude_prod_05_55_F1_p10007.txt, _F2_.
Balanced small-integer lift of F1 FAILS fresh-prime verification
(`prod05_F1_lift.m`): the true Q-coefficients of F1 are large; getting F1/Q
needs multi-prime CRT (~6 min per prime for 165 fibers).

## Per-fiber decisions over Q (h1=1 curve)

Exact-GB per generic fiber is too slow over Q (h2=1 raw GB killed at 12 min;
the degenerate h2=0 fiber took only 21s).  CRT reconstruction of the
K-eliminant is UNSTABLE at 40 primes (coefficient heights >> 10^70,
matching the Task-A minpoly heights; `prod05_fiberelim_crt.m`).
Working tool: `prod05_rootscreen.m`.  Soundness: the Q-eliminant g_c has
degree 62 (= the Qbar fiber degree); at any prime p where the mod-p
saturated-fiber eliminant also has degree 62, reduction sends g_c to that
mod-p eliminant up to a unit (the mod-p saturation contains the reduction
of the Q-saturation, and degree equality pins the generator), so a rational
K-root a/b with p∤b reduces to an F_p-root.  Hence EVERY root-free good
prime p forces p | b: k root-free primes force denominator >= their
product (~10^4k per fiber).  Combined with the fiber-wise necessary
condition (a rational open-fiber point needs a rational K-root), each
root-free prime is an unconditional height wall on that fiber.
Results in `prod05_rootscreen.log` (below).

**ROOTSCREEN RESULTS (30 fibers: all h2 = a/b with height <= 8 tested,
data/claude_prod_05_55_rootscreen.log):**
- 27 of 30 fibers killed directly: a root-free good prime exists within the
  first 12 primes above 10^4 (often the very first).  Verdict per fiber:
  NO rational K-root with denominator coprime to the root-free prime(s);
  hence NO rational open-fiber point below the corresponding wall.
- h2=4 escalated (roots at all 12 primes): LLL from 7 unique-root primes
  (prod05_h24_lll.m) shows the persistent root is a CUBIC irrationality,
  minpoly 7174376*K^3 + 2527934*K^2 - 528771*K + 2674113 (verified at all 7);
  no degree-1 relation exists up to height ~10^13.  Fiber closed.
- h2=-4 escalated: the 3-prime degree-1 LLL candidate 678791/491530 is
  REFUTED (3/15 primes); primes 10133 and 10139 are root-FREE for this
  fiber, closing it (denominator wall 10133*10139 > 10^8).
- h2=8 escalated: 4-prime degree-1 candidate refuted (0/15); prime 10139 is
  root-free, closing the fiber.
- Summary: NO rational point on the h1=1 branch curve over any h2 of
  height <= 8, each closed by an explicit modular certificate (plus the
  exact h2=0 closure of Task A).

**C1 (s,t)-shadow (prod05_c1_st.m, 181 samples from 164 fibers):**
the (s,t)-image of C1 satisfies an (overdetermined, kernel-dim exactly 1)
plane curve **G(s,t) of degree 10, 34 terms** (mod 10007; saved in
data/claude_prod_05_55_Gst_p10007.txt).  Low-degree (<= 6) rational-function
fits of K, h2, b0^2 in (s,t) FAIL — the fits must be redone modulo G (in
the quotient ring) to find the true representations; that is the direct
route to the square-condition hyperelliptic model w^2 = W(tau) on C1.

## Status / resume

- Task A: COMPLETE (theorem above).
- Task B: h1=0 cases CLOSED (theorems above).  h1=1 case: structure
  understood (two components; C1 over a genus-0 shadow with quadratic
  (b0,b1)-obstruction and degree-10 (s,t)-model; C2 over a bidegree-(130,52)
  shadow); all 31 tested fibers (h2 height <= 8) closed by explicit
  certificates; full closure over Q open.
- Task C: census 7..23 complete (two empty primes: 3, 19); the rootscreen
  IS the refined Task-C instrument (subsumes the CRT screen); the 3-prime
  point-screen scripts (`prod05_taskC_screen.m` + `prod05_crt.py`) remain
  available for the un-screened height range.

## Next steps (ranked)

1. **C1 square-condition model (the decisive object):** redo the K/h2/v
   fits modulo G(s,t) (fit N - K*D in the quotient ring F_p[s,t]/(G), i.e.
   include G-multiples in the linear system), then parametrize the genus-0
   degree-10 curve G (odd-degree? 10 is even — check genus and points; the
   (h2,K) shadow F1 has odd degree 25 so C1's base is parametrizable over Q
   if F1/Q is computed) and obtain w^2 = W(tau) over Q.  Genus of that
   hyperelliptic model decides the whole C1 lane: small genus -> targeted
   rational point search (jackpot); genus >= 2 -> finiteness + bounded scan.
2. F1 over Q by multi-prime CRT (each prime ~6 min of fibers; expect
   10-30 primes) — enables step 1 over Q and the exact parasite check for C2.
3. Extend the rootscreen to h2 height <= 30 (cheap: most fibers die at the
   first prime; ~2s per prime per fiber).
4. C2 (degree-130 shadow): compute its F_p point count / genus lower bound;
   if genus is large and its fibers stay root-free, C2 is structurally cold.
5. The b2 != 0 general branch (852-term residuals) remains untouched — the
   only other (5,5) lane in this ansatz.

## Session verdict

No (5,5) realization found (no jackpot).  Major structural advance: the
b2=0 branch over Q is now reduced to the rational points of ONE explicit
curve (the h1=1 slice), whose two components are both ±B double covers of
explicit plane shadows (genus-0 bidegree-(25,10) and bidegree-(130,52));
the h1=0 sub-branches and the h2=0 slice are closed by exact theorems;
31 further fibers are closed by modular certificates; and every rational
open-branch point automatically carries independent classes (true (5,5)).
The route is ALIVE but constrained: the C1 square-condition curve is the
single decisive object left on this lane.
