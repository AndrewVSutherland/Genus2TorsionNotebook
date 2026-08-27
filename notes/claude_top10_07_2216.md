# [2,2,16] (top-10 #7): both charts tested — twisted-family second halving is clean-negative to d<=1000; d=-2 fibers have rank 2 (Y(Q) infinite) but still 0 field passes

*(Claude agent, 2026-07-17. Scripts and raw outputs in session scratchpad
`top10_07_2216/`: `validateA.gp`, `gen_twisted.c`, `gen_twisted2.c`,
`halve8_twisted.m`, `stage2_B1000.m`, `stage3_untwisted.m`, `goodred16.m`,
`fiber_descent.gp`, `nearmiss_tors.m`, `simple_survivor.m`, `*_out.txt`.
Total ~8 CPU-min, single-threaded. No repo file modified.)*

## 1. Dossier and validation (all reproduced)

Target `(2,2,16)`, order 64. Chart 1 = halving the marked order-8 class on
`M_1(8,2,2)` (norm surface `Y`, 0/69 field passes, suspected Brauer-Manin on the
2-cover `Z`); chart 2 = order-8 second halving on the newly discovered *twisted*
`(2,2,2,8)` family (`claude_tier1_item1_2248_descent.md` sec. 4), never tried.

- **Validation A** (`validateA.gp`, instant): near-miss `(u,v)=(-49/15,25/24)`:
  `alpha,beta,gamma,delta` all squares, common class `d=-1`, field condition
  `d(s2-(r-+1)^2)` fails for BOTH signs, and `t = 49/40` (the S_3-orbit claim).
  Exactly as stated in `claude_tier2_2216_halving_sweep.md` /
  `claude_next_2216_normsurface.md`.
- **Validation B**: regenerated the twisted family
  `{abcd, a(a+b)(c-a)(d-a), b(a+b)(c-b)(d-b), c(c+d)(c-a)(c-b), d(c+d)(d-a)(d-b)}`
  all squares: exactly **52 tuples with d<=200** (dossier's count), of which
  **25 are gcd-primitive** — scalings `(ka,kb,kc,kd)` give Q-isomorphic curves
  (`x -> k^2 x` pulls out `k^10`, a square), so 25 distinct curves.

## 2. Main test A (chart 2): first-ever order-8 halving test — complete negative, with structure

Method (`halve8_twisted.m`): odd quintic `y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)`,
x-T map injective on `J(Q)/2J(Q)`; for every element `g` of `TorsionSubgroup`
of order 8 (32 per curve = all 16 twist classes `E+t`, `t in J[2]`, times +-),
compute the exact delta vector of its Mumford `U` (with the `f'(e)` correction,
conventions identical to the 0/600-mismatch-validated `claude_twist_sweep.py`);
`delta` trivial <=> rational order-16 point <=> torsion >= `(2,2,2,16) > (2,2,16)`.

Results on the 25 primitive curves (0.4 s calibration on 2 tuples, 2.5 s full):

- **all 25 have torsion exactly `[2,2,2,8]`** (extends the dossier's 8-curve
  verification to the whole d<=200 family);
- **0 halving hits; best delta = 1/5 square components** over all 25 x 32
  classes — never even 2/5 (cf. 3/5 attained in the tor2244 twist sweep);
- independent **local sieve**: `(2,2,2,16) <= J(F_p)` (i.e. 16 | largest
  invariant factor; four even factors are automatic with split `f`) FAILS at
  some good prime `p <= 89` for every curve — two independent methods agree.

**Forced-bad-reduction profile** (`goodred16.m`, enumerating ALL 6-rational-
Weierstrass-point curves over `F_p`, both twists, as in `goodred_profile.m`):
`(2,2,2,16)` is **NOT forced-bad at any p in {11,...,31}** — e.g. 60/168 curves
over `F_11` admit it (invariants `(2,2,2,16)` realized on the nose). So unlike
`(2,2,2,12)` there is no new CRT constraint, the 25 local exclusions are
curve-specific, and the family is not structurally dead. Measured conditional
pass rates P(16|inv_max given 8|inv_max): p=11: 60/60=1.00, 13: 60/70, 17:
120/170, 19: 240/380, 23: 780/1020, 29: 540/1380, 31: 900/2100 — a usable
density model for sizing the production search.

**Scale-up to d<=1000** (`gen_twisted2.c`: for fixed `a<b<c`, `abcd` square
forces `d = core(abc)*m^2` with `core(abc)` computed by pairwise-gcd core
composition — O(B^3) triples instead of O(B^4); 4 s single-threaded; B=200
output byte-identical to the naive sieve): **154 primitives**. Staged test
(`stage2_B1000.m`, 1.8 s): local-16 sieve kills **153/154** (first-excluding-
prime histogram spread 13..157, long tail); the unique survivor
`(318,371,396,462)` passes local-16 at ALL good p<200 yet delta-fails (best
1/5) — and is **geometrically split**: fails the 12th-power certificate at 40
primes, L-poly factorization shape `[2,2]` at 21/25 and `[2]` at 4/25 good
primes. So every geometrically simple member to d<=1000 is locally excluded
below 200.

**Route-3 bonus** (`stage3_untwisted.m`, 7 s): same staged test on the
UNTWISTED `(2,2,2,8)` family (`data/tor2228.txt`, 619 tuples, d<=16384) — also
never halving-tested: 4 local survivors `(14,116,259,2146)`,
`(31,124,836,3344)`, `(89,161,6319,11431)`, `(206,444,1545,3330)`, **all
delta-fail at 1/5, 0 hits**.

Net: ~800 curves across both known `(2,2,2,8)` components, all 16 twists each:
**zero order-8 halvings, and the best delta never exceeds 1/5**. That uniform
1/5 ceiling (against 3/5 occurring freely in the analogous tor2244 sweep) is
the single most suggestive new datum: it smells like a parity/Brauer-type
obstruction on order-8 halving in the full-2-torsion odd model, and should be
checked symbolically before any big enumeration is paid for.

## 3. Main test B (chart 1): the three d=-2 fibers have rank 2 — Y(Q) is infinite, and still 0 passes

`fiber_descent.gp`: fiber of `Y` at fixed `a` (chart `u=1/(1+da^2)` etc.,
`d=-2`): parametrize the conic `g^2 = (3+2dA) + (2d+d^2A)b^2` through the known
point, substitute into `h^2 = -d(4+3dA) - d(3d+2d^2A)b^2` -> genus-1 quartic
`F(x)`; `ellfromeqn` + `ellrank` (PARI 2.18):

| fiber | quartic Jacobian | rank (lower=upper, unconditional) | torsion |
|---|---|---|---|
| a=5/4  | `2809/324 x^4 - ...` | **2** | `[4,2]` |
| a=11/9 | `2809/324 x^4 - ...` | **2** | `[4,2]` |
| a=17/12| `2209/144 x^4 - ...` | **2** | `[4,2]` |

So **Y(Q) is provably infinite** (new — the prior state was 69 sporadic points).
Point harvest (`hyperellratpoints` to height 3e5, 15 s): 96+42+36 quartic
points -> 14 distinct open Y-points (dedupe by unordered `{u,v,t}`), field
condition `d(s2-(r-+1)^2)` square: **0/14** — consistent with and extending the
0/69.

**Twist gap found and closed at the near-miss locus**: all prior chart-1 sweeps
tested only the marked class `Q`; `(2,2,16)` needs ANY of the 16 order-8
classes `+-Q+t` divisible, and a rational half is automatically order-16
torsion — so exact `TorsionSubgroup` decides all twists at once.
`nearmiss_tors.m` (1.2 s) on all 5 near-miss curves (which are exactly the 5
least-obstructed Y-curves, e.g. `(-8/17,-81/161)` = chart `(-2; 5/4,11/9,53/12)`):
**torsion exactly `[2,2,8]` in every case** — no twisted halving either. The
twisted classes define 7 OTHER norm surfaces `Y_t` (delta(Q)*delta(t)
conditions) that have never been enumerated; that is now the freshest untouched
locus on chart 1.

## 4. Strategy (ranked)

1. **Symbolic parity audit, then chart-2 scale-up.** First (cheap, ~1 session):
   compute `delta(E+t)` symbolically on the twisted family — `E` the half of
   `D0+T_AB` — and test whether the observed 1/5 ceiling is a theorem (a
   product/norm relation forcing >=4 nonsquare components for every twist). If
   yes: route provably empty, write it up, stop. If no: run `gen_twisted2` to
   `B = 10^4` (O(B^3): ~1-2 core-h, embarrassingly parallel), staged local-16
   sieve to p<500, delta on survivors. Expected primitives ~10^3; the measured
   local pass rates say survivors will be rare and precious.
2. **Chart-1 twisted norm surfaces.** Derive the x-T criteria for the 7
   nonzero twists `Q+t` on `M_1(8,2,2)` (delta(t) explicit from the four
   rational roots {1,u,v,t}), reduce mod the S_3 symmetry (~2-3 genuinely new
   condition sets), conic-parametrize each twisted `Y_t` as in
   `claude_next_2216_normsurface.md`, and sweep (sweep2.c-style, minutes).
   Note: the `M_1(8,2,2)` configuration (four rational Weierstrass points + one
   quadratic pair <=> rational 2-torsion `(2,2,2)`) is the ONLY one that can
   give EXACT `(2,2,16)`; chart 2 can only overshoot to `(2,2,2,16)`.
3. **Chart-1 Brauer-Manin certification on the rank-2 fibers.** The fiber-
   restricted cover is now completely explicit: `Z_f: {h^2 = F(m), w^2 = X+(m)}`
   with both `F` and `X+ = d(s2-(r-1)^2)` rational functions of the conic
   parameter `m` (r, s2 rational on Y). Compute local invariants of the
   quaternion class `(X+, d)` along each fiber and sum against the rank-2
   Mordell-Weil group (generators known from `ellrank`); either certify the
   obstruction (fibers provably empty — the main chart is then morally dead) or
   locate where the invariant vanishes and aim the point search there. Also:
   group-law point generation (map `E -> quartic` inverse) to push 0/14 to
   0/10^3 cheaply if more evidence is wanted.

## 5. Verdict

Still no realization and no positive signal, but the picture sharpened on both
charts: chart 2 is pipeline-complete, cheap to scale, NOT locally obstructed as
a family (new fact), yet uniformly delta-dead so far (new fact: the 1/5
ceiling, ~800 curves, all twists — audit before scaling); chart 1 now has
provably infinitely many Y-points (rank-2 fibers, new fact) with the field
condition failing on every point ever tested, and the near-miss cluster is
dead in ALL 16 twists (new fact) — the Brauer-Manin hypothesis is strengthened
and is now concretely computable fiber-by-fiber. Highest-information next
moves: the symbolic parity audit (could kill chart 2 in one stroke) and the 7
twisted norm surfaces (the only untouched locus that can yield exact
`(2,2,16)`).

## 6. Exact commands

```text
gp -q validateA.gp                                   # dossier near-miss reproduction
gcc -O2 -o gen_twisted gen_twisted.c -lm && ./gen_twisted 200    # 52 tuples, 25 primitive
gcc -O2 -o gen_twisted2 gen_twisted2.c -lm && ./gen_twisted2 1000  # 154 primitives, 4 s
magma -b halve8_twisted.m        # 25 curves x 32 order-8 classes: 0 hits (2.5 s)
magma -b goodred16.m             # (2,2,2,16) not forced-bad at 11..31 (~2 min)
magma -b stage2_B1000.m          # local-16 sieve + delta on survivors (1.8 s)
magma -b stage3_untwisted.m      # same on data/tor2228.txt (7.4 s)
gp -q fiber_descent.gp           # ellrank=2 x3 fibers, 0/14 field passes (15 s at H=3e5)
magma -b nearmiss_tors.m         # all 5 near-misses: torsion exactly [2,2,8] (1.2 s)
magma -b simple_survivor.m       # (318,371,396,462) fails 12th-power cert; split
```
