# Assessment of the GPT-5.6-Sol `[2,22]` lanes (fork PR #9), with Lane-A measurements

Date: 2026-07-24.  Input: the "Strategic reset" section of
`notes/order222_from_order11.md` + `code/order222_three_root_sieve.m`
(cherry-picked from the fork, commit a6791c8).  Cross-checked against
`notes/claude_generic_222_2214_plan_2026_07_23.md` (the current [2,22]
program of record: routes B1/B2', Flynn lock, extended-DB census).

## Verdict per lane

### Lane C (mine the database) — ALREADY CLOSED

Addendum 2 of the plan note records the full extended-DB census
(6,216,959 curves): torsion `[11]`x171, `[22]`x37, `[33]`x1,
`[2,22]`x1 (= 19044.h.2, RM(sqrt5)); NO `[44]`, `[2,2,11]`, `[66]`.
Any exact group containing `[2,22]` would appear in that census; only the
known RM witness does.  Nothing further to mine.

### Lane B (deform the RM seed transversally to Humbert) — GOOD IDEA, WRONG CHART

The note instructs putting seed 19044.h.2 "in the same three-root chart as
Lane A" (`x(x-1)(x-r)(cubic)`, type `[1,1,1,3]`).  **That is impossible**:
verified fresh this session (`scratchpad/validate_222_funnel.m`), the
completed square of 19044.h.2 has factor type `[1,1,2,2]`, and so does the
second known `[2,22]` curve (corrected BLP C4) — each has only TWO rational
Weierstrass x-coordinates.  The deformation must therefore run in the
`x(x-1)*q1(x)*q2(x)` chart (both quadratics irreducible), which is exactly
the target shape of the existing routes B1 (CF-backward on
V11 ∩ [1,1,2,2], anchored at 19044.h.2) and B2' (BLP resultant surface's
[1,1,2,2]-locus, anchored at C4corr).

What Lane B genuinely ADDS to B1/B2' is machinery, not a new locus: impose
`[11](u,v) = 0` by denominator-cleared Cantor arithmetic on the chart with
Mumford variables, compute the Jacobian matrix at the RM seed mod several
good primes, look for tangent directions transverse to the seed's Humbert
equation, Hensel-lift, and rationally reconstruct — verifying every
reconstructed point exactly.  A transversality certificate at the seed would
prove the level scheme is not locally trapped in Humbert-5, which is
precisely the question B1 needs answered.  This is the highest-value item in
the PR; implement it on the `[1,1,2,2]` chart.

### Lane A (three-root integer-box sieve) — the RECORDED LOSING PATTERN, run as a cheap measurement only

The chart builds in the cheap condition (2-rank 2) and sieves for the
expensive one (11-torsion, codim 2) — exactly the "blind 2-rank-gated scan"
recorded as hopeless at commit 34cfb59 ("odd torsion is codim 2"; the
project's winning pattern is the reverse: build the odd part algebraically,
parametrize the 2-condition).  Also note both known `[2,22]` curves live in
the OTHER 2-rank-2 stratum (`[1,1,2,2]`, odd class at infinity, CF order
11); the `[1,1,1,3]` chart's known occupant is the generic `[2,14]` witness
1416.b.  The PR's own stop/go rule only asks for survivor-rate measurements
at two box sizes, which is fair and cheap — done below.

## Lane-A execution record (this session)

Funnel validation: the `ReductionGCD44` filter passes both known `[2,22]`
curves with gcd exactly 44 (sharp), and `TorsionSubgroup` returns `[2,22]`
for both.  Pipeline sound.

Script fixes (house conventions + calibration), committed:

- added `SetColumns(0)`, `SetMemoryLimit` (MemGB param), `PROGRESS`/
  `SEARCH_DONE` markers, and the missing `quit;` (without it a `-b` run
  never exits);
- extended the prime list `[3..19]` -> `[3..41]`: the B=6 run's single
  false positive `(r,a,b,c) = (6,2,2,2)` survived on a 3-prime gcd (88)
  because its disc killed 4 of the 7 original primes — backup primes close
  that hole at negligible cost (they only engage for near-survivors);
- `rlo`/`rhi` shard parameters for parallel boxes.

Measurements:

- **B=6** (`results/order222_three_root_sieve_B6.log`, 6m41s, 7-prime
  version): 18,744 curves tested, 1 gcd44 survivor = the false positive
  above (exact torsion `[2,2]`).  Zero true candidates.
- **B=12** (12-prime version, 3 r-shards,
  `results/order222_three_root_sieve_B12_shard{1,2,3}.log`): ~290k curves.
  RESULT TO BE APPENDED BELOW WHEN THE SHARDS FINISH.

Stop/go recommendation: unless B=12 (or a possible B=24) produces an exact
`[2,22]` or an anomalous cluster of gcd44 survivors with genuine 11-part,
Lane A stops at the measured rates and compute goes to B1/B2' plus the
Lane-B transversality machinery on the `[1,1,2,2]` chart.

## Lane-B anchor data exported (this session)

`code/claude_222_laneB_seed_export.m` (+ `results/claude_222_laneB_seeds.log`)
normalizes both RM witnesses into the `[1,1,2,2]` chart.  **Both land in the
MONIC chart** `y^2 = x(x-1)(x^2+p1*x+p2)(x^2+p3*x+p4)` — the leading unit
after normalization is a perfect square in both cases (1/16 and 27000^2):

- 19044.h.2:  `(p1,p2,p3,p4) = (-3, 8, 4, 27)`; order-22 generator is a
  rational-point class, Mumford `u = x-3` (P = (3,-192) on the 16-scaled
  model); 11-part `u11 = x^2-3x, v11 = 64x`.
- BLP C4corr: `(p1,p2,p3,p4) = (-34/15, 32/45, 31/15, -2/9)`; order-22
  generator `u22 = x^2+(11/3)x-2/9`; 11-part `u11 = x^2-(3/5)x+2/15`.

So the Lane-B ambient is A^4 with coordinates `(p1,p2,p3,p4)`, unknowns
extended by Mumford `(u1,u0,v1,v0)`, conditions `u | f - v^2` plus
`[11](u,v) = 0` by denominator-cleared Cantor arithmetic.  The incidence
variety W = {(p,D) : 11D = 0} is 4-dimensional and generically finite
(120:1) over the chart, so near a smooth seed the p-adic sheet is a full
neighborhood and "transversality to Humbert" is about finding NEW rational
points of W near the seed (Hensel lift + simultaneous rational
reconstruction, verified exactly) — the same computational class as route
B1(b)'s ideal interpolation, now with two explicit rational anchors.

## B=12 result

Three r-shards ([-12..-5], [-4..4], [5..12]), 12-prime filter, ~35 min wall
each (3 parallel Magma jobs):

```text
shard1: tested 108864, gcd44 survivors 0, exact checks 0
shard2: tested  95256, gcd44 survivors 0, exact checks 0
shard3: tested 108864, gcd44 survivors 0, exact checks 0
TOTAL : tested 312984, gcd44 survivors 0
```

Measured survivor rates: B=6 (7-prime filter): 1/18744, the FP documented
above; B=12 (12-prime filter): 0/312984 — the extended prime list eliminates
even the disc-degenerate FP channel.  Zero true candidates at either box.

**Lane A concluded per the stop/go rule.**  Two box sizes measured, survivor
rate at the noise floor, no anomaly to chase; consistent with the codim-2
thinness that killed the 34cfb59 scans.  A B=24 box (~5M curves, ~9 CPU-h)
would raise the expected FP count to O(1) without changing the structural
picture and is NOT recommended.  Compute goes to B1/B2' and the Lane-B
transversal machinery on the monic [1,1,2,2] chart (anchor data above).
