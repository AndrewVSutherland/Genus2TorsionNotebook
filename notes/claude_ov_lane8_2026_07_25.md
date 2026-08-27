# Lane 8 — the HLP split-torsion landscape, the deformation, and Howe's f70

Date: 2026-07-25. Overnight campaign, lane 8. Three sessions (two killed by an
API outage); this note is written from the committed logs of all three.

Commits: `db8137a` (f70 is a (3,3)-glue; corrected+extended landscape census),
`5ada083` (dimension of the HLP split-torsion locus, target by target),
`5ba80da` (the deformation experiment: lines), `fd285a5` (Z/63 sieve to height
1000), `08b01d4` (Z/63 sieve to height 2000), `f8ac77d` (the deformation
experiment: arcs), and this note.

---

## 0. What the lane was for

The Howe–Leprévost–Poonen gluing construction reaches large rational torsion
groups on genus-2 Jacobians — Z/35, Z/45, Z/63, Z/70, Z/7×Z/7, Z/5×Z/10,
Z/2×Z/24, Z/3×Z/12 — but *only on split Jacobians*, which are worthless for this
project (a result requires a **geometrically simple** J). Lane 8 asked three
things:

1. **f70.** Howe's order-70 curve: is it geometrically split, and therefore does
   Z/35 fall out of it or not?
2. **The deformation.** Take a split HLP curve with the target torsion, and
   deform it *off* the split locus keeping the torsion. Does the torsion locus
   have a component **transverse** to the split locus?
3. **Z/63**, which would be a new **cyclic record** (current record Z/40).

Plus: interpret the landscape census — does the near-absence of absolutely
simple candidates mean these orders are systematically non-simple, or is it a
chart artifact?

All four are answered below. Two of the answers are negative and one of the
negatives is structural (it closes the deformation route, not just one attempt).

---

## 1. f70 — CLOSED. The answer is **NO**.

**The curve.** `f70 = 22x^5 + (697/144)x^4 - (645/4)x^3 + (1045/4)x^2 - 162x + 36`,
integral model `144·f70 = 3168x^5 + 697x^4 - 23220x^3 + 37620x^2 - 23328x + 5184`.
Exact Magma `TorsionSubgroup` = `[70]`, order 70; 2-torsion `[2]`; factor type
`[2,3]`. (`results/claude_ov_lane8_f70.log`, `..._verify.log`.)

**It is a (3,3)-glue.** J(f70) is isogenous over Q to E1 × E2 with

| | Weierstrass coefficients | conductor | torsion |
|---|---|---|---|
| E1 | `[1,0,0,-45,81]` | 66 | Z/10 |
| E2 | `[1,0,0,-5774401,5346023177]` | 858 | Z/7 |

- `#J(F_p) = #E1(F_p)·#E2(F_p)` at **all 299 good p < 2000**, zero disagreements
  (`results/claude_ov_lane8_f70split.log`).
- Upgraded in this session from the single value `#J(F_p)` to the **full local
  L-factor identity** `chi_J(T) = (T² − a_p(E1)T + p)(T² − a_p(E2)T + p)`.
  **Established: agreement at all 1227 good `p < 10009`, zero disagreements**
  (`results/claude_ov_lane8_f70close.log`). This is the complete local statement
  of `J ~ E1 × E2` at each of those primes, not just an order match, and it is
  4× the prime range of the overnight `#J` check. The job was launched with
  `PMAX = 20000` and **stopped by PID at p = 10009** once the result was
  established, rather than left as an orphan process; the conclusion does not
  depend on the extension.
- **Not (2,2):** E1's 2-division polynomial factors 1+2, E2's is an irreducible
  cubic, so E1[2] ≇ E2[2]. **Is (3,3):** `a_p(E1) ≡ a_p(E2) mod 3` for every good
  p tested and for **no other** ℓ ≤ 37; independently, `psi_3(E1)` and `psi_3(E2)`
  are both irreducible quartics whose splitting fields both have degree 24 and
  are **isomorphic** (Magma `IsIsomorphic`), i.e. Q(E1[3]) = Q(E2[3])
  (`results/claude_ov_lane8_f70subfield.log`).
- Corroborating: chi_p is reducible over Z at **all 91 good p < 500**, with
  factor-degree type `[2]` or `[2,2]` only — never irreducible
  (`results/claude_ov_lane8_verify.log`).

**Consequence, stated plainly.** Howe's order-70 curve is **not geometrically
simple**: its Jacobian is a (3,3)-gluing of two elliptic curves. `J(Q)_tors =
Z/70` does contain a Z/35 — it is the index-2 subgroup `2·J(Q)_tors` — but that
Z/35 lives on the *same split Jacobian*. **Z/35 does not fall out of f70.**

**Extra:** the deformation test of §2 was also run *on f70 itself*, for both
N = 70 and N = 35: 5,000 random line directions and 80,000 polynomial arcs of
degree 2 and 3, **all killed**. One cannot deform Howe's curve off its split
locus along any tested line or arc either.

**What is not proved.** The isogeny is certified by matching L-factors at
thousands of primes, not by an unconditional Faltings–Serre computation and not
by an explicit degree-3 map C → E1. I tried to get an unconditional certificate:
Magma's `Subfields` does not accept `FldFunFracSch`, and on
`AlgorithmicFunctionField(F)` it computes only subfields of the degree-2
extension `A/Q(x)` (it returns the trivial one); no endomorphisms package
(Costa–Mascot–Sijsling–Voight) is installed on this box. Given the three
independent strands above plus f70's provenance in the *split-Jacobian*
literature (HLP, Forum Math. 12 (2000); Howe, BLMS 47 (2015), Thm 2.1 is
explicitly a (3,3)-gluing), the yes/no is settled.

---

## 2. The deformation — the lane's real target

### 2a. First-order deformation theory is **vacuous** here. (Structural result.)

In characteristic 0 the group scheme `J[N]` is **finite étale** over the base.
A rational point of order N on the anchor therefore lifts **uniquely** over every
Artinian thickening and over the completed local ring of moduli at the anchor.
Consequences:

- There is **no local obstruction** to deforming off the split locus. The
  tangent/normal computation carries no information.
- The "torsion locus" is **not a subvariety of M_2 at all**. It is the image of
  the rational points of the finite étale cover
  `W_N = {(C,P) : P ∈ J(C), ord P = N} → M_2`,
  every component of which is 3-dimensional and unramified over M_2. So the
  split locus preimage `π^{-1}(L) ∩ W_N` is a *divisor* in a 3-fold, and
  "transverse component" can only mean an **arithmetic** statement about
  W_N(Q) off that divisor.

This reproduces, for these anchors, exactly what
`notes/m612_hlp_deformation.md` found for [6,12] in July: the 21×28 marked
incidence Jacobian has rank 21, kernel 7 = 3 moduli + 4 group directions, i.e.
the marked cover is étale over M_2 at the HLP seed. **The deformation idea
therefore does not reduce the problem — it is equivalent to the original
problem.** That is the honest headline of this lane.

### 2b. What *can* be decided, rigorously and cheaply

Let `C_t : y² = f0(x) + t·g(x)`. If the order-N section extended over `Q(t)` —
i.e. if the line lay inside the torsion locus **as a family** — it would
specialise at *every* `t0 ∈ F_p` of good reduction, so `N | #J_{t0}(F_p)` for all
such `t0`. Therefore:

> **One residue `t0 ∈ F_p` with good reduction and `N ∤ #J_{t0}(F_p)` PROVES that
> the line `f0 + t·g` carries no Q(t)-rational order-N family.**

**Positive control** (`results/claude_ov_lane8_transverse_control.log`): on the
contact-7 and contact-9 charts, which *do* carry a marked rational class,

| p | contact-7: 7\|#J | contact-9: 9\|#J |
|---|---|---|
| 11 | 103/103 = 1.0000 | 9/9 = 1.0000 |
| 13 | 148/148 = 1.0000 | 11/11 = 1.0000 |
| 17 | 257/257 = 1.0000 | 15/15 = 1.0000 |
| 19 | 327/327 = 1.0000 | 16/16 = 1.0000 |
| 23 | 485/485 = 1.0000 | 20/20 = 1.0000 |
| 29 | 790/790 = 1.0000 | 24/24 = 1.0000 |

A line inside the torsion locus gives density **exactly 1**. The test has power.

**The split divisor through the anchor.** Every recovered HLP model is an even
sextic `f0 = a x⁶ + b x⁴ + c x² + d`, i.e. bielliptic = (2,2)-split. From the
sl2 + scaling action on binary sextics, the tangent space to the bielliptic
locus at f0 is `{even sextics} + span(u1,u2)` with

```
u1 = f0'              = (2c, 4b, 6a)     in odd coordinates (x, x^3, x^5)
u2 = 6x·f0 − x²·f0'   = (6d, 4c, 2b)
```

which is 6-dimensional in the 7-dimensional coefficient space: **codimension 1**,
with explicit normal functional `u1 × u2` (printed per anchor in the logs). A
direction g is **transverse** iff `oddpart(g) · (u1 × u2) ≠ 0`.

### 2c. The measurement (all good p < 300)

`code/claude_ov_lane8_transverse.gp`, one process per anchor. (The predecessor's
`code/claude_ov_lane8_transverse.m` is the same experiment in Magma; it is
committed for the record but was ~100× too slow — 47–140 ms per L-polynomial
versus 0.3–4 ms in PARI — so the run was done in GP.)

| anchor | N | transverse `f0+t·x^{1,3,5}` | in-locus `f0+t·x^{0,2,4}` | baseline (random sextics) | in-locus / baseline |
|---|---|---|---|---|---|
| Z/5×Z/10 | 50 | **0.01139** | 0.05026 | 0.01247 | 4.0 |
| Z/7×Z/7  | 49 | **0.00855** | 0.03152 | 0.00668 | 4.7 |
| Z/63     | 63 | **0.02667** | 0.10825 | 0.02546 | 4.3 |
| Z/45     | 45 | **0.03745** | 0.14874 | 0.03560 | 4.2 |
| Z/35     | 35 | **0.03661** | 0.12013 | 0.03611 | 3.3 |
| Z/2×Z/24 | 48 | **0.04214** | 0.26911 | 0.04544 | 5.9 |
| Z/3×Z/12 | 36 | **0.04037** | 0.18373 | 0.03781 | 4.9 |
| f70 (N=70) | 70 | 0.02892 | 0.02936 | 0.02071 | 1.0 |
| f70 (N=35) | 35 | 0.03574 | 0.03727 | 0.03602 | 1.0 |

(24k–25k line points and ~33k baseline points per anchor; the anchor itself
passes the compatibility test at *all* 55–59 good p < 300 in every case, which is
the sanity check that the pipeline is looking at the right curve.)

**Reading.** The transverse density equals the baseline to within noise at every
anchor: *the anchor exerts no influence in the transverse direction.* The
in-locus lines sit 3–6× above baseline, which is exactly the (2,2)-split
point-count statistics (`#J = #E1·#E2` makes `N | #J` much more likely for
composite N) — so the measurement does resolve a real structural effect when
there is one. For f70, whose model is a quintic with no even structure, both
line families are generic and both equal baseline, as they should.

### 2d. Direction kill — 156,695 lines, **zero survivors**

For each bielliptic anchor: **every** primitive sign-normalised odd direction
`c1·x + c3·x³ + c5·x⁵` with `|ci| ≤ 15` that is transverse to the bielliptic
locus — **12,385** of them, an exhaustive box — plus **5,000** random directions
in the full 7-dimensional coefficient space with coefficients in [−20,20].
For each f70 target: 5,000 random directions. Totals:

- 7 bielliptic anchors × 17,385 + 2 f70 targets × 5,000 = **156,695 lines**
- **156,695 killed, 0 survivors.**
- The hardest line needed only `p ≤ 13` and `t ≤ 6`. Most die at `p ∈ {3,5,7}`,
  `t = 1`.

### 2e. Arcs, not just lines — 639,979 more families, still zero survivors

The criterion of §2b applies verbatim to **any** polynomially parametrised
family `f_t = f0 + t·g1 + t²·g2 + … + t^d·g_d`, so
`code/claude_ov_lane8_arcs.gp` extends the test from lines to arcs. Two classes:

- **(A) transverse at first order** — `g1` outside the tangent space to the
  bielliptic locus. **359,986** arcs (9 anchor/target pairs × degrees 2 and 3 ×
  20,000 each).
- **(B) tangent-then-leaving** — `g1` drawn from the *full* tangent space to the
  bielliptic locus (even part arbitrary, odd part in span(u1,u2)), with a later
  term transverse. These arcs are tangent to the split locus and leave it only
  at order ≥ 2, so **a line test cannot see them at all**. **279,993** arcs
  (7 bielliptic anchors × degrees 2 and 3 × 20,000 each).

**All 639,979 killed, zero survivors**, hardest residue needed anywhere
`p ≤ 17`, `t ≤ 7`. Two controls are printed per run and pass everywhere: the
anchor itself passes at every good `p < 200`, and the **constant** arc
`f_t = f0` is *not* killed — so the killer is not trivially killing everything.

**Running total for the lane: 156,695 lines + 639,979 arcs = 796,674
Q(t)-families through the eight anchors, none of which carries the torsion.**

### 2f. Does the fibre genus/rank data already answer it?

**No — it answers a different question, and I want to be explicit about that.**
`5ada083` computed the genus and rank of the fibres of the *split-torsion locus*
`S_{M,N} : w² = Δ_M(s)·Δ_N(u)`. That decides whether the **split** curves with the
target torsion are Zariski-dense in the bielliptic surface. It says nothing about
transversality, because it never leaves the bielliptic locus. The transversality
question needs §2a (étale ⇒ vacuous first order) plus §2b–2d (arithmetic line
test). The two results are complementary, not substitutes.

### 2g. Verdict on the deformation

- Geometrically there *is* a transverse direction — trivially, because the cover
  is étale, so this is not information.
- Arithmetically, **no line and no low-degree arc through any of the eight
  anchors carries the torsion**: 156,695 lines (including an exhaustive box)
  and 639,979 arcs of degree 2 and 3, of which 279,993 are tangent to the split
  locus and therefore invisible to any line test. 796,674 families, 0
  survivors.
- This matches the only previous deep probe in the repo: for [6,12],
  `notes/m612_hlp_deformation.md` found that on the best transverse slice the
  marked order-3 support cover is irreducible of degree 40 and **genus 51**
  (Picard–Lefschetz + Riemann–Hurwitz, with an irreducible degree-40 resolvent
  over F_101(t) certifying irreducibility in characteristic 0), and the order-4
  support cover has degree 120 and genus 181 if connected. Faltings then gives
  finitely many rational t on the slice, and a height-10^6 sieve on the `1+x`
  slice left only `t = 0`, the original split point.
- **Conclusion: the "deform the HLP curve off the split locus" route is closed
  for linear slices, and it is closed for a structural reason (the cover is
  étale, so nothing local can help; the covers restricted to a slice are curves
  of large genus).** A future attempt must use a two-parameter transverse
  *surface* on which one contact block can acquire a section, not a line.

---

## 3. The split-torsion locus, target by target (recap of `5ada083`)

The HLP gluing locus for the pair (M,N) is carried by `S_{M,N} : w² =
Δ_M(s)·Δ_N(u)`; fixing s gives the fibre `w² = core(Δ_M(s))·Δ_N(u)`.

Squarefree-core degrees of `Δ_N(t)` for `E_N^t`:

| N | 5 | 6 | 7 | 9 | 10 | 12 | 26 | 28 |
|---|---|---|---|---|---|---|---|---|
| deg core Δ_N | 3 | 2 | 5 | 7 | 3 | 4 | 0 | 0 |

| target | (M,N) | fibre genus | measured | split-torsion locus over Q |
|---|---|---|---|---|
| Z/2×Z/24 | (26,28) | — | gluing condition **vacuous** | full 2-dim (s,u) plane |
| Z/3×Z/12 | (12,6) | 0 | 146–348 small points per fibre, 7/7 s | 2-dim, dense |
| Z/5×Z/10 | (10,10) | 1 | rank ≥ 1 in **15/15** tested s | 2-dim, dense |
| Z/35 | (7,5) | 1 | rank ≥ 1 in **10/13** | 2-dim, dense |
| Z/45 | (9,5) | 1 | rank ≥ 1 in **6/10** | 2-dim, dense |
| Z/7×Z/7 | (7,7) | 2 | Faltings | **finitely many** partners u per s |
| Z/63 | (7,9) | 3 | Faltings | **finitely many** partners u per s |

Also measured there: `E_N^t[2]` is reducible (1+2) exactly for N even (6,10,12)
and split (1+1+1) for N = 26,28, so for those the 2-torsion match is *equivalent*
to the discriminant match; for N odd (5,7,9) it is strictly stronger — measured
3/12 for Z/63 in a height-6 box, 12/12 for Z/35, 20/20 for Z/5×Z/10, 30/30 for
Z/3×Z/12, 0/0 for Z/45 and Z/7×Z/7.

All seven anchors were rebuilt and verified exactly
(`results/claude_ov_lane8_verify.log`): exact `TorsionSubgroup` equal to the
target in every case, both elliptic subcovers identified with conductor and
torsion, `#J(F_p) = #E1·#E2` at 88–92 good p < 500, and chi reducible at every
one of them — i.e. all seven are certified **split**, as expected.

---

## 4. Z/63 — the record target

Z/63 would be a new **cyclic** record (current record Z/40, LMFDB
`g2c` label with `geom_end_alg = Q`).

**What exists:** a split Z/63 curve, HLP Forum Math. 12 (2000) eq. (4),
`y² = 897x⁶ − 197570x⁴ + 79136353x² − 146398496`. Exact torsion `[63]`;
subcovers E1 = conductor 1482, torsion Z/9, rank 1 and E2 = conductor 42978,
torsion Z/7, rank 1; `#J = #E1·#E2` at 89 good primes; chi reducible at all 56
good p < 300. Split, so useless as a result — but it is a valid anchor.

**What was searched, and how far:**

| chart | parameters | box | points | primes | survivors | outcome |
|---|---|---|---|---|---|---|
| contact-7 | a, b ∈ Q | \|n\|,d ≤ 100 (overnight) | 1.48×10⁸ | 12 | 0 | — |
| contact-7 | a, b ∈ Q | \|n\|,d ≤ **1000** (this session) | **1.480×10¹²** | **16** | 1 | **false positive**: a = 628/779, b = 994/671 has exact torsion `[7]` |
| contact-7 | a, b ∈ Q | \|n\|,d ≤ **2000** (this session) | **2.368×10¹³** | **16** | 7 | **all 7 false positives**, exact torsion `[7]` in every case |
| contact-9 | a ∈ Q | \|n\| ≤ 3×10⁴, d ≤ 3×10³ (overnight) | 1.09×10⁸ | 12 | 0 | — |
| contact-9 | a ∈ Q | \|n\| ≤ 3×10⁵, d ≤ 3×10⁴ (this session) | **1.0943×10¹⁰** | 12 | 8 | **all 8 false positives**: exact torsion `[9]` in every case |

The eight contact-9 survivors were
`a = 201289/2665, −259396/8151, −30283/18088, 31882/20213, −61952/20655,
−186209/23001, −5155/24682, −132715/29497`; each was rebuilt on an integral
model and given an exact Magma `TorsionSubgroup`, which returned `[9]` (order 9,
equal to the marked class order) for all eight
(`results/claude_ov_lane8_c63verify_c9all.log`). Expected chance survivors at 12
primes and per-prime density ≈ 0.2 is `1.09×10¹⁰ × 0.2¹² ≈ 45`; observing 8 is in
the right ballpark.

The contact-7 extension is a **10⁴× enlargement** of the overnight box. The
rewritten sieve (`code/claude_ov_lane8_c63sieve.c`: precomputed coprime-pair
residue tables, staged prime-by-prime filtering, pthreads) reproduces the
overnight result exactly (12175 pairs → 148,230,625 points → 0 survivors) in
0.77 s, then covers 1.48×10¹² points in ~2 minutes on 24 threads. Expected
chance survivors at 16 primes and per-prime density ≈ 0.14 is
`1.48×10¹² × 0.14¹⁶ ≈ 0.02`, so the single survivor observed is consistent with
noise, and it verified as torsion exactly Z/7.

The seven contact-7 survivors at height 2000 are
`(a,b) = (354/377, 1931/1301), (319/527, −961/1064), (−1623/671, 1489/1073),
(628/779, 994/671), (−1878/901, 1391/899), (1833/989, −538/221),
(1462/1311, 1563/1859)`; **all seven** have exact Magma `TorsionSubgroup` `[7]`,
order 7, equal to the marked class order
(`results/claude_ov_lane8_c63verify_c7h2000.log`,
`..._c7h2000b.log`). Observing 7 rather than the naive expectation ~0.4 means the
true per-prime density is a little above 0.14 (the measured chart densities in
`results/claude_ov_lane8_c7_landscape.log` run 0.108–0.158), not that anything
survived.

**Result: no Z/63 curve on the contact-7 chart with both parameters of naive
height ≤ 2000, and none on the contact-9 chart with `a = n/d`, `|n| ≤ 3×10⁵`,
`d ≤ 3×10⁴`. Every one of the 16 survivors of the three sieves was rebuilt on an
integral model and killed by an exact Magma `TorsionSubgroup`.**

**No local obstruction.** The corrected census (§5) finds compatible
*absolutely simple* Weil classes for Z/63 at 59 of 61 primes p ≤ 293. The only
prime with **no** compatible class at all is p = 3 — so **every** genus-2 curve
over Q with a rational point of order 63 has **bad reduction at 3**. (Consistent:
the HLP Z/63 anchor has `3^15 ‖ disc`.)

**Transverse deformation from the split anchor:** 17,385 lines, all killed (§2d).

---

## 5. The landscape census — artifact, or systematically non-simple?

Both, and they must be separated.

**(a) The small-p "near-absence" WAS an artifact.** The overnight
`claude_ov_lane8_landscape.c/.gp` enumerated *curves* at p = 3,5,7,11,13 only
(O(p³) work), and used the wrong condition for the non-cyclic targets (it tested
`N | #J` where `(Z/n)² ≤ J(F_p)` actually needs n-rank ≥ 2, i.e.
`(T−1)² | chi mod n` together with `n² | chi(1)`). At p = 3 the whole Weil region
holds 63 isogeny classes and at p = 5 it holds 129 — "no absolutely simple
compatible class" there is a statement about a tiny box, not about the target.

The corrected census `code/claude_ov_lane8_census2.gp` enumerates **isogeny
classes** over the exact Weil region (O(p) classes per prime, so p reaches 293)
with the correct rank-2 condition. Result
(`results/claude_ov_lane8_census2.log`, 61 primes 3 ≤ p ≤ 293):

| target | primes with NO compatible class | primes that CAN certify simplicity |
|---|---|---|
| Z/35 | none | 60 / 61 (smallest 5) |
| Z/45 | {3} | 59 / 61 (smallest 7) |
| Z/63 | {3} | 59 / 61 (smallest 7) |
| Z/70 | {3} | 59 / 61 (smallest 7) |
| Z/7×Z/7 | none | 54 / 61 (smallest 23) |
| Z/5×Z/10 | {3} | 53 / 61 (smallest 11) |
| Z/2×Z/24 | {3} | 59 / 61 (smallest 7) |
| Z/3×Z/12 | none | 56 / 61 (smallest 13) |
| Z/40 (control, simple known) | {3} | 59 / 61 |
| Z/30 (control, simple known) | none | 59 / 61 |
| Z/48 (control) | {3} | 59 / 61 |

Compatible **and absolutely simple** classes exist at 53–60 of 61 primes for
*every* target, and the controls Z/40 and Z/30 — orders for which geometrically
simple curves are known — behave identically. **There is no local obstruction to
a geometrically simple curve for any of these orders.** (The Weil-admissibility
test is a *superset* of the isogeny classes of abelian surfaces, so "no
compatible class at p" is safe and "a compatible abs-simple class exists" is an
upper bound.)

The one genuine small-p fact: **p = 3 admits no Weil class compatible with Z/45,
Z/63, Z/70, Z/5×Z/10, Z/2×Z/24** (nor Z/40, Z/48), so every curve over Q with
those torsion groups has bad reduction at 3. Z/35, Z/7×Z/7, Z/3×Z/12 and Z/30
are not so constrained.

**(b) The *global* near-absence is REAL, not an artifact.** Queried against
LMFDB's extended/alpha genus-2 table `g2c_curves_new` (6,216,959 curves, torsion
orders up to 96), grouping by torsion subgroup and geometric endomorphism
algebra:

| torsion | # curves | `geom_end_alg` | geometrically simple |
|---|---|---|---|
| Z/35 | 6 | Q×Q (6) | **0** |
| Z/2×Z/24 | 76 | Q×Q (76) | **0** |
| Z/3×Z/12 | 18 | Q×Q (17), M₂(Q) (1) | **0** |
| Z/5×Z/10 | 1 | Q×Q (1) | **0** |
| Z/6×Z/6 | 176 | Q×Q (168), M₂(Q) (8) | **0** |
| Z/48 | 5 | Q×Q (5) | **0** |
| Z/45, Z/63, Z/70, Z/7×Z/7 | **0** | — | — (do not occur at all) |
| **Z/30 (control)** | 94 | Q×Q (88), **Q (6)** | **6** |
| **Z/40 (control)** | 10 | Q×Q (9), **Q (1)** | **1** |
| Z/39 | 2 | — | 2 |
| Z/2×Z/2×Z/14 | 2 | — | 2 |

So: **every** curve in a 6.2-million-curve database carrying one of lane 8's
torsion groups is geometrically split (`Q×Q` or `M₂(Q)`), while the controls
Z/30 and Z/40 do produce `geom_end_alg = Q`. In the smaller certified table
`g2c_curves` (|disc| < 10⁶, 66,158 curves) the largest torsion order occurring at
all is 39.

**Answer to the question as posed.** The near-absence of absolutely simple
candidates in the *local* census was a chart/small-p artifact and is now
corrected: there is no local obstruction. But the near-absence of absolutely
simple *curves* is not an artifact — it is a genuine, sharply measured global
phenomenon: for these orders every known example, in the literature and in a
6.2M-curve database, is split. Combined with §2, the picture is coherent: these
orders are reachable only through gluing constructions, which produce split
Jacobians by design, and the natural attempt to escape (deform off the split
locus) is obstructed arithmetically, not locally.

---

## 6. Numbers at a glance

- f70: torsion `[70]`; **full L-polynomial** identity with E(66)×E(858) at all
  **1227** good p < 10009, zero disagreements (the overnight `#J`-only check was
  299 primes < 2000); (3,3) not (2,2); Q(E1[3]) = Q(E2[3]), both degree 24.
- Deformation: **156,695** lines + **639,979** polynomial arcs of degree 2 and 3
  = **796,674** Q(t)-families tested across 9 anchor/target pairs,
  **0 survivors**; hardest kill needed p ≤ 17, t ≤ 7. The arc set includes
  279,993 families *tangent* to the split locus, which no line test can see.
- Positive control: density **1.0000** on contact-7 and contact-9 at 6 primes.
- Transverse density ≈ baseline at all 9 anchors; in-locus density 3.3–5.9×
  baseline.
- Z/63: contact-7 swept to height 2000 = **23,681,372,055,201 points**, 16-prime
  filter, 7 survivors, **all 7 verified false positives** (exact torsion `[7]`);
  contact-9 swept to **10,942,785,477 points**, 12-prime filter, 8 survivors,
  **all 8 verified false positives** (exact torsion `[9]`). 16 survivors in all,
  16 killed exactly.
- Census: compatible abs-simple Weil classes at 53–60 of 61 primes for every
  target; p = 3 forced bad reduction for 5 of the 8 targets.
- LMFDB alpha: **0 of 6,216,959** curves with a lane-8 torsion group is
  geometrically simple; controls Z/30 (6) and Z/40 (1) are.

---

## 7. Resume commands

```bash
# polynomial ARCS through an anchor (degree 2 or 3), both classes
ANCHOR=3 DEG=3 NARC=20000 gp -q -s 512M -f code/claude_ov_lane8_arcs.gp \
    > results/claude_ov_lane8_arcs_a3_d3.log 2>&1

# the deformation experiment (one process per anchor 1..9; ANCHOR=-1 = control only)
ANCHOR=3 PMAX=300 BOX=15 NRDIR=5000 gp -q -s 512M -f code/claude_ov_lane8_transverse.gp \
    > results/claude_ov_lane8_transverse_a3.log 2>&1

# f70: full L-factor identity, extend PMAX
PMAX=100000 gp -q -s 1G -f code/claude_ov_lane8_f70close.gp \
    > results/claude_ov_lane8_f70close.log 2>&1

# Z/63 sieve, contact-7 chart, EXTEND THE BOX (arguments: chart HN HD nthreads nprimes)
gcc -O3 -march=native -pthread -o /tmp/c63sieve code/claude_ov_lane8_c63sieve.c
/tmp/c63sieve 7 4000 4000 24 16 > results/claude_ov_lane8_c63_c7_h4000.log 2>&1
#   H=1000 -> 1.48e12 points in ~105 s;  H=2000 -> 2.37e13 in ~28 min (both done).
#   H=4000 -> 3.8e14 points, ~7.5 h on 24 threads; pair table 1.95e7 entries and
#   the per-thread candidate buffers are 2*8*NPAIR bytes each (~7.5 GB total at
#   24 threads), still fine on this box.  Beyond H ~ 6000 the O(NPAIR^2)
#   enumeration should be sharded by outer index across separate processes.

# Z/63 sieve, contact-9 chart (streaming, 1 parameter) -- use the OLD binary,
# the new one builds an O(HN*HD) pair table and will not fit
gcc -O3 -march=native -o /tmp/c79sieve code/claude_ov_lane8_c79sieve.c
/tmp/c79sieve 9 1000000 100000 > results/claude_ov_lane8_c63_c9_h1M.log 2>&1

# exact verification of any survivor
code/claude_magma_slot.sh -b CHART:=7 AS:="628/779" BS:="994/671" \
    code/claude_ov_lane8_c63verify.m > results/claude_ov_lane8_c63verify_c7.log 2>&1
code/claude_magma_slot.sh -b CHART:=9 AS:="201289/2665" \
    code/claude_ov_lane8_c63verify.m > results/claude_ov_lane8_c63verify.log 2>&1

# the corrected landscape census
gp -q -f code/claude_ov_lane8_census2.gp > results/claude_ov_lane8_census2.log 2>&1
```

Unfinished at the time of writing:

- `results/claude_ov_lane8_f70close.log` — the PMAX=20000 L-factor run was
  stopped by PID at p = 10009 (1227/1227 agree, 0 disagree) so as not to leave
  an orphan process; the log is committed at that point. Nothing depends on the
  remainder. To extend, re-run with a larger PMAX using the command above (cost
  per prime grows roughly linearly in p: p < 10^4 took ~1 h on one core).
Everything else finished. The contact-9 sweep completed
(`SEARCH_DONE c9 HN=300000 HD=30000 checked=10942785477 survivors=8`) and all
eight survivors were verified as false positives.

---

## 8. What I did **NOT** check

1. **No unconditional non-simplicity certificate for f70.** No explicit degree-3
   map C → E1, no Faltings–Serre, no endomorphism-algebra computation. The
   claim rests on L-factor matching at thousands of primes plus the literature
   provenance.
2. **Arcs were tested only to degree 3, and only by sampling.** Lines are
   covered exhaustively in an odd box `|ci| ≤ 15`; degree-2 and degree-3 arcs
   are covered by 20,000 random samples per class per anchor, not exhaustively.
   A transverse component that is a rational curve of degree ≥ 4 through the
   anchor, or one of degree ≤ 3 outside the sampled coefficient box `[−20,20]`,
   is not excluded. Nor is a component that meets the anchor only in a
   non-reduced or higher-multiplicity way.
3. **The line test excludes Q(t)-rational families only.** If the order-N section
   were defined over a degree-d extension `K(t)/Q(t)`, the mod-p density would be
   the proportion of Frobenius classes with a fixed point, not 1, and the test
   would not see it. The measured *equality with baseline* is evidence for full
   monodromy, but it is evidence, not proof.
4. **No genus computation for the marked covers of these anchors.** For [6,12]
   the analogous cover on a slice was shown irreducible of genus 51; here the
   analogous covers have degree ~(N⁴)/2 (e.g. 748,800 for N = 35) and the
   computation was not attempted. The mod-p density is the surrogate.
5. **Only the contact-7 and contact-9 charts were searched for Z/63.** Other
   normalisations of an order-7 or order-9 class (non-Weierstrass-difference
   supports, sextic models, the A(12)/M(2,12) family of charts) were not tried.
   A Z/63 curve whose order-7 class is not of contact type would be invisible to
   both searches.
6. **Weil-admissibility is a superset test.** The census's "compatible
   abs-simple class exists" is an upper bound: it does not exhibit a curve, only
   an isogeny class that is not excluded.
7. **No strict two-prime simplicity certificate was produced for anything** —
   because nothing geometrically simple was found. Every curve this lane
   *verified* is certified **split**, which is the opposite outcome.
8. **The LMFDB alpha table's `is_simple_geom`/`geom_end_alg` were taken as
   given**, not recomputed.
