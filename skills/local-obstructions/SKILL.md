---
name: local-obstructions
description: Decide whether a torsion target is locally OBSTRUCTED (impossible) versus merely globally thin (hard) — p-adic solvability, Sha[2]/ELS-torsor signatures, and reading finite-prime evidence. WHEN a route yields nothing and you must tell "impossible" from "hard": no p-adic point compatible with the torsion structure, everywhere-locally-soluble but no global point (Sha[2]), many smooth mod-p points but no rational hit, or a finite filter that kills the good-open chart. Trigger words: local obstruction, ELS, Sha[2], p-adic solvable, good-reduction, thin set, boundary-confined.
---

# Local obstructions: impossible vs merely thin

## When to use this

A search or cover has produced **no rational point**. Before you burn more
compute or declare the target hopeless, you must decide which of two very
different situations you are in:

- **Locally obstructed / impossible** — there is *no* p-adic (or real) point
  at some prime compatible with the required torsion structure, or the
  descent torsors are everywhere-locally-soluble but carry no global point
  (a Sha[2] obstruction). In the first case there is a *finite certificate*
  of impossibility; in the second the certificate is a nontrivial Selmer/Sha
  class, not a mere point count.
- **Merely thin / hard** — the object is nonempty everywhere locally, has
  many smooth points modulo small primes, but the rational points form a
  *thin set* (a covering condition, a higher-dimensional arithmetic problem,
  or a positive-genus cover of low rank). Low-height search misses them; that
  is a statement about *difficulty*, not about *existence*.

Getting this call right decides whether to (a) stop and record an obstruction,
(b) hand off to `component-boundary-analysis` for structural understanding, or
(c) push height / build a better chart.

## What a local obstruction actually is

**Type 1 — no compatible p-adic point (good-reduction obstruction).**
Fix a prime `p` of good reduction. Enumerate the smooth points of the relevant
variety over `F_p`, and test whether any of them satisfies the *torsion-square
conditions* that the target imposes. If **none** does, then no rational point
can have good reduction at `p`; a rational point (if any) must reduce into the
bad/boundary locus at `p`.

The worked instance is `A(2,2,4,4)` on the K3 surface
`(ab+ac+ad+bc+bd+cd)^2 = 4abcd` (`notes/a2244_local_obstructions.md`). For a
fixed 2+2 partition `{i,j}|{k,l}`, the A(2,2,4,4) condition over `Q_p` is that
the four oriented cross-differences

```text
a_i^2 - a_k^2,  a_i^2 - a_l^2,  a_j^2 - a_k^2,  a_j^2 - a_l^2
```

all lie in the **same** `Q_p^*/Q_p^{*2}` squareclass (same valuation parity
*and* same Legendre symbol of the unit — the sign matters, since `-1` is a
nonsquare at both `11` and `23`). The finite-field enumeration found

```text
p=11: 240 good-reduction roots, 0 with any surviving A2244 partition
p=23: 5280 good-reduction roots, 0 with any surviving A2244 partition
```

So `p=11` and `p=23` are **genuine good-reduction obstructions**: any rational
point must have bad reduction at both. That is a real, certifiable local fact
(`notes/a2244_local_obstructions.md`, "Structural good-reduction obstructions"
and "p-adic residue-class refinement").

**Type 2 — Sha[2] / ELS-torsor obstruction.** Here every 2-descent torsor
(cover) is *everywhere locally soluble* (ELS) — it has `Q_v`-points for all
places `v` — yet has no global rational point. There is no single prime that
kills it; the obstruction is a nonzero class in `Sha[2]` (or a residual
2-cover after descent). The diagnostic is a *descent computation*, not a point
count. In the `m18_416` `[4,16]` work:

- For the `R = -25/4` fiber, ordinary elliptic two-descent on the minimal
  fiber `E`, after removing torsion and all free generators, returns
  `TwoDescent residual num_covers = 0` and `MordellWeilShaInformation = [3,3]`
  — i.e. the obstruction is **not** an unaccounted `Sha(E)[2]` class for that
  elliptic fiber; the real remaining object is a genus-5 square-condition
  cover (`notes/agent_m18_416_R25_4_SB_descent_followup.md`).
- The `R = -8` ELS fiber survives every rank-zero quotient certificate; all
  three `V4` elliptic quotients have rank 1, torsion `[2,2]`, and an MW-sieve
  on both natural rank-one models **does not close** through `p ~ 89`
  (`notes/agent_m18_416_R8_mwsieve_attempt.md`). This is the ELS-but-hard
  case: everywhere locally soluble, no cheap contradiction, needs a genuine
  Chabauty/cover-descent input.

The ELS language: a fiber is *ELS* when local solubility holds at every tested
place; being ELS is exactly what makes a target *not* Type-1-obstructed, and
pushes the question into descent / global territory.

## How the a2244 Python sieves compute the evidence

These are the reference tools for producing the finite evidence you need.
They are Python (an exception to Magma-only; they do pure residue arithmetic).

`code/a2244_padic_signature_sieve.py` — the candidate-level **signed** sieve.
For each integer tuple and each 2+2 partition it computes the oriented
cross-difference squareclass at `p^k` and keeps only tuples with a **common
finite-local partition at both obstruction primes**. Key routines:

- `squareclass_qp_signed(n, p)` returns `(valuation parity, Legendre(unit))`,
  or `"0"` — the *signed* class, because `squareclass_mod_pk` and the exact
  `Q_p` version both track the unit's Legendre symbol, not just parity.
- `partition_status_mod_pk(...)` returns `killed` / `deep` / `resolved_ok`:
  `killed` = two distinct squareclasses appear (impossible for that
  partition), `deep` = some cross-difference is `0 mod p^k` (undetermined,
  needs a deeper lift), `resolved_ok` = all cross-differences share one class.
- `real_partitions_for_sorted_positive(t)` — for a sorted positive curve
  tuple the only real-compatible partition is `12|34`; the real place is an
  extra necessary filter on top of the finite primes.

Run it as (defaults target `p=11,23`, depth 3):

```text
python3 code/a2244_padic_signature_sieve.py <tuples.txt> \
  --primes 11 23 --max-depth 3 \
  --output data/....txt --survivors data/..._survivors.txt
```

The depth ladder is the crux of the impossible-vs-hard call: on the 118
height-10000 forced-boundary tuples, `k=1` left 28 finite-common tuples (9
real-compatible), `k=3` left 14 (4 real), and depth 4 + the real sign
condition left **0** before the exact `M(2,2,4,8)` sieve
(`notes/a2244_local_obstructions.md`, "Candidate-level signed p-adic signature
sieve"). Deepening the p-adic precision is how you distinguish a shallow
`deep` (which may still lift) from a true `killed`.

`code/a2244_component_adic_analysis.py` — the **component-wise** analyzer. It
separates the ten raw K3 boundary components `Z1..Z4, E12..E34` and the three
partitions, enumerates all mod-`p` roots, exhaustively lifts mod-`p` ambiguous
classes to `p^2`, and searches for a **smooth `p^3` resolved witness** for each
(component, partition) pair. The verdict priority is

```text
smooth_p3_resolved > modp_resolved > p2_resolved > deep_only > killed
```

Its headline result is *negative for a component-wise obstruction*: for both
`11` and `23`, all 30 (component, partition) pairs are `smooth_p3_resolved` —
no single boundary component, even with a chosen partition, is ruled out
locally. Read this as: **the boundary labels alone are too coarse to prove
impossibility**; the obstruction, if any, needs the same partition
simultaneously at both primes plus the real sign plus the full cover.

`code/a2244_small_prime_congruences.py` — the **small-prime** congruence
diagnostic (`p <= 7`). For these primes there are fewer than four nonzero
square residues in `F_p`, so *every* primitive tuple has bad reduction; the
script measures whether the finer squareclass condition still yields a useful
congruence filter (`p=2` depth 5, `p=3` depth 3, `p=5` depth 2, `p=7` depth
1). Use it to learn that `5` and `7` are **filters/degeneracy behavior**, not
good-reduction obstructions (`notes/a2244_local_obstructions.md`,
"Structural good-reduction obstructions": the small primes "cannot support
four distinct nonzero squared branch residues").

## The KEY interpretive rule

This is the single most important judgment in this skill. Read the two
patterns and match your evidence to one of them.

**Pattern A — thin/hard, NOT impossible.** *Many smooth nonboundary points
modulo small primes, but no low-height rational point.* This is a genuine
thin / higher-dimensional arithmetic problem, not geometric emptiness. The
template is the `[3,6]` cover in `notes/contact6_m36.md` ("Component analysis
of the [3,6] cover"): the nonboundary cover is nonempty mod every tested prime

```text
p=5:  generic 4    p=13: generic 90   p=23: generic 224
p=7:  generic 6    p=17: generic 98   p=29: generic 350
p=11: generic 40   p=19: generic 310  p=31: generic 1534
```

yet the low-height rational searches find only boundary points, and the
`G2core`/`G3core` (deg 18 / deg 10) are irreducible. The note's own verdict:
"the obstruction is not geometric emptiness ... a genuine higher-dimensional /
thin arithmetic problem." **Do not call this impossible.** The same pattern
appears in `m18_416`: no local obstruction up to `p=101`, dense p-adic
survivor counts, and MW sieves that do not close
(`notes/agent_m18_416_R8_mwsieve_attempt.md`) — "useful negative information:
the remaining case is not a small-prime or plain finite-MW-sieve artifact."

**Pattern B — a real obstruction on the open chart.** *A finite filter KILLS
the good-open chart at some prime.* Then the target phenomenon is confined to a
boundary / nonsimple locus — an actual obstruction on the open (simple) chart.
The template is the `[3,12]` halving in `notes/contact6_m36.md`
("Component-wise p=5 boundary analysis"): the mod-5 filter gives

```text
p=5: allowed_312 0, bad 19, good rank counts [<1,6>].
```

i.e. **zero** good-open residues admit `[3,12]`; any rational `[3,12]` must
reduce to the bad/boundary locus mod 5. Boundary analysis then showed the
`[3,12]` phenomenon is **confined to the split `Rinf+Z0` boundary**, and every
simple-certified boundary lift exact-tested to `[12]`, not `[3,12]`. That is a
real obstruction *for simple examples on this chart* — even though `[3,12]`
itself is not globally impossible (the split, nonsimple examples exist).

The distinction in one line: **finite-nonemptiness with a killed *good-open*
chart is an obstruction on the open locus; finite-nonemptiness with a healthy
good-open chart but no rational point is thinness.**

## How to run the point-count evidence to make the call

1. **Establish good-open vs boundary emptiness.** Enumerate smooth
   *nonboundary* points of the cover modulo a range of small primes. If the
   good-open count is `0` at some prime (Pattern B, e.g. `allowed_312 0`), you
   have a good-open obstruction — go to `component-boundary-analysis` to find
   which boundary stratum the phenomenon is confined to. If the good-open
   counts are healthy and growing (Pattern A), the cover is not empty; the
   problem is thinness.
2. **Run the signed depth ladder** (`a2244_padic_signature_sieve.py` style):
   test squareclass status at `p, p^2, p^3, ...`. A partition that stays
   `killed` across the ladder is a true local obstruction at that prime; a
   partition that is only ever `deep` needs a deeper lift before you may
   conclude anything. Track the **real place** as an extra necessary
   condition.
3. **Check component-wise resolvability** (`a2244_component_adic_analysis.py`
   style): if every (component, partition) pair is `smooth_p3_resolved`, no
   single component is a local obstruction — impossibility, if it exists, is a
   *joint* condition across primes + real + the full cover.
4. **If everywhere-locally-soluble (ELS), switch to descent.** Point counts
   will not settle a Sha[2] case. Run a 2-descent / rank-zero-quotient
   certificate on the relevant elliptic/`V4` fibers (see the `m18_416` ELS
   notes) and, if the quotients have positive rank, an MW sieve or Chabauty.
   A rank-zero killer quotient (e.g. `R=-29/8`, `R=-25/4`) is a genuine
   emptiness certificate; a non-closing MW sieve (e.g. `R=-8`) is *not* —
   record it as ELS-but-open.

## Pitfalls

- **Calling a dense-mod-p cover "impossible."** Many smooth nonboundary points
  mod small primes with no low-height hit is **Pattern A (thin/hard)**, never
  a geometric-emptiness proof. Cite the `[3,6]` generic counts before you make
  any impossibility claim; if the good-open chart is nonempty, do not say
  impossible.
- **Reading a `deep` status as `killed`.** In the signed sieve, `deep` means a
  cross-difference vanished mod `p^k` — the class is *undetermined*, not
  killed. Every mod-`p` ambiguous root-partition in the a2244 study had at
  least one deep `p^2` lift; the `p^2` refinement did **not** prove
  impossibility (`notes/a2244_local_obstructions.md`, "p-adic residue-class
  refinement"). You must deepen the lift, not conclude.
- **Dropping the sign.** Squareclass is *signed*: at `p=11, 23` the value `-1`
  is a nonsquare, so parity alone is wrong. Use the signed `(parity,
  Legendre(unit))` class (`squareclass_qp_signed`), not just valuation parity.
- **Forgetting the real place.** A partition can be finite-locally common at
  all primes and still be killed by the real sign condition (the sorted
  positive tuple only allows `12|34`). The a2244 depth-4 survivor was killed
  precisely by the real place. Always intersect with the real-compatible
  partitions.
- **Mistaking a non-closing MW sieve for an obstruction.** An MW sieve that
  exceeds its CRT cap without closing (e.g. `R=-8`, classes truncated at
  `200001`) proves **nothing** about emptiness — it is ELS-but-open, needing a
  global input (Chabauty/cover-descent). Do not record it as impossible.
- **Confusing "boundary-confined" with "globally impossible."** In the
  `[3,12]` case the phenomenon exists on the split `Rinf+Z0` boundary — the
  target is real, just not on the *simple open* chart. "Confined to a
  nonsimple boundary" is the correct, narrower claim.

## See also

- `component-boundary-analysis` — the next step once you have a good-open
  obstruction: eliminate variables, factor components, saturate spurious
  branches, and localize the phenomenon to a boundary stratum at the
  obstructing prime.
- `pell-cf-order` — exact `D_infinity` order; used to certify torsion order
  before you ask whether a target is locally attainable.
- `finite-prefilters` — the cheap `J(F_p)` necessary conditions whose *failure*
  on the good-open chart is the Pattern-B signal.
- `target-playbook` — where an "impossible on this chart" verdict sends you to
  a different route.
- `g2-torsion-lab` — the hub; read it first.
