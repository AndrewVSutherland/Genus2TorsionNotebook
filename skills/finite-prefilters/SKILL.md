---
name: finite-prefilters
description: Cheap necessary conditions over F_p that gate the expensive exact stages of a torsion search — k|#J(F_p) divisibility filters, subgroup-compatibility filters, residue sieves, how to calibrate their strength, and the false-positive accounting that budgets the exact stage. WHEN designing a search funnel, when a prefilter passes "too many" or "too few" candidates, when deciding how many primes to test, or when a filter killing the good-open chart might actually be a local obstruction. Trigger words - prefilter, 3|#J(F_p), divisibility filter, necessary condition, good prime, residue sieve, false positives, PreP.
---

# Finite prefilters

## When to use this

Load this when building or debugging the cheap-filter stage of a torsion
search. Exact `TorsionSubgroup(J)` over `Q` is orders of magnitude more
expensive than a handful of point counts over small finite fields, so every
search funnel (see `running-torsion-searches`) interposes finite necessary
conditions between the parametrization and the exact stage. This skill gives
the standard filters, their verbatim implementations, and — most importantly —
how to reason about their strength and their false positives.

## The principle

For a prime `p` of **good reduction** with `p` odd and coprime to the torsion
order `k`, reduction `J(Q)_tors -> J(F_p)` is injective on prime-to-`p`
torsion. Hence a rational `k`-torsion point forces

```text
k | #J(F_p)   at EVERY good prime p (p coprime to 2k).
```

One failing good prime **kills** the candidate outright. Passing any number of
primes proves **nothing** — the condition is necessary, never sufficient.
These two asymmetric facts drive the whole design: filters are for *rejecting
cheaply*, and every survivor still owes the full exact bill.

## The standard divisibility filter (verbatim)

From `code/agent_a2_24_wsplit_3tors.m` (the 14-prime production version; the
earlier `code/agent_a2_24_composite8x3.m` used the 10-prime
`[7,11,13,17,19,23,29,31,37,41]`):

```magma
PreP := [7,11,13,17,19,23,29,31,37,41,43,47,53,59];

function ThreeTorsionPrefilter(fInt)
    for q in PreP do
        PF := PolynomialRing(GF(q));
        fp := PF![GF(q)!co : co in Coefficients(fInt)];
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        if #Jacobian(HyperellipticCurve(fp)) mod 3 ne 0 then return false; end if;
    end for;
    return true;
end function;
```

Two details are load-bearing:

- **Bad reduction gives NO information — `continue`, never `return false`.**
  A prime where `f mod p` drops degree or is not squarefree says nothing about
  rational torsion; failing the candidate there is a correctness bug, and
  "passing" it there is not evidence.
- The filter runs on the **integral model** `fInt` (denominators cleared), so
  the mod-`p` reductions are meaningful. See `magma-lab-conventions`.

## Strength calibration and false positives (the load-bearing lesson)

For a curve **without** rational 3-torsion, `3 | #J(F_p)` still holds at a
positive density of primes — the density is a Chebotarev average over the
mod-3 Galois image. An AND over `n` primes therefore has residual pass rate
roughly `(density)^n`: small, but **never zero**, and the survivors are not
random — they are systematically the curves with small mod-3 image, which look
maximally like they have rational 3-torsion without having it.

The production numbers to remember (`notes/agent_a2_24_composite.md`, the
"3-torsion is the wall" section, W-split H=26 run):

```text
~1,000,000 curves passed the 2-rank-2 stage
     1,709 passed the 14-prime 3-divisibility prefilter   (~1 in 580)
         0 had genuine rational 3-torsion (exact stage: histogram [2,8]:1708)
```

Read that table until it hurts: a **1/580 concentration** and **100% false
positives** are simultaneously true. The filter did its job (it made the
exact stage affordable); it just cannot do the exact stage's job. Two rules
follow:

1. **Budget exact-stage time for every survivor.** The funnel's cost model is
   `(#survivors) x (cost of TorsionSubgroup)`, not zero.
2. **A prefilter-passing curve is a candidate, never a hit.** Log it as such;
   see `validate-and-record-a-hit`.

Choosing the prime count: each extra prime multiplies the residual density by
the per-prime average (empirically a few tens of percent for 3-torsion), and
costs one `#Jacobian` point count (~`q` field operations at prime `q`). The
lab's production choice — 10 to 14 primes, smallest first — is a good default:
cheap primes reject first, and beyond ~14 primes the marginal rejection rarely
pays for the added cost.

## Don't filter what the construction guarantees

On the A(8) chart the visible marked class `D8` forces `8 | #J(F_p)` at every
good prime automatically — an "8 | #J" prefilter there is pure waste. The same
holds for any named chart's marked order (see `named-charts-reference`).
Filter only the conditions the construction does NOT already impose (for
`[2,24]` on A(8): the 3-part, since order-8 and 2-rank come from the chart and
the W-split).

## Subgroup-compatibility filters (stronger than divisibility)

Bare divisibility `k | #J(F_p)` is the weakest necessary condition. When the
target has structure (e.g. `[3,12]`), require instead that `J(F_p)` **contain
a subgroup compatible with the full target** — i.e. its invariant factors
admit the target's. This is the filter used by the `[3,12]` search
`code/contact6_m36_halveD_m312_search.m` (compatibility of `J(F_p)` with
`[3,12]` at every good `p != 2,3`), and its per-prime allowed-residue counts
are the `allowed_312` numbers in `notes/contact6_m36.md`.

**The diagnostic edge case:** if at some prime the count of good-open-chart
residues compatible with the target is **ZERO** (the `[3,12]` case:
`p=5: allowed_312 0`), that is no longer a filter — it is a **Pattern-B local
obstruction on the open chart**: every rational example must reduce to the
boundary at that prime. Stop filtering and switch to `local-obstructions` /
`component-boundary-analysis`. A filter that rejects everything is telling you
something structural.

## Residue sieves (the heavier variant)

For targets where divisibility filters are too weak, the lab uses mod-`p`
**allowed-residue enumeration**: enumerate the parameter residues mod `p` that
survive the target's conditions, then only lift/scan parameters in those
residue classes (CRT across several primes). This inverts the funnel — instead
of testing candidates, you enumerate the survivors.

Grounded examples (read the notes for the full pattern):

- `notes/m2224_plus3.md`, `code/m2224_plus3_search.m` — first genuine
  good-reduction obstruction at `p=13`; a residue-driven `p=13`-boundary
  enumerator plus point-count primes through `73` reduced each search box to
  ~6 tuples, all then killed by cubic-contact filters.
- `notes/m3222_plus3.md`, `code/m3222_plus3_search.m` — the open finite
  condition is EMPTY at `p=7,11,13`, so any example must be boundary at all
  three; strengthening the necessary condition through `p=73` killed every
  height-20 parameter.

The pattern to copy: find the obstructing prime(s), enumerate boundary/allowed
residues there, lift residue-by-residue, and strengthen with point-count
primes until the box is empty or a survivor emerges.

## Cost ordering in the funnel

Order the stages by unit cost (see `running-torsion-searches` for the full
funnel): free algebraic conditions (degree, discriminant, factor-type/2-rank)
first; then divisibility prefilters (smallest primes first — a rejection at
`q=7` costs ~7 operations, one at `q=59` costs ~59); then any marked-class
exact order check (a few Jacobian additions); `TorsionSubgroup` last, and the
simplicity certificate only on actual hits. Instrument every stage with a
counter and print them in the `PROGRESS` lines so the funnel's shape is
visible in the log.

## Pitfalls

- **Failing a candidate at a bad-reduction prime.** Bad reduction = no
  information = `continue`. Returning `false` there rejects good candidates
  (and has silently narrowed searches before being caught in review). Check
  the `Degree ... lt 5 or not IsSquarefree` guard is a skip, not a kill.
- **Treating survivors as hits.** 1709 survivors, 0 genuine — the prefilter
  concentrates, it never confirms. Every claim goes through exact
  `TorsionSubgroup` + the simplicity certificate
  (`validate-and-record-a-hit`, `simplicity-certificates`).
- **Filtering the chart's own guarantee.** `8 | #J(F_p)` on A(8) is a
  tautology; testing it wastes a point count per prime per candidate.
- **Using the rational model instead of the integral one.** Reductions of a
  non-integral `f` are meaningless or crash; clear denominators first
  (`magma-lab-conventions`).
- **Reading "filter rejects everything" as "search harder."** Zero compatible
  good-open residues at a prime is a local obstruction (Pattern B), not a
  filtering problem. Diagnose with `local-obstructions` before widening the
  scan.
- **Testing `p | k` primes.** For `k`-torsion, primes dividing `k` (and `2`)
  are excluded from the injectivity statement; the lab's `PreP` lists start at
  7 for 3-torsion precisely to stay clear. Keep that property when editing the
  prime list.

## See also

- `running-torsion-searches` — the funnel this skill's filters slot into, and
  the stage-counter instrumentation.
- `two-rank-and-factor-types` — the free algebraic filter that precedes these.
- `local-obstructions` — when a filter's zero-count is an obstruction, not a
  filter (Pattern A vs Pattern B).
- `named-charts-reference` — what each chart already guarantees (so you don't
  filter it).
- `simplicity-certificates`, `validate-and-record-a-hit` — the exact stages
  the survivors still owe.
- `magma-lab-conventions` — integral models, and the cost discipline for
  expensive intrinsics.
- `g2-torsion-lab` — hub.
