---
name: g2-torsion-lab
description: Orientation hub for the geometrically-simple genus-2 Jacobian torsion project — the goal, repository map, target-naming convention, the shared methodology, the state of play, and an index of every other lab skill. Read this FIRST whenever you start work on this repo (searching for large rational torsion, debugging a Magma search, certifying simplicity, or picking up a torsion target like [2,24], [6,6], [4,16], Z/35, Z/48).
---

# The genus-2 simple-Jacobian torsion lab

## What this project is

We hunt for genus-2 curves `C: y^2 = f(x)` over `Q` (with `deg f in {5,6}`)
whose Jacobian `J = Jac(C)` is **geometrically simple** (no isogeny factors
over `Qbar`, i.e. `End(J_Qbar)` is an order in a field) **and** has a
prescribed large **rational torsion subgroup** `J(Q)_tors`. The mathematical
point: torsion of *simple* abelian surfaces over `Q` is far less understood
than for elliptic curves, and explicit simple examples with large torsion are
the raw material for the paper (Balakrishnan, Najman, Shnidman, Sutherland;
an early draft lived in `paper/main.tex`).

Two things must hold for a result, and **both must be certified**:
1. `TorsionSubgroup(J)` equals the target invariant factors, exactly; and
2. `J` is geometrically simple (a *local Frobenius certificate* or the
   Sage/Lombardo endomorphism test — see `simplicity-certificates`).

A `[6,6]` curve that turns out to cover an elliptic curve is worthless here;
simplicity is the whole game.

## Repository map

```text
code/     ~300 Magma (.m) + some Python (.py) scripts. One script per
          experiment; naming mirrors notes/. Shared helpers in
          code/torsion_cover_lab_utils.m (the TC_* functions).
notes/    ~80 markdown design/status notes. THIS is the project memory —
          each route has a note documenting the chart, equations, searches,
          and conclusions. Read the note before touching the code.
results/  Captured Magma run logs (*.log, *.txt).
data/     Larger run outputs and generated symbolic-equation dumps.
paper/    main.tex and supporting TeX; the write-up of confirmed results.
```

Convention: a route `foo` usually has `code/<foo>.m`, `notes/<foo>.md`, and
`results/<foo>*.log` together. To pick up any target, **find its note first**
(`ls notes/ | grep <target>`), read it, then read the referenced code.

## Target naming

`m` + concatenated invariant-factor digits denotes the torsion group:
`m2220 = Z/2 x Z/2 x Z/20`, `m36 = Z/6 x Z/6`, `m2226/m2228/m2248` the
`[2,2,26]/[2,2,28]/[2,2,48]` lines, `m3222` the `[2,2,16]` worker line,
`m18_416` the `M_1(8,4)` chart for `[4,16]`. Cyclic targets: `Z35`, `Z48`,
`a2_24` (= `A(2,24)` = `Z/2 x Z/24`). Constructions: `contactN` = the
order-`N` contact family.

## The methodology (how examples are built)

Every route is a composition of a few reusable moves. Learn these as verbs:

1. **Contact** — force a marked rational torsion class of order `n` by writing
   `f = h^2 - c*(x-r)^n` (marked class from `h - y`), or an *independent*
   rational 3-torsion class via `f = h3^2 + kappa*q3^3`. → `contact-torsion-constructions`.
2. **Named charts** — normalized parametric families that carry a marked
   high-order class for free: `A(8)` (order 8), `A(12)` (order 12), `M(2,12)`,
   `M_1(8,4)`. → `named-charts-reference`.
3. **Halving / doubling** — algebraic covers that double a class's order
   (`n -> 2n`); beware the "halving wall". → `halving-and-doubling`.
4. **2-rank control** — pick the factor type of `f` ([1,2,2], [1,1,2,2], ...)
   to build in extra rational 2-torsion. → `two-rank-and-factor-types`.
5. **Compose coprime orders** — a rational order-`a` class plus a rational
   order-`b` class with `gcd(a,b)=1` gives order `ab` for free (e.g.
   `24 = 8 x 3`). → `target-playbook`.

Then the **search funnel** turns a family into examples:
parametrize → cheap **finite prefilter** (`J(F_p)` divisibility) →
exact `TorsionSubgroup(J)` → **simplicity certificate**. →
`running-torsion-searches`, `finite-prefilters`, `simplicity-certificates`.

When a route yields nothing, decide *impossible vs merely hard* with
`local-obstructions` and `component-boundary-analysis` before moving on.

## State of play (confirmed simple examples exist; see notes for the frontier)

Geometrically simple examples that have been found and documented include
`[2,2,20]` (`notes/how_we_found_2220_examples.md`), `[6,6]`
(`notes/contact6_m36.md`), and cyclic `Z/24` (`notes/agent_a2_24_composite.md`).
Active/open frontier lines (see the `main_four_target_*` notes for the live
priority list): `Z/35`, `Z/5 x Z/5`, `Z/48`, `[4,16]` (`m18_416`), `[2,24]`
(`a2_24`), `[2,2,16]` (`m3222`), and the `[2,2,2n]` families. Treat any
"found/open" status as a snapshot — **always re-read the target's note and the
newest `main_four_target_*` pass for current truth** before acting.

## Skill index (load the one that matches your task)

Construction:
- `contact-torsion-constructions` — build a marked torsion class by contact
  (`f=h^2-c(x-r)^n`) and independent 3-torsion (`f=h3^2+kappa q3^3`).
- `halving-and-doubling` — double a class's order; the halving wall.
- `named-charts-reference` — A(8), A(12), M(2,12), M_1(8,4) formulas.

Verification:
- `simplicity-certificates` — certify geometric simplicity (Frobenius D4
  certificate; Lombardo test). **Required before claiming any result.**
- `two-rank-and-factor-types` — 2-rank from factor type; targeting extra 2-torsion.
- `validate-and-record-a-hit` — the checklist + how to document a hit.

Running & pruning:
- `magma-lab-conventions` — Magma idioms, the TC_* utils, and the recurring
  pitfalls. **Read before writing/debugging any script.**
- `running-torsion-searches` — the funnel, memory caps, ≤3 jobs, background
  jobs, monitoring, height enumeration.
- `finite-prefilters` — cheap `J(F_p)` necessary conditions.

When stuck:
- `local-obstructions` — is the target locally impossible or just thin?
- `component-boundary-analysis` — eliminate/factor a cover, saturate spurious
  branches, compute genus/rank, boundary analysis at the obstructing prime.
- `pell-cf-order` — exact order of `D_infinity` via continued fractions.

Strategy:
- `target-playbook` — decompose a target into a route; find prior work.

## Operating rules (do not violate)

- Magma only. At most **3 concurrent Magma jobs**; cap memory
  (`SetMemoryLimit(3*10^9)` is typical). The machine has OOM'd at 6 jobs.
- Long runs: use `-b` batch mode, redirect to a `results/` or `data/` log,
  run in the background, and poll the log — **`-b` buffers stdout, so a run
  killed by timeout shows nothing**. Never rely on partial stdout of a piped
  `-b` run.
- Commit locally with a clear message; **the user pushes to GitHub** (agent
  push fails on auth). Co-author line: `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.
- Never claim a curve is a result until BOTH exact torsion and a simplicity
  certificate are in hand.

## See also
Start here, then load the specific skill for your task from the index above.
