---
name: target-playbook
description: Decompose a genus-2 rational-torsion target into a concrete route — coprime composition N=a*b, halving/doubling (the halving wall), contact for the prime-power part plus independent torsion, 2-rank via factor type, and the "+3" independent-3-torsion move — then map the mXXXX naming and find prior work in notes/. WHEN starting any new target ([2,24], [6,6], [4,16], Z/35, Z/48, Z/5xZ/5, [2,2,2n], order-60) or choosing a strategy before writing a search script.
---

# Target playbook: from a torsion target to a route

## When to use this

Load this the moment you are handed a torsion target `T` (e.g. `[2,24]`,
`[6,6]`, `Z/35`, `[4,16]`) and have to decide *how* to build it before you
write any Magma. This skill answers three questions:

1. Which **decomposition** turns `T` into reusable moves (coprime split,
   halving, contact + independent torsion, 2-rank engineering, `+3`)?
2. What is `T` **called** in this repo (`mXXXX` naming), and where is the
   **prior work**?
3. What is the **decision checklist** to go from `T` to a script to run?

Every strategy below is grounded in a real note. Read the cited note before
committing to a route — the notes record which variants are already dead.

## The five decomposition moves (each with a real example)

The methodology is a small algebra of moves (see the hub
`g2-torsion-lab`). Here is how to *choose* among them for a target.

### (a) COPRIME COMPOSITION `N = a*b`, `gcd(a,b)=1` — PREFERRED

Build a rational order-`a` class and an **independent** rational order-`b`
class. Because `gcd(a,b)=1`, they automatically compose to a class of order
`ab` — no halving, no high-degree unit. This is the preferred route whenever
`N` has a nontrivial coprime factorization.

- `24 = 8 x 3`: the `A(8)` chart carries a free order-8 class `D8`; any curve
  on it that *also* has a rational 3-torsion class has `Z/8 + Z/3 = Z/24` for
  free. This sidesteps both the `8->16` halving wall and the degree-24 Pell
  unit. Validated: two geometrically simple cyclic `Z/24` curves (D4 Frobenius
  certificates), e.g. `(r,p,t)=(5,-5/2,-9/2)` with `chi_17` irreducible.
  Cite `notes/agent_a2_24_composite.md`, code `code/agent_a2_24_composite8x3.m`.
- `20 = 4 x 5`: the quintic-contact-5 family carries an order-5 class; the
  explicit divisor `H=[x^2+2x/(t+1),(t+2)x+1]` with `2H=[x-1,0]` supplies the
  order-4 part, giving order 20 for free. Cite
  `notes/how_we_found_2220_examples.md`, `notes/contact5_order40_family.md`.
- `35 = 5 x 7`: contact-7 family plus an independent rational 5-torsion class.
  Cite `notes/agent_Z35_next_route.md` (the simultaneous contact7/contact5
  point equations).

**Key caveat (grounded):** coprime composition is cheap on paper but the two
sublocii can *barely intersect*. For `[2,24]` on `A(8)`, the order-8 + 2-rank-2
locus ([2,8]) is dense (~90% of the `W`-split family), but rational 3-torsion
on that same locus was absent through ~1M samples — "3-torsion is the wall".
No local obstruction was found, so `[2,24]` exists; it is a sparse global
point. Cite `notes/agent_a2_24_composite.md`. The lesson: when the coprime
factors live on nearly-disjoint sublocii, parametrize their **intersection**
directly (e.g. cubic-contact 3-torsion on the `W`-split family) rather than
scanning blindly.

### (b) HALVING / DOUBLING `n -> 2n`

An algebraic cover that doubles the order of a class: order `4->8`, `6->12`,
`[4,16]`. Powerful but **frequently obstructed** — the "halving wall". Halving
solutions often exist only on a split/nonsimple boundary, which is worthless
here (simplicity is the whole game).

- `Z/12 x Z/2 -> Z/12 x Z/4` on the one-parameter `M(12)` line `a=(1-r)/4`:
  obstructed. Good mod-7 residues `r=3,4` have `6D` divisible by 2 but the
  independent 2-torsion class not divisible; a combined residue sieve through
  height 300 leaves only the singular `r=-1`. Cite `notes/m12_z12x2_halving.md`.
- `[4,16]` second-halving on `M_1(8,4)`: exhaustively pushed. Through height
  150, on all six smooth `Q_7` strata, every verified rational first-half fails
  the final `P_R` second halving; a tangent-congruence sieve at
  `{17,23,29,47}` kills essentially all first-cover survivors. Local points
  exist (six smooth strata), so it is a height/global-geometry problem, not a
  `p=7` obstruction. Cite `notes/m18_416_component_and_smooth_strata_2026_07_02.md`.
- `[3,6] -> [3,12]` by halving the marked contact-6 class: the algebraic chart
  equals the standard `M(2,12)` chart, but through height 40 every simple
  boundary lift exact-tested to `[12]`, and the `[3,12]` phenomenon is confined
  to the split `Rinf+Z0` boundary (nonsimple). Cite `notes/contact6_m36.md`.

**Decision rule:** prefer coprime composition over halving whenever a coprime
factor exists. Reach for halving only when the target's prime-power part forces
it (e.g. reaching order `16` from `8`, or `[4,16]`), and expect to fight the
wall. Diagnose impossible-vs-hard with the boundary/obstruction skills before
spending a long height push.

### (c) CONTACT for the prime-power part + INDEPENDENT torsion

Use a **contact** construction (`f=h^2-c*(x-r)^n`) to force a marked class of
order equal to the prime-power part, then impose **independent** torsion (often
cubic-contact 3-torsion `f=h3^2+kappa*q3^3`) for the coprime part.

- `[6,6]` from the contact-6 chart `h6=1+a x+b x^2+x^3`, `f=h6^2-(x-1)^6`
  (marked order-6 class `D`, `P=(1,h6(1))`) plus an independent cubic-contact
  3-torsion class represented by `q=x^2+U x+v^2`, `h3^2-f=m^2 q^3`. Restricting
  to the `[1,2,2]` factor-type locus (2-rank 2) gives `[6,6]`. **Confirmed
  simple** example: `a=133/39, b=-7/13`, `L=29/16, U=-9/4, v=5/2`, torsion
  `[6,6]`, simple certificate `p=23`, `L_p=529 T^4-26 T^2+1` irreducible. Cite
  `notes/contact6_m36.md`, code `code/contact6_m36_core_slice_search.m`.
- `Z/48 = 16 x 3`-style: start from a verified order-16 (`A(16)`) point and
  add rational 3-torsion via the cubic-contact diagnostic
  `h^2-f=Lambda*q^3`, `q=x^2+U x+V`. The exact test is
  `TorsionSubgroup(...)` exponent divisible by 3. Current status: the order-16
  candidates all exact-test to `[16]` (`has3=false`); 3-torsion is the wall.
  Cite `notes/agent_Z48_cubic_contact_route.md`, code
  `code/agent_Z48_cubic_contact_route.m`.

Contact families are the standard carriers for order `5,6,7,8` prime-power
parts — see `contact-torsion-constructions` and `named-charts-reference` for
the explicit formulas.

### (d) 2-RANK via FACTOR TYPE

Engineer the factorization of `f` to produce the `[2,...]`/`[2,2,...]` prefix.
The rational 2-rank equals (number of rational-over-`Q` irreducible factors of
`f`) minus 1 (see `two-rank-and-factor-types`). Choose the factor type to build
in extra rational 2-torsion.

- `[2,2,20]`: on the contact-5 order-20 family, force the residual quartic
  `f/(x-1)` to have a rational linear factor (factor type `[1,1,2]`), which
  gives 2-rank 2 on top of the order-20 class. The `z=-1/7` specialization
  (`t=-8233/7225`) has type `1+1+2` and torsion `[2,2,20]`, geometrically
  simple (Lombardo true; Frobenius `X^4+2X^3+14X^2+142X+5041` at `p=71`
  irreducible under the 12th-power transform). Cite
  `notes/how_we_found_2220_examples.md`.
- `[6,6]`: restrict the contact-6 chart to factor type `[1,2,2]` (both
  quadratic factors of `f=x*((b+3)x^2+(a-3)x+2)*(2x^2+(b-3)x+(a+3))`
  irreducible) to get 2-rank 2 without falling into the many-rational-root
  extra-root locus. Cite `notes/contact6_m36.md` (core-cover section).

**Caveat (grounded):** forcing extra *rational roots* can push you onto a
split/nonsimple locus. In the contact-6 chart the extra-root type `[1,1,1,2]`
gave `[6,6]` curves that were all nonsimple (local `L`-polys factor into
quadratics); the `[1,2,2]` irreducible-quadratic locus is what yielded the
simple example. Cite `notes/contact6_m36.md`.

### (e) "+3" ROUTES — add independent rational 3-torsion

The `mXXXX_plus3` notes multiply an existing chart's target by 3 by imposing an
independent rational 3-torsion class (cubic-contact `h3^2-f=m^2 q^3`), filtered
by the necessary condition `3 | #J(F_p)` at every good `p != 3`.

- `M(2,2,2,4)+3` (aiming past `[2,2,2,4]`): first genuine good-reduction
  obstruction at `p=13`; residue-driven `p=13`-boundary enumerator + point-count
  primes to 73 leave only ~6 tuples per box, all killed by cubic-contact
  filters. Cite `notes/m2224_plus3.md`, code `code/m2224_plus3_search.m`.
- `M_1(8,2,2)+3` (i.e. `[2,2,8]+3`): open finite condition empty at
  `p=7,11,13`, so any example must be boundary at all three; strengthening the
  necessary condition through `p=73` kills every height-20 parameter. Cite
  `notes/m3222_plus3.md`, code `code/m3222_plus3_search.m`.

The `+3` move is the same cubic-contact machinery as move (c); the difference
is you are bolting 3-torsion onto an *already-built* even-torsion chart rather
than co-designing both. It is powerful but the good-prime `3 | #J(F_p)`
condition is very restrictive and repeatedly forces boundary reductions.

## How to PICK a route

1. **Read the target's note FIRST.** `ls notes/ | grep <target-digits>`, then
   read it. The notes record which variants are already dead and which branches
   are still live. Do not re-run a search the note already reports as cold.
2. **Check the live priority.** Read the newest `main_four_target_*` pass
   (currently `notes/main_four_target_fifth_pass_2026_07_02.md`) for the live
   worker assignments, the priority order, and "results landed". As of that
   pass the four active lines and their owners are: `Z/35`
   (compressed-state obstruction automaton, highest priority), `Z/5 x Z/5`
   (`b2=0` global/local sieve, second), `Z/48` (bounded cubic-contact
   production, background), `A(2,24)` (alternate scout, low). Treat any
   status as a snapshot; the newest pass wins.
3. **Prefer coprime composition** (move a) over halving (move b) whenever `N`
   has a coprime factorization. Halving is the fallback for a forced
   prime-power jump.
4. **Identify which chart carries the prime-power part.** `A(8)` for order 8,
   `A(12)` for order 12, `M_1(8,4)` for `[4,16]`, contact-`n` for a marked
   order-`n` class. See `named-charts-reference` for formulas.
5. **Add the coprime/2-rank part with the matching move:** cubic-contact
   3-torsion for a factor of 3; factor-type engineering for the `[2,...]`
   prefix; an independent order-5/7 contact class for those primes.

## The mXXXX naming (find work by name)

`m` + concatenated invariant-factor digits denotes the torsion group. Learn
these so `ls notes/ | grep` finds prior work:

```text
m2220   Z/2 x Z/2 x Z/20        ([2,2,20], found & simple)
m36     Z/6 x Z/6               ([6,6], found & simple)
m3222   Z/2 x Z/2 x Z/16        the [2,2,16] worker line
m2224   [2,2,2,4] base (m2224_plus3 = that base + 3)
m2226   [2,2,26]   m2228 [2,2,28]   m2248 [2,2,48]
m244    M(2,4,4) base (m244_to_248 = route toward [2,4,8])
m18_416 the M_1(8,4) chart for [4,16]
a2_24 / A(2,24)   Z/2 x Z/24    ([2,24])
Z35     Z/35        Z48  Z/48
z5x5    Z/5 x Z/5
contactN  the order-N contact construction
```

Cyclic targets use `Z<n>`; a `mXXXX_plus3` note is that base chart plus an
independent rational 3-torsion class.

## How to FIND prior work

```bash
# 1. Find the target's notes by its digits / name.
ls notes/ | grep 2220        # or 36, 416, Z35, Z48, z5x5, a2_24, plus3, ...

# 2. Read the note; it names the chart, the equations, and the code/results.
#    Convention: route foo has code/foo.m, notes/foo.md, results/foo*.log.

# 3. Read the referenced script and its captured log.
ls code/ | grep <route>
ls results/ data/ | grep <route>
```

The notes are the project memory: each records the chart, the exact equations,
every search that was run with its counts, and the conclusion (often negative,
deliberately, so you do not repeat it). Read the note before touching code.

## DECISION CHECKLIST: target `T` -> a concrete script

1. **Name it.** Write `T` in `mXXXX`/`Z<n>` form (see naming above).
2. **Find prior work.** `ls notes/ | grep <digits>`; read the note(s) and the
   newest `main_four_target_*` pass. If a route is already recorded as dead,
   do not repeat it — pick a different move or a live boundary branch.
3. **Factor `N`.** If `N=a*b` with `gcd(a,b)=1`, plan a coprime composition
   (move a): pick a chart that carries order-`a` for free and add an
   independent order-`b` class. This is the default.
4. **Locate the prime-power carrier.** Choose the named chart / contact family
   for the largest prime-power factor (`A(8)`, `A(12)`, `M_1(8,4)`, contact-`n`)
   — see `named-charts-reference`.
5. **Add the coprime part:** cubic-contact 3-torsion (`+3`, move e / c) for a
   factor of 3; a factor-type restriction (move d) for the `[2,...]` prefix;
   an independent order-5/7 contact class for those primes. If no coprime
   split exists, fall back to halving (move b) and expect the wall.
6. **Impose 2-rank if needed** (move d): pick the factor type of `f`
   (`[1,2,2]` for 2-rank 2, etc.) — cross-link `two-rank-and-factor-types`.
   Prefer irreducible-quadratic factor types over forced extra rational roots,
   to stay off the split/nonsimple locus.
7. **Build the funnel:** parametrize the chart -> cheap finite prefilter
   (`3 | #J(F_p)` for a 3-part, `n | #J(F_p)` for order-`n`, factor-type
   residues for 2-rank) -> exact `TorsionSubgroup(J)` -> simplicity
   certificate. See `running-torsion-searches`, `finite-prefilters`,
   `simplicity-certificates`.
8. **Copy the closest existing script** (`code/<sibling_route>.m`) and adapt
   its filters — do not invent Magma idioms; see `magma-lab-conventions`.
9. **If nothing survives,** decide impossible-vs-hard with `local-obstructions`
   and `component-boundary-analysis` (list the boundary components, check
   first-order lifts) before any wider blind height push.
10. **Before claiming a hit:** certify BOTH exact torsion and geometric
    simplicity — see `validate-and-record-a-hit` and `simplicity-certificates`.
    A `[6,6]`/`[3,12]` curve that splits is worthless.

## Worked route selections (grounded)

- `[2,24]`: coprime `24 = 8 x 3` on `A(8)`; 3-torsion is the wall on the
  2-rank-2 sublocus. Next: parametrize the 3-torsion sublocus of the `W`-split
  family directly. `notes/agent_a2_24_composite.md`.
- `[6,6]`: contact-6 order-6 class + independent cubic-contact 3-torsion,
  restricted to factor type `[1,2,2]`. **Done, simple.** `notes/contact6_m36.md`.
- `[2,2,20]`: quintic-contact-5 (order-20 via the explicit `H`,`2H`) + factor
  type `[1,1,2]` for 2-rank 2. **Done, simple.** `notes/how_we_found_2220_examples.md`.
- `Z/48`: order-16 `A(16)` point + cubic-contact 3-torsion (`48=16 x 3`);
  3-part is the decisive exact filter, currently the wall.
  `notes/agent_Z48_cubic_contact_route.md`.
- `Z/35`: contact-7 + independent 5-torsion (`35=5 x 7`); currently a
  compressed-state boundary-automaton problem, not raw height.
  `notes/agent_Z35_next_route.md`, `notes/main_four_target_fifth_pass_2026_07_02.md`.
- `Z/5 x Z/5`: literal contact-5/contact-5 is **empty off the same-contact
  boundary** (the fixed quartic `5X^4-10X^3+10X^2-5X+1` is irreducible over
  `Q`), so use a `b2=0` global/local sieve from smooth `F_7/F_11` charts
  instead. `notes/agent_z5x5_contact5_contact5.md`,
  `notes/agent_z5x5_b2zero_global_sieve.md`.
- `[4,16]`: `M_1(8,4)` second-halving; halving wall, quiet through height 150.
  `notes/m18_416_component_and_smooth_strata_2026_07_02.md`.
- Order 60: `M(12)+5` or contact-5 `[2,20]+3`; both boundary-obstructed, no
  point found. `notes/order60_attempts.md`.

## Pitfalls

- **Reaching for halving when a coprime split exists.** `24=8x3`,
  `20=4x5`, `35=5x7`, `48=16x3` are coprime compositions — do NOT try to halve
  from order `n` to `2n` first. Halving hits the wall; composition is free.
  (`notes/agent_a2_24_composite.md`.)
- **Assuming coprime composition is automatically easy.** The two sublocii can
  be nearly disjoint (the `[2,24]` "3-torsion is the wall" result). If a blind
  intersection scan is cold, parametrize the intersection locus directly.
- **Forcing extra rational roots for 2-rank and landing on a split curve.**
  Extra-root factor types (`[1,1,1,2]`) gave nonsimple `[6,6]`; the
  irreducible-quadratic type `[1,2,2]` gave the simple example. Prefer
  irreducible quadratic factors. (`notes/contact6_m36.md`.)
- **Skipping the note and re-running a dead search.** The notes deliberately
  record negative results (e.g. `Z/5xZ/5` contact-5/contact-5 is empty; the
  `M(12)` halving line is obstructed at `p=7`). Always `grep notes/` and read
  before scripting.
- **Trusting a stale status.** "found/open" in the hub or an old note is a
  snapshot. Re-read the target's note AND the newest `main_four_target_*` pass
  for current truth and current worker ownership.
- **Claiming a hit on exact torsion alone.** Both exact `TorsionSubgroup` and a
  geometric-simplicity certificate are required; a split `[6,6]`/`[3,12]` is
  worthless here. (`notes/contact6_m36.md` found `[2,6,6]`/`[6,6]` curves that
  were all nonsimple before the simple one.)
- **Ignoring the good-prime necessary condition for the coprime part.** For
  rational `k`-torsion, `k | #J(F_p)` at every good `p != k` is a cheap, strong
  prefilter — the `+3` and order-60 routes are largely *decided* by it. Use it
  before any exact torsion call. (`notes/m3222_plus3.md`, `notes/order60_attempts.md`.)

## See also

- `g2-torsion-lab` — the orientation hub (read first).
- `contact-torsion-constructions` — the `f=h^2-c(x-r)^n` marked-class and
  `f=h3^2+kappa q3^3` independent-3-torsion formulas used by moves (a),(c),(e).
- `named-charts-reference` — `A(8)`, `A(12)`, `M(2,12)`, `M_1(8,4)` formulas;
  which chart carries which prime-power order.
- `halving-and-doubling` — move (b) mechanics and the halving wall.
- `two-rank-and-factor-types` — move (d): 2-rank from the factor type of `f`.
- `simplicity-certificates` — required before claiming any hit.
- `local-obstructions`, `component-boundary-analysis` — impossible-vs-hard when
  a route yields nothing.
- `running-torsion-searches`, `finite-prefilters` — build and run the funnel.
- `validate-and-record-a-hit` — the hit checklist and how to document it.
