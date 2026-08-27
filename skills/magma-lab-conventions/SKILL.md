---
name: magma-lab-conventions
description: Magma idioms and the recurring pitfalls of this lab - the Factorization unit trap, RationalReconstruction arity, function-field coercion failures (FldFunRatMElt), -b batch buffering, parameter passing, memory caps, multi-value returns, Mumford-coordinate scaling on integral models, and the cost discipline for expensive intrinsics. READ BEFORE writing or debugging ANY Magma script in this repo. Trigger words - Magma error, bad argument types, Factorization, RationalReconstruction, FldFunRatMElt, buffering, -b, assigned, StringToInteger, SetMemoryLimit, elt<J|...>, runtime error.
---

# Magma lab conventions

## When to use this

Load this before writing or debugging any `.m` script here. Every item below
is a convention or a landmine that has actually cost this project time —
several caused wrong mathematical conclusions that had to be retracted. The
format is: the trap, the fix, and the file that grounds it.

## 1. THE FACTORIZATION UNIT TRAP (read this one twice)

Magma's `Factorization(poly)` returns **monic** irreducible factors with
multiplicities and **silently drops the leading unit/constant**:

```text
Factorization(2*x^2 - 2)  ->  [ <x - 1, 1>, <x + 1, 1> ]     // the 2 is GONE
```

The product of the returned factors does **not** equal the input. This single
behavior caused THREE real bugs in this project: a lost factor `2` in a `c4`
quantity, a lost factor `-2` in a `C2` quantity, and — worst — a spurious
rank-0 quotient conclusion in the `[4,16]` `R=-8` work that had to be
publicly retracted (see the "Validation and a cautionary note" and ERRATUM
sections of `notes/agent_m18_416_R8_dA_quotients.md`; the `c4` factor-2 and
`C2` factor-(-2) incidents are recorded in
`notes/agent_m18_416_p7_blowup_notes.md`).

**The rule:** whenever the actual polynomial (not just its factor degrees)
matters, reconstruct the unit explicitly:

```magma
fac := Factorization(g);
unit := g div &*[ t[1]^t[2] : t in fac ];   // the dropped constant
assert unit * &*[ t[1]^t[2] : t in fac ] eq g;
```

or track `LeadingCoefficient(g)` separately. Degree-only uses (the `TwoRank`
function, factor-type tests) are safe — see `two-rank-and-factor-types`.

## 2. Return arities and multi-value returns

- `RationalReconstruction` returns **two** values, `ok, value` — not three.
  Ground: `okc, crat := RationalReconstruction(Integers(pk)!cc);` in
  `code/agent_a2_24_construct_lift.m`.
- A function returning `k` values must return `k` things on **every** path:
  `return false,_,_,_,_,_;` with the right number of `_` placeholders. A
  mismatch is a runtime error raised only when that branch executes — a search
  can run for hours before hitting it. (This bit the lab: an early-return for
  even-degree cases needed six placeholders, not five.)

## 3. Function fields Q(r,t): coercion failures and the clearing idiom

Over a multivariate rational function field `K<r,t> := RationalFunctionField(Q,2)`,
the intrinsics `LCM`, `LeadingCoefficient`, and friends **fail on
`FldFunRatMElt` arguments** with `Bad argument types`. Denominators must be
coerced into the underlying polynomial ring first. The working idiom
(verbatim pattern from `ClearToR3` in `code/agent_a2_24_d0_cover_ptsearch.m`
and `code/agent_a2_24_d0_saturate.m`):

```magma
NR := PolynomialRing(Q, 2);            // underlying poly ring of K
cs := Coefficients(poly);              // elements of K
D := NR!1;
for c in cs do if c ne 0 then D := LCM(D, NR!Denominator(c)); end if; end for;
// ... then per coefficient:
cK := cs[i]*(K!D);                     // now denominator-free in K
Nn := NR!Numerator(cK);                // land in the poly ring
```

Note both coercions: `NR!Denominator(c)` before `LCM`, and `K!D` before the
multiply. Skipping either reproduces the `Bad argument types` error.

## 4. Batch mode: buffering, columns, quit

- `magma -b` **buffers stdout**. A `-b` run killed by timeout or signal shows
  **nothing** — not "the output so far", nothing. Never diagnose a killed
  `-b` run from its empty log; rerun it to completion or to an explicit
  checkpoint. Corollary: long runs are launched with
  `nohup magma -b ... > logfile 2>&1 &` and monitored via the log file
  (see `running-torsion-searches`).
- `SetColumns(0);` at the top of every script — otherwise Magma wraps long
  lines and breaks greppable one-line log markers.
- End every script with `quit;` — without it `-b` sits at the prompt and the
  job never exits.

Ground: the headers/footers of `code/agent_a2_24_wsplit_3tors.m`,
`code/agent_a2_24_composite8x3.m` (and essentially every search script).

## 5. Parameter passing and memory caps

`magma -b Name:=value script.m` passes **strings**. Every parametrized script
uses the guard idiom (verbatim from `code/agent_a2_24_ztors_sample.m`):

```magma
if not assigned H then H := 18; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
```

House limits (from `AUTHORING.md`, learned the hard way): **at most 3
concurrent Magma jobs** (the machine has OOM'd at 6), default `MemGB:=3` per
search job.

## 6. Jacobian points: construction, try/catch, and the scaling rule

- `J![u, v]` and `elt<J | u, v>` **throw** on invalid Mumford data. In search
  loops, always wrap:

  ```magma
  ok8 := true;
  try
      D8 := J![g8, v8];
      if 8*D8 ne O or 4*D8 eq O then ok8 := false; end if;
  catch ee ok8 := false; end try;
  ```

  (verbatim from `code/agent_a2_24_wsplit_3tors.m`).

- **The integral-model scaling rule.** Searches clear denominators:
  `fInt = L^2 * f` (the `IntModel` function returns both `fInt` and `L`).
  On the scaled model the substitution is `y -> L*y`, so the `y`-parts
  (Mumford `v`-polynomials) of divisors **scale by `L`**:

  ```magma
  fInt, Lden := IntModel(f);
  v8 := (-Lden*ellBase) mod g8;        // scaled model: v8 = Lden*(unscaled v8)
  ```

  (`code/agent_a2_24_ztors_sample.m`). Over a **finite field with the
  unscaled chart model**, do NOT multiply by the clearing scalar — the
  finite-field probe `code/agent_a2_24_locus_geometry.m` documents this
  explicitly (`v8 := (-ellBase) mod g8` with a comment). Mixing the two
  conventions makes the marked class silently fail its order check, which
  looks exactly like "the locus is empty".

## 7. Small solver idioms (copy these, don't invent)

- **Rational roots of a rational-function equation:**
  `Roots(Numerator(E))` after building `E` over a function field — the
  `SplitWPvals` pattern in `code/agent_a2_24_wsplit_3tors.m` (guard
  `En eq 0 or Degree(En) lt 1` first).
- **Tiny 0-dimensional solves:** build the coefficient equations and call
  `Variety(ideal<R | eqs>)` — the `ExtractContact` pattern in
  `code/agent_a2_24_contact_extract.m`.
- **Force a quantity nonzero in an ideal (Rabinowitsch):** add a fresh
  variable `w` and the generator `w*kappa - 1`; if the ideal becomes `(1)`,
  the system has no solution with `kappa != 0` —
  `code/agent_a2_24_d0_saturate.m`, and see `component-boundary-analysis`
  for what that means mathematically.
- **Dedupe + sort:** `Sort(Setseq(Seqset(vals)))` — the `HeightRationals`
  builders in every search script.
- **Formatted printing:** `printf` with `%o` for everything
  (`printf "HIT r=%o t=%o\n", rv, tv;`).

## 8. Cost discipline for expensive intrinsics

- `TorsionSubgroup(J)` over `Q` is the expensive endgame call. Precede it
  with: finite prefilters (`finite-prefilters`), and cheap exact order checks
  on a known divisor (`k*D` loops / `8*D8 ne O and 4*D8 eq O`-style tests) —
  a few Jacobian additions each.
- `EllipticCurve(C, pt)`, `Jacobian(C)`, `GenusOneModel(C)` on high-degree
  **singular plane models** (e.g. a degree-12 model of a genus-1 cover) stall
  for hours or fail outright ("Failed to obtain a genus one model"). Do not
  wait on them: determine the rank/points **empirically** by direct rational
  point search on the plane model instead. Ground:
  `notes/agent_a2_24_d0_derivation.md` ("Magma's EllipticCurve/Jacobian on
  the deg-12 model is very slow, so the rank was settled by direct
  rational-point search"), implementation
  `code/agent_a2_24_d0_cover_ptsearch.m`; method details in
  `component-boundary-analysis`.
- Rational (fraction) arithmetic **slows several-fold at high height**; a
  scan's early rate overestimates its late rate (see
  `running-torsion-searches` for budgeting).

## 9. Shared helpers: torsion_cover_lab_utils.m

`code/torsion_cover_lab_utils.m` holds the shared `TC_*` functions — reuse
them instead of re-implementing:

```text
TC_SumInts, TC_MakeMonic, TC_HeightRationals, TC_FactorDegreeMults,
TC_ContainsPoint, TC_GoodReductionPolynomial, TC_PointCountGate,
TC_NormalizeTorsionInvariants
```

(Older scripts inline their own copies of `HeightRationals`/`TwoRank` etc.;
new scripts should prefer the utils where they fit, but NEVER change the
utils' semantics in place — dozens of scripts load them.)

## Pitfalls (the compressed checklist)

- `Factorization` dropped your unit constant (**§1** — three real bugs, one
  ERRATUM). Reconstruct the unit whenever the polynomial itself matters.
- `RationalReconstruction` has arity 2, not 3 (**§2**).
- `return` paths with the wrong number of `_` placeholders (**§2**).
- `LCM`/`LeadingCoefficient` on `FldFunRatMElt` — coerce via the underlying
  polynomial ring (**§3**).
- Diagnosing a killed `-b` run from its empty log (**§4** — buffering).
- Missing `SetColumns(0)` (wrapped log lines break monitors) or missing
  `quit;` (job never exits) (**§4**).
- Forgetting `StringToInteger` on passed parameters — comparisons against a
  string silently misbehave (**§5**).
- More than 3 concurrent Magma jobs, or no `SetMemoryLimit` (**§5** — OOM).
- Unwrapped `J![u,v]` in a loop — one invalid candidate kills the whole run
  (**§6**).
- Wrong Mumford scaling: `Lden` on the integral model, NO `Lden` on the
  unscaled finite-field model (**§6** — failure masquerades as an empty
  locus).
- Waiting on `EllipticCurve`/`GenusOneModel` for a big singular plane model
  (**§8** — use point search).
- Editing `torsion_cover_lab_utils.m` semantics in place (**§9**).

## See also

- `running-torsion-searches` — launching, sharding, monitoring, and killing
  the runs these conventions serve.
- `two-rank-and-factor-types` — the degree-only `Factorization` use that IS
  safe.
- `finite-prefilters` — the cheap gates before `TorsionSubgroup`.
- `component-boundary-analysis` — the point-search fallback for genus-1
  covers, and Rabinowitsch saturation.
- `validate-and-record-a-hit` — the ERRATUM convention for retracting a wrong
  claim.
- `g2-torsion-lab` — hub.
