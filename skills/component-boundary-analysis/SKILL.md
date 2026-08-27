---
name: component-boundary-analysis
description: When a cover yields no rational points, turn a failed search into structural understanding — eliminate a linearly-appearing variable to get the defining equations, factor into components, saturate spurious branches (Rabinowitsch), compute genus/rank, and do boundary analysis at the obstructing prime. WHEN a route dead-ends and you must know whether the core is thin, whether "has a rational torsion point" is a covering condition rather than a subvariety, and which boundary stratum the phenomenon is confined to. Trigger words: eliminate variable, resultant, component, Rabinowitsch, saturate spurious branch, covering condition, boundary stratum, genus, rank-0.
---

# Component and boundary analysis of a failed cover

## When to use this

A parametrized cover produced no rational hit (or only boundary hits). Rather
than a wider blind height search, extract the *structure*: what are the
defining equations, do they factor, is the "core" thin or does it carry a
rational component, is the torsion condition even a subvariety of the base, and
— if a prime obstructs the good-open chart — which boundary stratum does the
phenomenon live on. This is the tool that converts "search found nothing" into
"here is why, and here is the only place left to look."

Pair it with `local-obstructions`: that skill tells you *whether* you are in a
thin (Pattern A) or good-open-obstructed (Pattern B) situation; this skill
gives you the equations and the boundary decomposition to act on it.

## 1. Cover elimination + component analysis

Template: `notes/contact6_m36.md`, "Component analysis of the `[3,6]` cover"
(`code/contact6_m36_component_analysis.m`).

The recipe when one variable appears **linearly** in one of the defining
equations:

1. **Solve for the linear variable.** On the generic branch `v^3 != 1` the
   first cubic-contact equation `F1=0` has coefficient
   `D = 8*L^2*(v^3 - 1)` in `a`, so it determines `a = N/D`.
2. **Substitute and clear denominators** to get numerator equations. Multiply
   through by the right power of `D` so the result is polynomial:

   ```text
   G2 = numerator( D^2 * F2(a = N/D) ),
   G3 = numerator( D   * F3(a = N/D) ).
   ```

3. **Remove known boundary factors.** Compute the gcd of the `Gi` and divide
   out any factor that is a *known* geometric-degeneration locus. Here the
   only common factor over `Q` was `L^2` (the built-in `L^2 = 1/m^2` scale).
   Then check the residual against the explicit nonboundary product

   ```text
   L * v * (v^3-1) * (U^2-4*v^2) * (b+3) * (numerator of a+b+2),
   ```

   which lists every geometric-degeneration locus (`m=0`, `v=0`, the `v^3=1`
   branch cut, the `U^2=4v^2` degeneracy, the degree-drop `b=-3`, the
   `h6(1)=a+b+2=0` locus). No further boundary factor divided the cores.
4. **Test the remaining core.** If it is **irreducible and high-degree**, the
   generic branch is thin (no rational component to exploit). If it
   **factors**, you have found a rational component to feed forward.

**Worked outcome (the diagnostic to compare against).** After stripping `L^2`,
the two cores were irreducible over `Q`:

```text
G2core: total degree 18, terms 116
G3core: total degree 10, terms 39
```

with the special branch `v=1` giving an irreducible degree-6 `H` and an
irreducible degree-18 `Jcore` (gcd 1). Because the cores are irreducible and
high-degree — and because the finite-field counts show many smooth nonboundary
points — the conclusion is **thin, not empty**: "the generic branch is not
collapsing to an obvious boundary or rational factor." That is the signature
you are matching your own cover against.

## 2. The covering-condition subtlety (do NOT skip)

Template: `notes/agent_a2_24_d0_derivation.md`,
`code/agent_a2_24_d0_saturate.m`.

The trap: "has a rational torsion point of order `n`" is a **covering
condition**, NOT a hypersurface on the base. Over `Qbar` *every* curve in the
family carries the geometric `n`-torsion, so if you eliminate the
torsion-witness variables (the coefficients of the auxiliary cubic `h3`, the
contact point, etc.) you get **no condition on the base parameters** — the
elimination ideal over `Q(params)` is trivial. Eliminating naively "finds" the
whole base and tells you nothing.

**The Rabinowitsch test to detect this.** Add a variable `w` and the relation
`w*kappa - 1` to force `kappa != 0` (here `kappa` is the genuine-3-torsion
coefficient; `kappa=0` is the perfect-square / spurious locus giving order 8,
not 24). Then form the ideal over `Q(params)` and check whether it collapses to
`(1)`. From `code/agent_a2_24_d0_saturate.m`:

```magma
kap := f6 - m^2; bN := 2*m*n - f5;
Eq4 := 3*kap*(f4 - n^2 - 2*m*j) - bN^2;
Eq3 := 27*kap^2*(f3 - 2*m*l - 2*n*j) + bN^3;
I := ideal<S | Eq4, Eq3, w*kap - 1>;   // Rabinowitsch: forces kap<>0
printf "Ideal is (1) [no genuine soln generically]? %o\n", I eq ideal<S|1>;
```

The ideal `<Eq4, Eq3, w*kappa-1>` over `Q(r,t)` collapsed to `(1)`. The
interpretation: **genuine rational 3-torsion is NOT cut out by any hypersurface
in the base parameters — it is a covering (thin-set) condition.** The genuine
`Z/24` curves are the rational points of the *cover*

```text
X = { (r,t,m) : Eq4 = Eq3 = 0, kappa = f6 - m^2 != 0 },
```

a space curve, not a locus in `(r,t)`.

**Then materialize the cover explicitly.** Take a resultant to get a plane
model of `X`, and study *that* curve's genus and rational points. Eliminating
`r` (or `t`) from `X` gave a degree-12 plane curve of **genus 1** (both
projections agreed), with the known example `B` at `(r,t,m) = (1/3,-1,2/9)`.
Now the question "how many simple `Z/24` on this slice" becomes "what is the
rank of this genus-1 cover" — a concrete, answerable question.

## 3. Boundary analysis at an obstructing prime

Template: `notes/contact6_m36.md`, "Component-wise `p=5` boundary analysis"
(`code/contact6_m36_halveD_p5_boundary_analysis.m`).

Use this when `local-obstructions` has shown the good-open chart is **killed**
at some prime `p` (Pattern B, e.g. `allowed_312 0` at `p=5`): every rational
point must reduce to the boundary, so enumerate and classify the boundary.

1. **Factor the discriminant over `F_p`** to enumerate the boundary
   components. Clearing denominators in the affine `(z,r)` chart over `F_5`,
   the sextic discriminant factored as

   ```text
   r^6 * (r+1)^21 * z^2 * (z+1)^11 * (z-1)^11 * H(z,r),
   ```

   giving affine boundary components `R0: r=0`, `Rminus1: r=-1`, `Z0: z=0`,
   `Zplus/Zminus: z=+-1`, and the residual high component `H`.
2. **Classify the full `P^1 x P^1(F_p)` residue classes.** Rational points may
   also have denominator divisible by `p`, so work in
   `P^1_z x P^1_r(F_5)`: `36` classes total, `6` good-open and `30` boundary.
   The known split example reduces to `Rinf + Z0`.
3. **Check whether the target phenomenon is confined to a split/nonsimple
   boundary.** Run the height search restricted to the surviving boundary
   classes, then exact-test. The height-40 run left `48` boundary residue
   survivors, `45` simple survivors, **`0` simple `[3,12]` hits** — every
   simple survivor exact-tested to `[12]`. The only `[3,12]` cases all sat on
   `Rinf + Z0` and Magma certified them split (automorphism group order 4, two
   degree-2 subcovers). **Conclusion: the `[3,12]` phenomenon is confined to
   the split `Rinf+Z0` boundary** — a real obstruction for *simple* examples on
   this chart, even though `[3,12]` itself is not globally impossible.

That "confined to a split/nonsimple boundary" verdict is the deliverable of
boundary analysis: it names the exact locus where the phenomenon lives and
certifies it is off the simple open chart.

## 4. Distinguishing finite-nonemptiness from thinness

The finite-field counts alone do NOT decide impossibility (see
`local-obstructions`). In the `[3,6]` core analysis the nonboundary cover was
nonempty mod every tested prime (generic counts 4, 6, 40, 90, ... up to 1534 at
`p=31`), yet the cores were irreducible high-degree. The right reading:
**finite-nonemptiness + irreducible high-degree core + no low-height rational
point = a genuine thin / higher-dimensional problem**, not geometric emptiness.
Contrast the `[3,12]` case, where the *good-open* count itself was `0` mod 5 —
that is a good-open obstruction, and the interesting content is the boundary
decomposition of §3.

## Pitfalls

- **Naive elimination of the torsion witnesses.** Eliminating the auxiliary
  variables and finding "no condition on the base" is the *expected* behavior
  of a covering condition — it is NOT evidence that every curve works. Always
  run the Rabinowitsch `w*kappa - 1` saturation to force the genuine-torsion
  coefficient nonzero; if the saturated ideal is `(1)`, the condition is a
  cover, not a subvariety (`code/agent_a2_24_d0_saturate.m`).
- **Spurious resultant branches.** Resultants pick up spurious factors — e.g.
  the `kappa=0` locus, which here gives torsion `[8]`, not `[24]`, and appears
  as a rank-1 elliptic curve `y^2 = x^3 + x^2 - 4x` that is *not* the genuine
  cover (`notes/agent_a2_24_d0_derivation.md`). Compute `g0 = GCD(P, Q)` of the
  two eliminant polynomials and divide it out **before** taking the resultant,
  then keep only the factor through your known example:

  ```magma
  g0 := GCD(P4, P3); if TotalDegree(g0) ge 1 then P4 := P4 div g0; P3 := P3 div g0; end if;
  Rt := Resultant(P4, P3, tt);
  // then select the factor with your known point on it:
  for g in Factorization(Plane) do
      if Evaluate(g[1],[1/3,2/9]) eq 0 and TotalDegree(g[1]) ge 3 then Grm := g[1]; end if;
  end for;
  ```

  (`code/agent_a2_24_d0_cover_ptsearch.m`, lines building `g0`, `Rt`, `Grm`.)
- **`EllipticCurve` / `Jacobian` on a high-degree singular plane model.**
  Magma's `EllipticCurve` / `Jacobian` on the degree-12 singular plane model of
  a genus-1 cover is pathologically slow (`notes/agent_a2_24_d0_derivation.md`:
  "Magma's EllipticCurve/Jacobian on the deg-12 model is very slow"). **Do not** call it
  to get the rank. Instead determine the rank *empirically* by direct
  rational-point search over a low-height auxiliary coordinate: fix small-height
  `m0`, solve `Grm(r, m0) = 0` for rational `r`, and also scan `r0`
  (`code/agent_a2_24_d0_cover_ptsearch.m`). Finding only a **bounded** set of
  points is the signature of **rank 0 / an isolated example**:

  ```text
  rational points on the (r,m) model, height <= 40:  ONLY 4
    (0,0),(1,0),(2,0)  -- m=0, degenerate
    (1/3, 2/9)         -- B, the sole genuine point.
  ```

  A bounded 4-point set => rank 0; `B` is essentially the only genuine simple
  `Z/24` on the `d=0` slice, so there is no positive-dimensional family to
  intersect with the 2-rank-2 locus.
- **Not stripping *all* known boundary factors before judging the core.** If
  you leave a geometric-degeneration factor in, an "irreducible core" may be an
  artifact. Multiply out the full nonboundary product (as in §1) and confirm
  the residual shares no factor with it before calling the core irreducible.
  In the `m18_416` `[4,16]` projection, the repeated degree-5 factor of
  `Res_b(E1core,E0core)` was exactly the excluded `d4=0` boundary
  (`d4 = (R-1)*4*F5`), not a new geometric branch — ignoring it would have sent
  the search onto a degenerate locus
  (`notes/m18_416_component_and_smooth_strata_2026_07_02.md`).
- **Confusing "boundary-confined" with "globally impossible."** Boundary
  analysis that confines a phenomenon to a split/nonsimple stratum proves it is
  off the *simple open* chart, not that the group is unattainable in general.
  State the narrower claim.

## See also

- `local-obstructions` — decides thin (Pattern A) vs good-open-obstructed
  (Pattern B); run it first to know whether boundary analysis (§3) is even the
  right move.
- `pell-cf-order` — exact `D_infinity` order; the covering-condition insight
  (§2) applies equally to order-`n` classes certified there.
- `contact-torsion-constructions` — where the `f = h3^2 + kappa*q3^3`
  independent-3-torsion ansatz and its `kappa` come from.
- `simplicity-certificates` — the split/nonsimple test used to certify a
  boundary lift is off the simple locus.
- `g2-torsion-lab` — the hub; read it first.
