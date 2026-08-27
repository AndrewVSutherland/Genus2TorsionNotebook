---
name: two-rank-and-factor-types
description: Compute and ENGINEER the rational 2-torsion rank of a genus-2 Jacobian from the factorization type of f — the TwoRank function, the factor-type table (sextic AND quintic accounting), points at infinity, and the W-split parametrization that makes 2-rank 2 free on the A(8) chart. WHEN a target has a [2,...] prefix ([2,24], [2,2,16], [2,2,20]), when choosing the factor type of f, when TwoRank returns a surprising value, or when a 2-rank filter is the bottleneck of a search. Trigger words - 2-rank, factor type, [1,2,2], [2,2,2], rational 2-torsion, W-split, even subsets, TwoRank, Weierstrass points.
---

# 2-rank and factor types

## When to use this

Load this when a torsion target has a `[2,...]` or `[2,2,...]` prefix and you
must (a) *compute* the rational 2-rank of a given curve, or (b) *engineer* a
family whose members have the required 2-rank built in. The 2-rank is the
number of `Z/2` invariant factors of `J(Q)_tors`; it is controlled entirely by
the factorization type of `f` over `Q`, so it is the one torsion condition
that costs nothing to impose — *if* you know the accounting rules below.

## The TwoRank function (source of truth)

Verbatim from `code/agent_a2_24_wsplit_3tors.m` (identical copies in
`code/agent_a2_24_composite8x3.m`, `code/agent_a2_24_cf_search.m`, and other
search scripts):

```magma
function TwoRank(fp)
    degs := [Degree(g[1]) : g in Factorization(fp)];
    k := #degs; even := 0;
    for mask in [0..2^k-1] do
        ss := 0; for i in [1..k] do if (mask div 2^(i-1)) mod 2 eq 1 then ss +:= degs[i]; end if; end for;
        if ss mod 2 eq 0 then even +:= 1; end if;
    end for;
    return Ilog2(even) - 1;
end function;
```

It enumerates all subsets of the irreducible factors of `f`, counts those with
even total degree, and returns `log2(#even subsets) - 1`.

**Scope: this function implements the SEXTIC accounting.** Its callers in the
search scripts all pass sextic integral models `fInt`. For a quintic model it
undercounts by 1 — see "Degree-5 models" below.

## Why it works

For `C: y^2 = f(x)` with `deg f = 6`, the Weierstrass points are the six roots
of `f`, and over `Qbar`

```text
J[2] = { even-size subsets S of the six roots } / (S ~ complement of S),
```

which has `2^4 = 16` elements (rank 4). A class is rational iff its subset is
Galois-stable, and the Galois-stable subsets are exactly the unions of the
factor orbits — i.e. subsets of the irreducible factors of `f` over `Q`. The
even-total-degree subsets form an `F_2`-vector space under symmetric
difference; the quotient by `S ~ S^c` (the full set has degree 6, even, so
complementation is available) removes one dimension. Hence

```text
2-rank = log2(#even-degree subsets of factors) - 1,
```

which is the function above. (Subsets only swapped with their complement by
Galois would need odd size 3 and cannot represent 2-torsion classes, so no
rational class is missed — the factor-orbit count is exact.)

## The closed form and the factor-type table (degree 6)

Counting even-degree subsets gives a closed form. With `k` = number of
irreducible factors:

```text
2-rank = k - 2      if at least one factor has ODD degree,
2-rank = k - 1      if ALL factors have even degree.
```

(If some factor degree is odd, exactly half of all `2^k` subsets are even;
if all are even, every subset is.) The full table for sextics:

| factor type of `f` (deg 6) | k | 2-rank |
|---|---|---|
| `[6]` | 1 | 0 |
| `[1,5]`, `[3,3]` | 2 | 0 |
| `[2,4]` | 2 | 1 |
| `[1,1,4]`, `[1,2,3]` | 3 | 1 |
| `[2,2,2]` | 3 | **2** |
| `[1,1,2,2]`, `[1,1,1,3]` | 4 | 2 |
| `[1,1,1,1,2]` | 5 | 3 |
| `[1,1,1,1,1,1]` | 6 | 4 |

Anchors used across the lab: the A(8) W-split locus is `[2,2,2]` (2-rank 2,
see below); the `[2,24]` searches require `TwoRank(fInt) ge 2`
(`code/agent_a2_24_wsplit_3tors.m`); the CF-search family
`(x^2-1)(x^2+ax+b)(x^2+cx+d)` is `[1,1,2,2]` (2-rank 2,
`code/agent_a2_24_cf_search.m` header).

## Degree-5 models: the point at infinity is an extra factor

A quintic model has ONE Weierstrass point at infinity, and it is rational. The
correct accounting appends an extra degree-1 orbit to the factor multiset:

```text
quintic 2-rank = (k+1) - 2 = k - 1     (the appended 1 makes an odd degree present)
```

| factor type of `f` (deg 5) | effective multiset | 2-rank |
|---|---|---|
| `[5]` | `[5,1]` | 0 |
| `[1,4]`, `[2,3]` | `[.,.,1]` | 1 |
| `[1,2,2]` | `[1,2,2,1]` | **2** |
| `[1,1,3]` | `[1,1,3,1]` | 2 |
| `[1,1,1,2]` | `[1,1,1,2,1]` | 3 |
| `[1,1,1,1,1]` | `[1,1,1,1,1,1]` | 4 |

The load-bearing anchor: the contact-6 chart's `f = x*((b+3)x^2+(a-3)x+2)*
(2x^2+(b-3)x+(a+3))` (a quintic, `notes/contact6_m36.md`) has factor type
`[1,2,2]` and 2-rank **2** when both quadratics are irreducible — the locus
that carried the simple `[6,6]` example. Feeding this quintic to the coded
`TwoRank` returns 1, not 2: the function's sextic accounting misses the
infinity orbit. **When working on a quintic model, add 1 (equivalently, append
a `1` to the multiset) or pass an even-degree model.**

## Points at infinity of a sextic model

The two points at infinity of `y^2 = f(x)`, `deg f = 6`, are rational iff
`LeadingCoefficient(f)` is a square. This is exactly the

```magma
if IsSquare(LeadingCoefficient(fp)) then cnt +:= 2; end if;
```

term in the `CountCurve` point-count used by the simplicity certificate
(`code/agent_a2_24_composite8x3.m`); forgetting it corrupts `chi(T)` — see
`simplicity-certificates`. (For the 2-rank itself the infinity pair of a
sextic contributes no extra class: the six roots already generate `J[2]`.)

## Engineering 2-rank: the W-split parametrization (flagship result)

Blind scanning for 2-rank on a chart is hopeless: on the A(8) chart the
2-rank-2 locus has measured density about `3e-5` (31 per million) in a blind
`(r,p,t)` scan (`notes/agent_a2_24_composite.md`). The fix is to PARAMETRIZE
the factor-type condition so it holds identically.

On A(8), `f = q*(Q^2 + q)` with `q` a quadratic. Write `W := Q^2 + q`
(quartic). Factor type `[2,2,2]` (2-rank 2) holds when `W` splits into two
rational quadratics. Imposing the split ansatz reduces to `E = 0` with

```text
E(p) = beta^2*(2*d + a + beta^2)^2 - 4*(d^2 + c)*beta^2 - b^2,
```

which for fixed `(r, t, beta)` is **quadratic in the chart parameter `p`** —
so every rational root of `E` hands you a 2-rank-2 curve. Verbatim from
`code/agent_a2_24_wsplit_3tors.m`:

```magma
// W-split condition as a quadratic in p; return rational p-roots
function SplitWPvals(rv, tv, bv)
    Pp<pp> := PolynomialRing(Q);
    e := tv^2 - 2*pp*tv/rv; d := e + 2*pp - rv^2; lambda := rv/tv;
    a := rv^2 - lambda;
    b := 2*rv*pp - 2*lambda*(pp + rv*tv) + 2*rv*lambda;
    c := pp^2 + 2*pp*rv^2 - rv^4 - rv^3*tv - rv*pp^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pp + rv^2*tv - 3*pp*tv + rv*tv^2);
    E := bv^2*(2*d + a + bv^2)^2 - 4*(d^2 + c)*bv^2 - b^2;
    En := Numerator(E);
    if En eq 0 or Degree(En) lt 1 then return []; end if;
    return [rt[1] : rt in Roots(En)];
end function;
```

Measured impact (`notes/agent_a2_24_composite.md`): ~90% of W-split-
parametrized curves are genuinely 2-rank 2, vs `3e-5` blind — a factor of
~30,000 in search efficiency. In the H=26 production run, ~1M genuinely
2-rank-2 order-8 curves were produced this way. **This is the template move:
when a codimension-1 condition is the bottleneck, find the sub-parameter in
which it is low-degree and parametrize it away.** (The searches still call
`TwoRank(fInt) ge 2` afterwards as a cheap confirmation — parametrized
conditions can degenerate at special parameter values.)

## Pitfalls

- **Applying the coded `TwoRank` to a quintic model.** It implements the
  sextic accounting and undercounts a degree-5 model by exactly 1 (it misses
  the rational Weierstrass point at infinity). The contact-6 `[1,2,2]`
  quintic has 2-rank 2, not 1. Add 1 for quintics, or homogenize your reasoning
  to the table above.
- **Forcing extra rational roots to raise 2-rank.** Factor types with many
  linear factors (`[1,1,1,2]`, `[1,1,1,1,2]`, ...) push you toward
  split/nonsimple Jacobians: the contact-6 extra-root `[6,6]` curves were ALL
  nonsimple, while the irreducible-quadratics `[1,2,2]` locus carried the
  simple example (`notes/contact6_m36.md`; also `target-playbook`). Prefer
  irreducible quadratic factors — `[2,2,2]` over `[1,1,2,2]` over
  `[1,1,1,1,2]` — when the target allows.
- **Even sextics.** If the coefficients of `x^1, x^3, x^5` all vanish, the
  curve is bielliptic and `J` splits — worthless here regardless of 2-rank.
  The searches skip them: `&and[Coefficient(f,i) eq 0 : i in [1,3,5]]`
  (`code/agent_a2_24_composite8x3.m`). Do not "fix" a search by removing this
  guard.
- **The Factorization unit trap (cross-ref).** Magma's `Factorization` returns
  monic factors and silently drops the leading unit. For *counting degrees*
  (this skill) that is harmless; for *reconstructing* `f` or a quotient curve
  from the factors it has caused real bugs — see `magma-lab-conventions`.
- **2-rank is necessary, not sufficient, for the `[2,...]` prefix.** 2-rank 2
  guarantees two `Z/2` invariant factors ONLY if the rest of the torsion
  cooperates; the final arbiter is `TorsionSubgroup(J)`'s invariant factors.
  A `[24]` curve with 2-rank 1 is not `[2,24]`; see
  `validate-and-record-a-hit`.
- **Degenerate parametrized members.** On the W-split locus, special
  `(r,t,beta)` values collapse `q` or `W_i` (discriminant 0, degree drops).
  Keep the standard degenerate skips (`Degree(f) lt 5`, `Discriminant(f) eq 0`)
  before trusting the parametrization — see `running-torsion-searches`.

## See also

- `g2-torsion-lab` — hub; move (4) "2-rank control" in the methodology.
- `target-playbook` — move (d): choosing the factor type for a target's
  `[2,...]` prefix.
- `contact-torsion-constructions` — the contact-6 factorization whose
  `[1,2,2]` locus is the worked simple-`[6,6]` example.
- `simplicity-certificates` — the `CountCurve` infinity term; why extra-root
  loci produce split curves you must reject.
- `running-torsion-searches` — where the 2-rank filter sits in the funnel.
- `magma-lab-conventions` — the Factorization unit trap.
