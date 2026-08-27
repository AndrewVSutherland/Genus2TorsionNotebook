---
name: named-charts-reference
description: Reference formulas for the normalized parametric charts that carry a marked high-order rational class for free — A(8) (order 8), A(12) (order 12), M(2,12) (order 12), M_1(8,4) ([4,16]). WHEN working on order-8/12/16/24 targets ([2,24], Z/24, [3,12], [4,16]) and you need the exact chart polynomials and the visible marked class.
---

# Named-charts reference

## When to use this

A quick, cited lookup of the four normalized charts that come with a marked
high-order rational class built in, so you can specialize/search without
re-deriving the parametrization. Each entry gives the exact chart polynomials,
the parameters, the visible marked class, and the file that grounds it. Where
an exact formula is not reproduced here, the entry says which file to read.

Charts covered:

| Chart      | Marked order | Params   | Ground file |
|------------|--------------|----------|-------------|
| A(8)       | 8            | `r,p,t`  | `code/agent_a2_24_composite8x3.m` |
| A(12)      | 12           | `p,z,r`  | `code/agent_a12_224_descent_setup.m` |
| M(2,12)    | 12           | `z,r`    | `notes/contact6_m36.md` |
| M_1(8,4)   | [4,16] target| `R,w`    | `code/agent_m18_416_search_crt.m` |

---

## A(8) — order-8 chart (params `r`, `p`, `t`)

Source: function `A8f` in `code/agent_a2_24_composite8x3.m` (m=1 gauge, from
`notes.tex` / `search_A8_*`). Copy the formulas VERBATIM:

```text
e      = t^2 - 2*p*t/r;
d      = e + 2*p - r^2;
lambda = r/t;

u = p + r*t - 2*r;
v = e + r^2 - r*p - r^2*t + 3*p*t - r*t^2;

a = r^2 - lambda;
b = 2*r*p - 2*lambda*(p + r*t) + 2*r*lambda;
c = p^2 + 2*p*r^2 - r^4 - r^3*t - r*p^2/t
    - lambda*(r^2 + e) + 2*lambda*(r*p + r^2*t - 3*p*t + r*t^2);

Q = x^2 + d;
q = a*x^2 + b*x + c;
f = q*(Q^2 + q);
```

**Visible order-8 class.** `g8 = x^2 + u*x + v`, with the (base) `y`-part

```text
L = r*x + (p - r^2),
ellBase = -(q + Q*L),
v8 = (-Lden*ellBase) mod g8,
```

where `Lden` is the denominator-clearing scalar `L` from `IntModel` (in the
script, `fInt, L := IntModel(f)` and then `v8 := (-L*ellBase) mod g8`). The
marked class is `D8 = J![g8, v8]`; the script asserts `8*D8 = O` and
`4*D8 != O` (exact order 8). Note: the identifier `u`/`v` here are the two
coefficients of `g8`; `L`/`ellBase` build its `y`-part. Do not confuse them
with the halving `u, ell`.

**Coprime use.** A(8) carries order 8 for free; combine with any independent
rational 3-torsion (see `contact-torsion-constructions`) to get `Z/24 = 8 x 3`
(`code/agent_a2_24_composite8x3.m`), and with 2-rank >= 2, `[2,24]`.

**The `d = 0` slice** (from `notes/agent_a2_24_d0_derivation.md`): imposing

```text
d = e + 2*p - r^2 = 0   <=>   p = r*(r+t)/2
```

gives `Q = x^2` (so `f = q*(x^4 + q)`) and makes `x = 0` a free rational point
(`f(0) = q(0)^2 = c^2`). Curve B (`r=1/3, p=-1/9, t=-1`), one of the two known
simple `Z/24`, lives on this slice. This slice is tractable but the genuine
`Z/24` cover on it is genus-1 rank-0 — see the note and
`contact-torsion-constructions` for why `[2,24]` resists this route.

---

## A(12) — order-12 chart (params `p`, `z`, `r`)

Source: `code/agent_a12_224_descent_setup.m` (header + symbolic build,
asserted over `Q(p,z,r)`); also used in `code/agent_a12_224_funnel.m`
(function `ChartData`). Structure is `f = R*F` with `F = Q^2 + R*ell^2`. Copy
VERBATIM:

```text
s      = (z^2 - 4*p^2 + 1)/(2*z);
t      = (z^2 + 4*p^2 - 1)^2/(8*p^2*z);
mu     = ((s^2 - 1)*(2*p*r + 1) - p^2*(2*s*t - 4))/(4*p^3);
lambda = (4 - mu^2)*p^2/(s^2 - 1);

T1  = p*x + r;
R   = (T1^2 + x - 1)/lambda;
ell = s*x + t;
Q   = 2*T1 + mu*R;
F   = R*x^2 + 4*(R + x - 1)*(R - 1);      // = Q^2 + R*ell^2  (asserted)
f   = R*F;                                 // degree 6
```

**Visible marked class of order 12.** From the same file:

```text
P4  = [ Q_monic,        (R*ell) mod Q_monic       ];
P6  = [ (R + x - 1)_monic, (x*R) mod (R+x-1)_monic ];
P12 = P4 + P6,   of order 12.
```

(In `agent_a12_224_funnel.m`: `u4 = Q/lc(Q)`, `u6 = (R+x-1)/lc(R+x-1)`,
`P4 = J![u4, (L*R*ell) mod u4]`, `P6 = J![u6, (L*x*R) mod u6]` with `L` the
integral-model scalar; then `P12 = P4 + P6`.)

**`[2,24]` use.** Target `[2,24]` = some 2-torsion translate `D = P12 + T`
divisible by 2, plus 2-rank >= 2. The funnel
(`code/agent_a12_224_funnel.m`) restricts to the conic locus
`r = (v^2 - p^2)/p` (params become `(p, z, v)`) and runs a symbolic 2-descent
sieve (`N1 == N2 == p*(p+r) = v^2` identically) before exact
`IsDivisibleBy`/`TorsionSubgroup`. Read `notes/agent_a12_224_descent.md` and
the `agent_a12_224_*` scripts for the descent details.

---

## M(2,12) — order-12 chart (params `z`, `r`)

Source: `notes/contact6_m36.md`. Cross-listed in `halving-and-doubling` (it is
the halved-contact-6 chart). Copy VERBATIM:

```text
m = (1 - z^2)/(4*(r + 1));
T = m*x^2 - x + r;
Y^2 = (x - r)^2*(T + 1)^2 + 4*m*x^2*T*(T + 1);
```

The marked class has **order 12**. This is the algebraic halving cover of the
contact-6 marked order-6 class, in the centered coordinate. Add an independent
rational 3-torsion direction to target `[3,12]`
(`code/contact6_m36_halveD_m312_search.m`). Beware the halving wall at `p=5`
(see `halving-and-doubling`): the good open chart has `allowed_312 = 0` mod 5,
and the small `[3,12]` points found are split/nonsimple.

---

## M_1(8,4) — the `[4,16]` chart (params `R`, `w`)

Source: function `FamilyPolynomial(R, w)` in
`code/agent_m18_416_search_crt.m`; overview in
`notes/m18_416_component_and_smooth_strata_2026_07_02.md` and
`notes/agent_m18_416_search_notes.md`. The curve is `C: y^2 = x*A*B`. Copy
VERBATIM:

```text
t = (2*R^2 + (1 - w^2)*R - 2*w^2)/(4*(w^2 - 1));

A = x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B = (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);

f = x*A*B;
```

Auxiliary data used in the `[4,16]` search (from the same script), for the
second-halving covers:

```text
Y_R = -2*R*(R-1)^2*Qfac/(w^2 - 1),   Qfac = R^2 - (1/2)*R*w^2 + (1/2)*R - w^2;

PlusDisc  = -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)*(R*w - 3*R + 3*w - 1);
MinusDisc =  4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
              *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2);
```

`FirstCoverPossible(R, w)` is `IsSquareQ(PlusDisc) or IsSquareQ(MinusDisc)`;
the two-torsion classes halved are `T_x = [x, 0]` and `P_R = (-R, Y_R)`. The
target `[4,16]` requires halving both, and the exact second halving of `P_R`
is where the wall bites (through height 150, `P_R` never halves on the six
smooth `Q_7` strata). For the search machinery — the `p=7` boundary closure,
the mod-343 gate, aux primes `11,13,37,41`, and the tangent sieve
`17,23,29,47` — read `notes/m18_416_component_and_smooth_strata_2026_07_02.md`
and `notes/agent_m18_416_search_notes.md`; the live search is
`code/agent_m18_416_live_stratum_search.m`.

**Not reproduced here** (read the file): the explicit second-halving surface
equations `E1core, E0core`, the `p=7` blowup/level-2 lift, and the
tangent-candidate congruence — these live in
`code/agent_m18_416_p7_blowup_level2.m`,
`code/agent_m18_416_e1e0_reduction.m`, and the `search_crt` /
`live_stratum_search` scripts.

---

## Pitfalls

- **Overloaded `u`, `v`, `L`, `Q`, `t`, `r` across charts.** In A(8), `u, v`
  are the coefficients of the marked `g8` and `L = r*x + (p-r^2)`; in A(12),
  `t` is a chart parameter *and* an internal quantity `t =
  (z^2+4p^2-1)^2/(8p^2 z)`; in M_1(8,4), `t` is the internal
  `(2R^2+(1-w^2)R-2w^2)/(4(w^2-1))`, not a free parameter. Always copy the
  whole block for one chart; never mix a symbol's meaning across charts.
- **Two different order-12 charts.** M(2,12) (`m, T, Y` in `(z,r)`, the
  halved-contact-6 chart) and A(12) (`f = R*F`, `F = Q^2 + R*ell^2`, in
  `(p,z,r)`) are *different* families. M(2,12) is for `[3,12]`; A(12) is the
  `[2,24]` descent chart. Do not swap their formulas.
- **Integral model scalar.** The visible marked classes use a
  denominator-clearing scalar `L` (from `IntModel`/`IntegralModelPolynomial`)
  in the `y`-part, e.g. A(8) `v8 = (-L*ellBase) mod g8`, A(12)
  `P4 = J![u4, (L*R*ell) mod u4]`. Forgetting `L` gives a class on the wrong
  (non-integral) model. Copy the scalar handling from the script.
- **Degenerate/pole loci.** Each chart has parameter values that break it:
  A(8) needs `t != 0`, `r != 0`, `r*t != 1` (lambda pole / `a = 0`); A(12)
  needs `p, z != 0`, `s^2 != 1`, `lambda != 0`; M(2,12) needs `r != -1`;
  M_1(8,4) needs `w^2 != 1`. Guard these (the scripts do) or you get spurious
  degenerate curves.
- **M_1(8,4) is not a "free order-16" chart.** It gives `x*A*B` with the
  `[4,16]` *target*; the order-16 class is obtained only after halving both
  `T_x` and `P_R`, and the second halving is obstructed. Do not treat the raw
  chart as carrying a marked order-16 class the way A(8)/A(12) carry 8/12.
- **Cite, don't paraphrase.** These formulas are load-bearing. When you copy
  one into a script or the paper, keep the exact expression and cite the file;
  a transcription slip in `mu` or `lambda` silently produces a curve with the
  wrong torsion.

## See also

- `g2-torsion-lab` — hub.
- `contact-torsion-constructions` — the contact primitives behind these charts,
  and the independent 3-torsion to compose with A(8)/A(12) for `Z/24`/`[2,24]`.
- `halving-and-doubling` — M(2,12) is cross-listed there; the halving wall that
  governs M_1(8,4)'s second halving and the M(2,12) `[3,12]` boundary.
- `running-torsion-searches`, `finite-prefilters` — the funnel that turns a
  chart into examples.
- `simplicity-certificates` — required before claiming any chart hit is a
  result.
