# Lane 2 (`claude_ov_c7z`) — the order-112 decision on the contact-7 chart

**Date:** 2026-07-25 · **Machine:** claude-box (32 vCPU), all Magma via `code/claude_magma_slot.sh`.
Earlier stages of this lane also used aws-spot-5 (192-way, C sieve); both spot boxes were stopped
on 2026-07-25 and nothing here needs them.

**Target.** `[2,2,2,14]`, order **112** — would be a record for a geometrically simple genus-2
Jacobian over Q. Route: the **contact-7 chart** (the same chart that produced the six generic
`[2,2,14]` curves of commits `39b71da` / `683fff7`).

---

## 0. The setup, re-derived from scratch

All of the following is re-derived symbolically in `code/claude_ov_c7z_xdecide.m` PART 0
(nothing is inherited on trust from earlier sessions); log `results/claude_ov_c7z_xdecide.log`.

Normalized contact-7 chart:

```
h  =  1 - (7/2) x + A x^2 + B x^3 ,        f = ( h^2 + (x-1)^7 ) / x^2      (deg f = 5)
```

`y^2 = f(x)` carries a marked rational class of order 7 by contact. A **rational root** of `f`
is `r = 1 - v^2` where `v` is a root of the quintic

```
Q(v) = (v+1)^2 (c4 v^2 + c0) + v^5 - v^3 - v^2/2 ,     c4 = B+2 ,  c0 = 5/2 - (A+B).
```

Identically in `(c4,c0)`: `Q(-1) = -1/2`, `Q'(-1) = 3`, `Q'(0) = 2 Q(0)`. Substituting
`z := 1/(1+v)` and clearing denominators,

```
P(z) = prod_i (z - z_i) = z^5 - 6 z^4 + e2 z^3 - e3 z^2 + e4 z - 2
```

with the **universal relations** (verified by `assert` in the script)

```
e1 = 6 ,      e5 = 2 ,      e3 = 2 e4 - 2                       (2 free parameters e2, e4)
```

equivalently `sum z_i = 6`, `prod z_i = 2`, `sum 1/(1 - z_i) = 3` (the form used by the C sieve;
the two forms are equivalent because `sum 1/(1-z_i) = P'(1)/P(1)` and `P'(1) - 3P(1) = e3 - 2e4 + 2`).

`k` = number of rational `z_i` ⇔ factor type of the quintic ⇔ rational 2-rank.
**`k = 5` ⇔ torsion ⊇ `[2,2,2,14]`, order 112.** `k = 3` is the "three-root" locus that produced
the known `[2,2,14]` curves.

Forbidden (degenerate) `z`-values: `z = 1` (`v = 0`, `r = 1` a double root), `z = 1/2` (`v = 1`,
`r = 0`, divided out by the `x^2`), `z = 0` (`v = ∞`).

---

## 1. The order-112 surface `S`

`code/claude_ov_c7z_surface.m` PART A → `results/claude_ov_c7z_surface.log`.

`S = { e1 = 6, e5 = 2, e3 = 2 e4 - 2 } ⊂ A^4` (after eliminating `z5 = 6 - z1 - ... - z4`).
**`S(Q)` is exactly the set of order-112 configurations.** Measured:

| invariant | value |
|---|---|
| `dim S` | 2 |
| `deg S` (projective closure in `P^4`) | 20 |
| irreducible / reduced | true / true |
| projective closure | complete intersection of degrees **4 and 5** in `P^4` |
| `ω_S` | `O_S(4)` (Gorenstein — unconditional, holds regardless of singularities) |
| `dim Sing(S)` (projective) | 1 — 25 components (10 of dim 1, 15 of dim 0), all of degree 1 |
| `Sing(S) ⊂ {t=0}` ? | **false** |
| boundary `S ∩ {t=0}` | dim 1, degree 10, 10 components |
| affine `S` | dim 2, irreducible, reduced, `Sing` of dim 0 and degree 5 |
| (numerics of a *smooth* (4,5) c.i., for reference only) | `K^2 = 320`, `p_g = 69`, `q = 0` |

So `S` is a canonically-polarized-looking surface (`ω_S = O_S(4)` ample) — morally of general type,
hence Bombieri–Lang would predict `S(Q)` non-dense. **This is not used as an input to anything
below**; it is recorded because it says "expect finitely many, expect them to be sparse", and
because `dim Sing = 1` means the smooth-c.i. numbers above are *not* literally the invariants of `S`.

---

## 2. The slice fibration `z_1 = c`

Fix `z_1 = c`. The other four `z` are the roots of the pencil

```
P_m(w) = w^4 - e1 w^3 + m w^2 - (A m + B) w + e4 ,
e1 = 6 - c ,   e4 = 2/c ,   T = 3 - 1/(1-c) ,
A = (T-2)/(T-1) ,           B = ( T(c - 5 + 2/c) - 3c + 14 ) / (T-1).
```

Three curves sit over each slice:

* `Cab_c` — the quartic splits into two **rational quadratics**
  `(w^2 - a w + b)(w^2 - (e1-a) w + e4/b)`:
  `Cab_c : (e1-A-a) b^2 + (A a^2 - A e1 a - B) b + e4 (a - A) = 0` — a conic in `b` over the
  `a`-line, i.e. `y^2 = Δ_c(a) := (A a^2 - A e1 a - B)^2 - 4 (e1 - A - a) e4 (a - A)`.
* the **three-root cover** `u^2 = d1 = a^2 - 4b` over `Cab_c` (3 rational `z` ⇒ `[2,2,14]`-type);
* the **order-112 fibre** `X_c : u^2 = d1, v^2 = d2 = (e1-a)^2 - 4 e4/b` over `Cab_c`.

### 2.1 Genus table — measured, `code/claude_ov_c7z_surface.m` PART B

| slice | `genus(Cab_c)` | `genus(3-root cover)` | `genus(X_c)` |
|---|---|---|---|
| **`c = 2`** | **0** | **1** | **4** |
| every other `c` tested (18 of them) | 1 | 4 | **13** |

Tested `c ∈ {-1, 3, 4, 5, -2, 3/2, -1/3, 5/2, 7, 2/3, -5, 10, 6, -1/2, -1/9, -7/3, 22/7, 49/36}`.

### 2.2 `c = 2` is the **only** special slice — a symbolic proof, not a scan

`code/claude_ov_c7z_specialc.m` → `results/claude_ov_c7z_specialc.log`.
`disc_w(P_m)` has degree 5 in `m`; its discriminant in `m` (the "double discriminant") is a
rational function of `c` whose **numerator has degree 58** and factors as

```
(C-2)^15 · (C-1)^15 · (C^2 - 6C + 4)^2 · (C^3 - 12C^2 + 20C - 8)^2
        · (C^6 - 17C^5 + 93C^4 - 167C^3 + (3622/27) C^2 - (1388/27) C + 8)^3
```

**Its only rational roots are `C = 1` and `C = 2`**, and `C = 1` is a forbidden `z`-value.
Hence for **every** rational `c ∉ {0, 1/2, 1, 2}` the branch divisor of `w ↦ m` is non-degenerate
and the three-root fibre is the smooth genus-4 `(3,3)` curve.

Corroborated by brute force twice:
* `results/claude_ov_c7z_slicegenus.log` — **464** values of `c`: every one irreducible of
  genus 4, except `c = 2` (genus 1).
* `results/claude_ov_c7z_cscan.log` — **429** values of `c`: only `c = 2` flagged LOWGENUS,
  with factor type `[<3,3,1,1>]`.

**Consequence (Faltings).** For every rational `c ≠ 2` there are only *finitely many* three-root
configurations with `z_1 = c`. `c = 2` is the **unique** slice with infinitely many — its
three-root curve is `y^2 = x^3 - 12x`, conductor 288, `MW = Z/2 × Z` (rank 1, proven).
So the slice this lane decided is exactly the one slice that could have carried an infinite family.

---

## 3. MAIN RESULT — the `c = 2` slice is decided

Scripts `claude_ov_c7z_xdecide.m`, `claude_ov_c7z_chabauty.m`, `claude_ov_c7z_verify.m`,
`claude_ov_c7z_gap.m`; logs of the same names. Commits `63c7044`, `83bad19`, `0a929a1`, and this one.

For `c = 2`: `e1 = 4`, `e4 = 1`, `T = 4`, `A = 2/3`, `B = 0`, pencil
`P_m(w) = w^4 - 4w^3 + m w^2 - (2/3) m w + 1`.

```
Cab_2 :  -3ab^2 + 10b^2 + 2a^2b - 8ab + 3a - 2 = 0        plane cubic, GENUS 0
```

Parametrized from the point `(-22/9, 7/27)`:

```
a(T) = (-22/9 T^2 + 80/3 T - 40)/(T^2 + 4T - 12)     poles T = -6, 2
b(T) = ( 7/27 T^2 + 10/9 T -  8/3)/(T^2 - 2T)        poles T = 0, 2 ; zeros T = -6, 12/7
```

Square classes of `d1, d2` (extracted with the unit-safe `SqClass`, avoiding the
`Factorization`-drops-the-unit trap):

```
G1 = T (T-3) (T^2 - 24T + 36)                  [three-root cover  u^2 = G1]
G2 = (T-3) (7T-12) (37T^2 - 96T + 36)          [the other quadratic splits]
gcd(G1, G2) = T - 3
```

```
X : u^2 = G1(T),  v^2 = G2(T)          GENUS 4   (function-field tower, deg 4 over Q(T))
```

Its three `(Z/2)^2`-quotients:

| quotient | curve | data |
|---|---|---|
| `E1 : u^2 = G1` | `y^2 = x^3 - 12x` | conductor 288, `MW = [2,0]`, **proven** |
| `E2 : v^2 = G2` | `y^2 = x^3 - 12x` | conductor 288, `MW = [2,0]`, **proven** |
| `H : y^2 = G1 G2/(T-3)^2` | genus 2 | see below |

```
H : y^2 = 259 T^6 - 7332 T^5 + 37512 T^4 - 74304 T^3 + 60912 T^2 - 15552 T
      = 259 · T (T - 12/7) (T^2 - 24T + 36) (T^2 - 96/37 T + 36/37)
disc(H)         = 2^54 · 3^46
Jac(H) torsion  = [2, 6]
RankBounds      = 0 .. 1                                    [0.840 s]
```

**Rank is exactly 1**: the explicit class `P = (x^2 - 2x, 8x, 2)` has canonical height
`ĥ(P) = 0.50118239204717836371111998436553 > 0`. `P` is **saturated at every prime ≤ 50**
(same generator, same height, index 1). `HeightConstant(J) = 22.2342101128697832838282821298`.

**Chabauty.** `Chabauty(P)` returns **8 points**, second return value (the index-bound set)
**empty**; the answer is *identical* to `Points(H : Bound := 200000)`. Because `rank = 1`,
Magma's `Chabauty(P)` semantics ("points whose image lies in the **saturation** of `⟨P⟩`") give
**all** of `H(Q)`:

```
H(Q) = { (0:0:1), (12:0:7), (±): (2:∓16:1), (-6:∓11664:1), (6:∓3888:5) }     — 8 points
T-coordinates:  { -6, 0, 6/5, 12/7, 2 }
```

`lc(G_H) = 259` is not a square ⇒ **no rational point at infinity** on `H`; on `X` the fibre over
`T = ∞` needs `√259` ⇒ residue field `Q(√259)`, **not rational**.

**Pull back to `X`.** For `T ≠ 3`, `y := uv/(T-3)` is a rational point of `H`, so `X(Q) ⊆` the
preimage of `H(Q)` together with the `T = 3` fibre (where `G1G2 = (T-3)^2 G_H` loses the factor and
the map is undefined). Checked (`claude_ov_c7z_gap.m`):

* `T = 6/5`: `G1 = -11664/625 < 0`, `G2 = -104976/625 < 0` — no lift.
* `T = 3`: `G1/(T-3)|_{T=3} = -81` (not a square) — the two places of `X` over `T = 3` have
  residue field `Q(√-1)`, **irrational**; and `G_H(3) = -59049 < 0`, so `H` has no rational point
  over `T = 3` either. The affine model's singular point `(T,u,v) = (3,0,0)` is the only rational
  point of the affine `X` that does not map to a rational point of `H`, and it is tested explicitly.

```
X(Q)  (affine model, 13 points)   T ∈ { -6, 0, 12/7, 2, 3 }
```

Independent naive recounts agree exactly: `|p|, q ≤ 1200` and `|p|, q ≤ 3000` both give
`T ∈ {-6, 0, 2, 3, 12/7}`.

**Every one of the 13 is degenerate:**

| `T` | why it is not a configuration |
|---|---|
| `-6` | pole of `a(T)` — a point at infinity of `Cab_2`, `a = ∞` |
| `0` | pole of `b(T)` — `b = ∞`, forcing a root `w = 0`, contradicts `e4 = 1` |
| `2` | pole of both `a(T)` and `b(T)` |
| `12/7` | `b = 0` ⇒ some root `w_i = 0` ⇒ `e4 = 0`, contradicts `e4 = 1` |
| `3` | `(a,b) = (2,1)`, quartic `= (w-1)^4`; the five `z` would be `(2,1,1,1,1)` — repeated roots, so `f` is not squarefree and there is no genus-2 curve |

### 3.1 The one gap in the earlier write-up, now closed

A rational parametrization `P^1 → Cab_2` reaches every rational *place*, but a rational **singular
point** of the plane model with conjugate branches would be a rational point of `Cab_2` *not* in the
image. `claude_ov_c7z_gap.m` measures this: `Sing(Cab_2)` is **one point**, `(a:b:1) = (2:1:1)`,
with **two places, both of residue degree 1** — a node with rational branches, hence in the image
(it is exactly the `T = 3` point). `d1 = d2 = 0` there and the quartic is `(w-1)^4`. **No rational
point of `Cab_2` is missed.**

### 3.2 Statement of the theorem — and its exact scope

> **Theorem (Lane 2, unconditional).** No order-112 configuration on the contact-7 chart has
> `z = 2` among its five `z`-values. Equivalently, no genus-2 curve in the normalized contact-7
> chart with torsion containing `[2,2,2,14]` has the rational Weierstrass point `r = 3/4` arising
> from the quintic root `v = -1/2`.
>
> Since `c = 2` is (§2.2) the **unique** rational slice whose three-root fibre is not of genus 4,
> this is exactly the statement: *the one slice of the contact-7 order-112 surface that carries
> infinitely many `[2,2,14]`-type configurations carries no order-112 configuration at all.*

§4 below adds a **second** decided slice by a completely different (rank-0) argument:

> **Corollary.** No order-112 configuration on the contact-7 chart has `2` **or** `2/3` among its
> five `z`-values. These are the two `z`-values `1/(1±v)` sitting over the single quintic root
> `v = ±1/2`, i.e. over the single Weierstrass point `r = 3/4`.

This is answer **(a)** in the orchestrator's taxonomy — two slices out of infinitely many —
with the (c)-flavoured addendum of §2.2: they are two slices, but `c = 2` is the *only
geometrically distinguished* one, so the decided part is exactly the part that could have carried
an infinite family. It is emphatically **not** (b): the contact-7 three-root surface as a whole is
**not** decided, and nothing here says order 112 is impossible on the contact-7 chart.

---

## 4. Second decided slice: `c = 2/3` (new, this session)

`code/claude_ov_c7z_cabrank.m` → `results/claude_ov_c7z_cabrank.log`,
`results/claude_ov_c7z_cabrank_c23.log`.

For 83 slices `c = n/d` (`d ≤ 4`, `|n| ≤ 12`, plus the 27 `z`-values occurring in the eleven
recorded three-root hits) I put `Cab_c` in the model `y^2 = Δ_c(a)`, built the elliptic curve and
computed `TorsionSubgroup` and `RankBounds`. **71 of the 83 completed** in 38 m 55 s; the run was
then stopped — the 12 missing slices are the huge-height `z`-values from the recorded hits
(`-61/450, 625/114, -41/124, -16/17, 289/124, 297/133, 361/525, 50/33, 189/155, 121/155, 49/36,
50/63`), on which `RankBounds` enters an unconditional class-group proof for an 11-digit
discriminant and does not return. The log records this explicitly and gives the per-slice resume
command. Result over the 71 completed slices (see §6 for the exact census):

* `c = 2` — `Δ_2 = (4/9)(a-2)^2 (a^2 - 4a + 5)`, **not squarefree** ⇒ `Cab_2` rational. (Control.)
* **`c = 2/3`** — `E = [0,1,0,16,180]`, **conductor 48**, torsion `Z/8`, **rank 0 .. 0**.
  `#E(Q) = 8`; naive search on `C` finds 8 points, so `C(Q)` is complete; every resulting
  `(a,b)` was tested and **none** has `d1` and `d2` both squares. **Slice decided.**
* every other slice: `Cab_c` is an elliptic curve of **rank ≥ 1** (mostly 1, several 2, a few 3).

`c = 2/3` is not an accident: `z = 2 ⇔ v = -1/2` and `z = 2/3 ⇔ v = +1/2` are the two `z`-values
sitting over the *same* root `r = 1 - v^2 = 3/4` of `f` (the involution `z ↦ z/(2z-1)`).
So the two decided slices are the `v ↔ -v` pair over `r = 3/4`.

**Why no third slice is decidable this way.** Every other `Cab_c` scanned has positive rank, so
`Cab_c(Q)` is infinite and the slice cannot be finished at the `Cab` level. One has to go up to the
genus-4 three-root cover, or to `X_c` of genus 13 — and Magma has no Chabauty for genus 4.

---

## 4b. Negative, measured: there is **no** small-prime local obstruction

`code/claude_ov_c7z_localobs.m` → `results/claude_ov_c7z_localobs.log`.

If a slice fibre had no `F_p`-point at a prime of good reduction it would be decided outright, for
free. Weil says this is only possible at small `p`: a genus-`g` curve over `F_p` certainly has a
degree-1 place once `p + 1 > 2g√p`, i.e. for `p ≥ 59` when `g = 4` and `p ≥ 683` when `g = 13`.
So I counted degree-1 places of the **smooth models** (function-field tower, which handles infinity
and the singularities of the plane model correctly) of

* `D_c` — the three-root cover `u^2 = d1` over `Cab_c` (genus 4, genus 1 at `c = 2`), and
* `X_c` — the order-112 fibre `u^2 = d1, v^2 = d2` (genus 13, genus 4 at `c = 2`)

over every prime `p < 100` at which the genus is preserved (a necessary condition for good
reduction; primes where the genus drops, where the plane curve becomes reducible or not
geometrically irreducible, or where a cover degenerates, are discarded).

```
31 slices  ×  17–23 genus-preserving primes each  =  594 (c, p) pairs
EMPTY counts: ZERO — every single one of the 594 has a degree-1 place on BOTH D_c and X_c.
```

Slices tested: `2, ±1, ±2, ±3, ±4, ±5, ±6, 7, 8, 9, 10, 11, 12, ±1/2, ±1/3, ±2/3, 3/2, 5/2, 4/3,
5/3, 7/3, 1/4, 3/4, 5/4`.

> **Conclusion.** Order 112 on the contact-7 chart is **thin, not locally obstructed**, at least at
> every prime below 100 on every slice tested. Nobody should look for a congruence proof of
> impossibility here; the obstruction, if any, is global.

---

## 5. The superseded sieve — what it actually proved

`code/claude_ov_c7z_splitsieve.c`, log `results/claude_ov_c7z_splitsieve_H3200.log`
(aws-spot-5, 192 threads, stopped 2026-07-25 ~12:20 UTC after 27 351 s).

Formulation: given a pair `(z1, z2)`, the residual cubic `t^3 - At^2 + Bt - C` is **forced**
(`A = 6 - z1 - z2`, `C = 2/(z1 z2)`, `T = 3 - 1/(1-z1) - 1/(1-z2)`, `u = 1-T`,
`B·u = T(1-A-C) - 3 + 2A`); order 112 ⇔ that cubic splits completely over Q. Funnel: real test
`sign(disc) > 0` (division-free doubles, cache-tiled) → Legendre chain on `disc` over 60 primes
(from 101 up) → survivors printed for exact re-verification.

```
value list                 12 453 069 rationals z = n/d, max(|n|,d) ≤ 3200
cold table                 8.52 GB
FINAL (height ≤ 2250)      1.895 e13 pairs,  1.012 e13 passed the real test (53.3884 %),
                           ZERO candidates,  27 351 s
```

Both filters are **necessary conditions**, so the negative is exhaustive inside the box:

> **No order-112 configuration has two of its five `z`-values both of naive height ≤ 2250.**
> (Equivalently: at most one `z_i` has height ≤ 2250.)

**Why it is superseded**: the order-112 locus is a *curve* on each slice (genus 4 at `c = 2`,
genus 13 elsewhere), hence Faltings-finite. Enlarging a search box cannot decide a Faltings-finite
set; a rank + Chabauty computation can, and did (for `c = 2`). Search was the wrong instrument;
it is kept only as the quantitative statement above.

**Caveats on the sieve, honestly:** (i) the real-sign test is done in **double precision** with a
`1e-9` relative cancellation guard — it errs toward *keeping* borderline pairs, but a true solution
whose discriminant is positive yet catastrophically cancelling could in principle be dropped;
(ii) a pair is only reported if at least 40 of the 60 primes gave `W ≢ 0 (mod p)`;
(iii) `disc = 0` (repeated roots) is excluded by the strict `> 0` test — but such a point is
degenerate anyway.

---

## 6. Every number, in one place

| quantity | value | source |
|---|---|---|
| universal relations | `e1=6, e5=2, e3=2e4-2` | `xdecide.log` PART 0 (asserted) |
| `dim / deg S` | 2 / 20, irreducible, reduced | `surface.log` |
| `ω_S` | `O_S(4)` (Gorenstein) | `surface.log` |
| `dim Sing(S)` | 1, 25 components, not at infinity | `surface.log` |
| double-discriminant numerator | degree 58; rational roots **{1, 2}** only | `specialc.log` |
| slices scanned for low genus | 429 (`cscan`) and 464 (`slicegenus`) — only `c=2` special | those logs |
| `genus(Cab_c, 3-root, X_c)` | `(0,1,4)` at `c=2`; `(1,4,13)` at 18 other `c` | `surface.log` |
| `E1 = E2` | `y^2 = x^3 - 12x`, cond 288, `MW=[2,0]` proven | `xdecide.log` |
| `H` | genus 2, disc `2^54 3^46`, torsion `[2,6]`, RankBounds `0..1` | `xdecide.log` |
| generator of `Jac(H)(Q)` | `(x^2-2x, 8x, 2)`, `ĥ = 0.50118239204717836371111998436553` | `verify.log` |
| saturation | index 1 at every prime ≤ 50 | `verify.log` |
| `HeightConstant(J_H)` | 22.2342101128697832838282821298 | `verify.log` |
| `#H(Q)` | 8 (Chabauty = naive `Bound 200000`) | `chabauty.log`, `verify.log` |
| `#X(Q)` (affine) | 13, `T ∈ {-6,0,12/7,2,3}`, all degenerate | `chabauty.log`, `verify.log` |
| `Sing(Cab_2)` | one node `(2:1:1)`, two places of degree 1 | `gap.log` |
| `T=3` fibre of `X` | residue field `Q(√-1)` — irrational | `gap.log` |
| `Aut(H)/Q` | order 4; `GeometricAutomorphismGroup = V4` (bielliptic) | `verify.log` |
| `L`-poly of `H` | **reducible at all 44 good primes in [5,200]** ⇒ `Jac(H)` splits | `verify.log` |
| slices with `rank(Cab_c) = 0` | exactly **`c = 2/3`** (cond 48, `Z/8`) out of 71 completed | `cabrank.log`, `cabrank_c23.log` |
| local point counts | 594 `(c,p)` pairs, 31 slices, `p < 100`, **0** empty | `localobs.log` |
| C sieve | `1.895e13` pairs, height ≤ 2250, **0** candidates, 27 351 s | `splitsieve_H3200.log` |

Note: `Jac(H)` **splits** — `H` is bielliptic and `L_p` is reducible at all 44 good primes tested.
That is fine and expected (`H` is a quotient of the `(Z/2)^2`-cover `X`); it is *not* a candidate
curve for the project, it is a proof device. No simplicity certificate is needed or claimed for it.

---

## 7. What I did **not** check — the honest list

1. **Only two slices of infinitely many are decided** (`c = 2`, `c = 2/3`). For every other
   rational `c`, `X_c` has genus 13, `Faltings` gives finiteness and nothing more. The contact-7
   order-112 *surface* is **not** decided. Anyone quoting this result must say "the `c = 2` slice"
   (or "`z = 2` is not among the five `z`-values"), never "contact-7 cannot produce order 112".
2. The local search of §4b covers **31 slices and primes `p < 100` only**, and "genus preserved"
   is used as a *proxy* for good reduction — I did not build integral models or verify good
   reduction properly. So §4b is evidence that no small-prime obstruction exists, not a theorem
   that none exists. No `p`-adic (as opposed to `F_p`) solvability, no Sha/descent, no
   Brauer–Manin was attempted, on `X_c` or on `S`.
3. `Sing(S)` has dimension 1 and is **not** confined to the hyperplane at infinity, so the
   "general type / Bombieri–Lang" reading of §1 is heuristic. No resolution, no Kodaira dimension
   of a resolution, no actual `p_g` of `S` was computed.
4. The rank census of §4 completed **71 of 83** slices. Ranks came from `RankBounds` with
   `lo = hi` in every reported case. The 12 unfinished slices (`-61/450, 625/114, -41/124, -16/17,
   289/124, 297/133, 361/525, 50/33, 189/155, 121/155, 49/36, 50/63`) are recorded in the log with
   their resume command; I do **not** know their ranks, so I cannot exclude a third rank-0 slice
   among them.
5. **Only `[2,2,2,14]` is addressed.** The other order-112 groups (`[2,2,28]`, `[4,28]`, `[2,56]`,
   `Z/112`) require 2- or 4-divisibility of a class rather than full splitting of the quintic, and
   are **completely untouched** by this lane.
6. **Only the contact-7 chart is addressed.** Any other construction of a rational 7-torsion class,
   and every non-contact route to order 112, is untouched.
7. `Automorphisms(X)` on the genus-4 function field did **not terminate** (~5 min, killed).
   `code/claude_ov_c7z_aut.m` (a structural branch-data bound on `Aut(X)`) was **written but never
   run** — no log exists. It is not needed for the theorem.
8. No search was run over the *undecided* slices this session; the only evidence there is the C
   sieve of §5 (which is a statement about pairs of small `z`, not about any single slice).

---

## 8. Exact resume commands

Everything below runs from the repo root; **never call `magma` directly**.

```bash
# finish / redo the rank census of the slices  (~15-25 min, one slot)
code/claude_magma_slot.sh -b MemGB:=6 code/claude_ov_c7z_cabrank.m \
    > results/claude_ov_c7z_cabrank.log 2>&1 &

# a single slice, e.g. c = 5/7
code/claude_magma_slot.sh -b MemGB:=6 ONLYN:=5 ONLYD:=7 code/claude_ov_c7z_cabrank.m \
    > results/claude_ov_c7z_cabrank_c57.log 2>&1 &

# genus table for more slices: edit `clist` at the bottom of surface.m, then
code/claude_magma_slot.sh -b MemGB:=16 code/claude_ov_c7z_surface.m \
    > results/claude_ov_c7z_surface.log 2>&1 &

# re-verify the whole c=2 decision from scratch  (~2 min, ~4 min, ~6 min)
code/claude_magma_slot.sh -b MemGB:=12 BH:=1200 code/claude_ov_c7z_xdecide.m  > results/claude_ov_c7z_xdecide.log  2>&1
code/claude_magma_slot.sh -b            code/claude_ov_c7z_chabauty.m         > results/claude_ov_c7z_chabauty.log 2>&1
code/claude_magma_slot.sh -b            code/claude_ov_c7z_verify.m           > results/claude_ov_c7z_verify.log   2>&1
code/claude_magma_slot.sh -b            code/claude_ov_c7z_gap.m              > results/claude_ov_c7z_gap.log      2>&1

# local point counts on the slice fibres  (~25 min for 31 slices, one slot)
code/claude_magma_slot.sh -b MemGB:=6 code/claude_ov_c7z_localobs.m \
    > results/claude_ov_c7z_localobs.log 2>&1 &

# the (never-run) automorphism bound for X
code/claude_magma_slot.sh -b code/claude_ov_c7z_aut.m > results/claude_ov_c7z_aut.log 2>&1 &

# the superseded C sieve (do NOT restart it; kept for reproducibility)
gcc -O3 -march=native -fopenmp -o /tmp/splitsieve code/claude_ov_c7z_splitsieve.c -lm
/tmp/splitsieve 3200 192          # 8.5 GB, ~8 h to height 2250 on 192 cores
```

### What it would take to finish the surface

Ranked by plausibility, none of them cheap:

1. **Mordell–Weil sieve on `Cab_c`, slice by slice.** `Cab_c` is elliptic of rank 1–3 with known
   generators; the condition "`d1(P)` is a square" is a quadratic-residue condition on
   `P = n_1 G_1 + … + T`. Sieving `(n_i)` mod many primes gives strong per-slice negatives and,
   with a Chabauty-type height bound on the genus-4 cover, could in principle be made complete.
   This is the concrete next move and it is not currently implemented.
2. **Genus-4 Chabauty / covering collections** on the three-root fibre — no Magma support; would
   need a hand-built descent to elliptic-curve Chabauty over a number field.
3. **A uniform argument on `S`** (a second fibration, or a `(Z/2)`-cover trick that works in the
   parameter `c`). Nothing found; note that the symmetric relations `Σz = 6, Πz = 2, Σ1/(1-z) = 3`
   are invariant under permuting the `z_i`, so every configuration lies on **five** slices at once
   — deciding *any one* of its `z`-values kills it. That is why deciding small-height slices is
   worth more than it looks, and why the pair-sieve of §5 is the right complementary statement.
4. ~~Local obstruction hunt~~ — **done and negative** (§4b): 594 `(c,p)` counts, none empty.
   A deeper local attack would have to be `p`-adic/Sha-flavoured rather than a point count.
