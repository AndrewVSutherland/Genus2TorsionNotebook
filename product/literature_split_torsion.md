# Torsion groups known to occur on geometrically split genus-2 Jacobians over Q — literature compilation

Compiled 2026-08-12 (literature-search session). Scope: finite abelian groups G for which a
genus-2 curve C/Q is **known** with J = Jac(C) geometrically split (Qbar-isogenous to a product
of two elliptic curves) and G ⊆ J(Q)_tors (or G = J(Q)_tors where marked exact).

**RESOLUTION (2026-08-13) of the flagged P&P order-48 item**: the models are
printed in Platonov, Russian Math. Surveys 69:1 (2014) §6 (free full text via
mathnet.ru, DOI 10.1070/RM2014v069n01ABEH004877):
f36 = (x+2)(3x²−6x+8)(3x³+6x²+3x+4), f48,1 = (x−2)(x+2)(x⁴−10x²−3),
f48,2 = (x²+x+2)(x⁴+3x³+21x²−27x+18); Thm 6.2 asserts exact cyclic torsion
36/48/48, and §7 proves NON-simplicity of J36, J48,1, J48,2 over Qbar via
explicit zero-divisors in End⁰ (so no threat to Elkies' simple cyclic record
40).  Independent Magma verification (product/logs/verify_pp_curves.log):
exact torsion [36]/[48]/[48] confirmed, L-polys reducible at 44/44 good
p < 200, and the curves are precisely LMFDB 4860.f1 (= f36, factors implied by
trace match), 1764.a1 (= f48,1; factors cond 21 [8] × cond 84 [6]) and
5292.c2 (= f48,2; factors cond 126 [6] × cond 42 [8]) — i.e. the minimal DB
split witnesses for [36] and [48] are P&P's own 2012 curves.

Evidence tags used throughout:
- **[L]** stated in the cited literature (theorem/section quoted from the actual source, all sources
  re-verified this session from the papers themselves);
- **[M]** verified in this session with Magma V2.29-9: `TorsionSubgroup` (provable) on an integral
  simplified model; "D4 certificate" = Leprévost's absolute-simplicity criterion (JTNB 7 (1995),
  Lemme 3.1.2: some good prime p with irreducible degree-4 Frobenius charpoly with Galois group of
  order 8 ⇒ J absolutely simple); "trace matching" = the degree-4 L-polynomial factors into two
  elliptic factors at *every* good p < 200, and the factor traces match a specific pair of LMFDB
  elliptic isogeny classes complementarily at every p < 100 (overwhelming evidence of a Q-isogeny
  to that product, but not a formal isogeny proof);
- **[DB]** LMFDB (g2c database; `torsion_subgroup` is the provably computed rational torsion of the
  Jacobian, `is_simple_geom`/`geom_end_alg` are the rigorously certified endomorphism data of
  Costa–Mascot–Sijsling–Voight).

"Split/Q" means J is Q-isogenous to a product of elliptic curves over Q; "Qbar only" means J is
simple over Q but geometrically split. Groups are written by invariant factors [d1,d2,...], d1|d2|...

---

## (a) Merged table

### a.1 Groups of order ≥ 19 (the headline entries)

| Group | Order | Split | Exact? | Example curve (y² + h·y = f) | Source (theorem) | URL |
|---|---|---|---|---|---|---|
| [2,2,4,8] | 128 | Q | contains | family, positive-rank elliptic surface (glue E^t_{2,8}, E^u_{2,8}; eq. (5)) | HLP Thm 1 (Table 1) + §3.7 | arxiv.org/abs/math/9809210 |
| [2,2,24] | 96 | Q | contains | family on rank-2 curve Y²=X³−1681X (glue E^t_{2,6}, E^u_{2,8}) | HLP Thm 1 + §3.7 | arxiv.org/abs/math/9809210 |
| [2,6,6] | 72 | Q | contains | family, pos.-rank elliptic surface (glue E^t_{2,6}, E^u_{2,6}, system (6)) | HLP Thm 1 + §3.7 | arxiv.org/abs/math/9809210 |
| [6,12] | 72 | Q | contains | family, pos.-rank elliptic surface (Δ-matching N=N'=12) | HLP Thm 1 + §3.5 | arxiv.org/abs/math/9809210 |
| [70] | 70 | Q | **exact [M]** | y²+(2x³−3x²−41x+110)y = x³−51x²+425x+179 (3-gluing of 858k1 [7-tors] and 66c2 [10-tors]) | Howe 2015, Thm 2.4 (+ Rem. 2.5, 2.6) | arxiv.org/abs/1407.2654 |
| [8,8] | 64 | Q | contains | 2-param family (glue two curves with (2,8)-structure, special 2-tors pts NOT identified) | HLP Thm 1 + §3.4 | arxiv.org/abs/math/9809210 |
| [2,4,8] | 64 | Q | contains | 2-param family (same, special 2-torsion points identified) | HLP Thm 1 + §3.4 | arxiv.org/abs/math/9809210 |
| [63] | 63 | Q | **exact** [L]+[M] | y² = 897x⁶ − 197570x⁴ + 79136353x² − 146398496 (glue Δ₇(−16/3)≡Δ₉(4); exactness via #J(F₅)=63) | HLP eq. (4), §3.6 | arxiv.org/abs/math/9809210 |
| [60] | 60 | Q | contains | family on rank-1 curve 900A1 (Δ-matching N=10, N'=12, t=1/3) | HLP Thm 1 + §3.5, eq. (3) | arxiv.org/abs/math/9809210 |
| [5,10] | 50 | Q | contains | family, pos.-rank elliptic surface (Δ-matching N=N'=10) | HLP Thm 1 + §3.5 | arxiv.org/abs/math/9809210 |
| [7,7] | 49 | Q | **exact** [L]+[M] | y² = x⁶ + 3025x⁴ + 3232987x² + 869675859 (glue Δ₇(7)≡Δ₇(−14/13); exactness via F₅) | HLP §3.6 | arxiv.org/abs/math/9809210 |
| [48] | 48 | Q | **exact [M]** | y²+(x²+x)y = x⁶−3x⁵−5x⁴+14x³+8x²−16x, and y²+(x²+x)y = x⁶−x⁵+5x⁴−11x³+10x²−6x+2; infinitely many, parameterized by rank-2 curve D: y²=x³+14x²+196x | Howe 2015, Thm 3.1–3.2, Cor 3.5, Rem 3.3 | arxiv.org/abs/1407.2654 |
| [48] | 48 | Q | **exact [M]** | two curves predating Howe: f48,1 = (x−2)(x+2)(x⁴−10x²−3) = LMFDB 1764.a1 (factors: cond-21 [8] × cond-84 [6]) and f48,2 = (x²+x+2)(x⁴+3x³+21x²−27x+18) = 5292.c2; models from Platonov RMS 69:1 (2014) §6, non-simplicity proven there §7, Magma re-verified — see the resolution block above | P&P 2012b Thm 3; Platonov RMS 69:1 §6–7 | doi.org/10.1134/S1064562412050304 |
| [2,24] | 48 | Q | contains | 2-param family (glue (2,6)-structure with (2,8)-structure) | HLP Thm 1 + §3.4 | arxiv.org/abs/math/9809210 |
| [45] | 45 | Q | contains | infinitely many; glue E₉^(−5) with 5-torsion curves; rank-2 aux. curve y²=Δ₉(−5)t(t²−11t−1), cond. 13838400 | HLP Thm 1 + §3.6 | arxiv.org/abs/math/9809210 |
| [40] | 40 | Q | contains | family, pos.-rank elliptic surface (Δ-matching N=8, N'=10) | HLP Thm 1 + §3.5 | arxiv.org/abs/math/9809210 |
| [3,12] | 36 | Q | contains | 2-param family (Δ-matching N=12, N'=6) | HLP Thm 1 + §3.5 | arxiv.org/abs/math/9809210 |
| [6,6] | 36 | Q | contains | 2-param family (glue two (2,6)-structures) | HLP Thm 1 + §3.4 | arxiv.org/abs/math/9809210 |
| [6,6] | 36 | Q | **exact [M]+[DB]** | X₀(28): y²+(x²+x)y = x⁶+3x⁵+6x⁴+7x³+6x²+3x+1 = LMFDB 196.a.21952.1; J ~ 14.a × 14.a | Magma model + LMFDB | lmfdb.org/Genus2Curve/Q/196/a/21952/1 |
| [36] | 36 | Q (trace matching [M]) | **exact [M]** | y²+(6x³−3x²−x+2)y = 3x³−4x²+2x (J ~ 90.c × 270.d); y²+(6x³+3x²−x+2)y = 2x²+2x (J ~ 54.b × 90.c, i.e. Z/9-curve × Z/12-curve) | P&P Dokl. Math. 86 (2012) Thm 2; Platonov–Zhgun–Petrunin Dokl. Math. 87 (2013) Thm 3 (models as printed in Howe 2015 Table 2); splitness is THIS session's observation, not claimed in the sources | doi.org/10.1134/S1064562412050304 |
| [35] | 35 | Q | contains | infinitely many; glue E₇^(−1) with 5-torsion curves; rank-2 aux. curve y²=−26t(t²−11t−1), cond. 54080 | HLP Thm 1 + §3.6 | arxiv.org/abs/math/9809210 |
| [2,2,8] | 32 | Q | **exact [DB]** | y²+(x³+x)y = −3x⁴+7x²−5 (LMFDB 360.a.6480.1, and one more curve) | LMFDB | lmfdb.org/Genus2Curve/Q/360/a/6480/1 |
| [30] | 30 | Q | contains | 2-param family (Δ-matching N=10, N'=6) | HLP Thm 1 + §3.5 | arxiv.org/abs/math/9809210 |
| [3,9] | 27 | Q | contains | 2-param family (Lemma 11 with E=E₉^t) | HLP Thm 1 + §3.6 | arxiv.org/abs/math/9809210 |
| [3,9] | 27 | Q | **exact [DB]** | y²+(x²+x+1)y = x⁶−3x⁵+5x⁴−6x³+x (LMFDB 2916.a.139968.1, cond. 54², Q×Q) | LMFDB | lmfdb.org/Genus2Curve/Q/2916/a/139968/1 |
| [27] | 27 | Q | **exact** [L]+[M] | C₂₇: y² = (2x³−15x²−3)(2x³−15x²+12x−3); explicit degree-2 cover X=(x−1)²/(2x³−15x²−3) to E₁: Y²=648X³+297X²+18X+1 (3-torsion), so Jac ~_Q E×E' with points of order 3 and 9 | Leprévost, JTNB 7 (1995), Thm 1.2.1 + §3.2 ("Les courbes C₂₁, C₂₄,₂ et C₂₇") | numdam.org/item/JTNB_1995__7_1_283_0/ |
| [27] | 27 | Q (trace matching [M]) | **exact [M]** | y²+(2x³+3x²−3x+2)y = 6x³+6 (Howe's "new" order-27 curve; J ~ 54.b × 162.d — splitness observed this session, not stated by Howe) | Howe 2015, §4 Table 2 ("New") | arxiv.org/abs/1407.2654 |
| [5,5] | 25 | Q | **exact [M]** | X₀(22): y² = −2x⁶−10x⁴+26x²+242 (conductor 121; J ~ X₁(11)×X₀(11) = 11.a×11.a, (2,2)-isogenous) | HLP Cor. 5 + Remark (model); torsion computed this session; consistent with cuspidal C(22)≅(Z/5)² and Ohta/Yoo (see §(b6)) | arxiv.org/abs/math/9809210 |
| [24] | 24 | Q | **exact** [L]+[M] | C₂₄,₂: y² = (x²−x+1)(x⁴−3x³+8x²−3x+1); bielliptic via (x,y)↦(1/x, y/x³); covers E₁: Y²=(X+3)(X²−X+4) (8-torsion) and E₂: Y²=(1+3X)(1−X+4X²) (3-torsion) | Leprévost, JTNB 7 (1995), Thm 1.2.1 + §3.2 | numdam.org/item/JTNB_1995__7_1_283_0/ |
| [24] | 24 | Q | exact [DB] | y²+(x³+1)y = x⁵+3x⁴+3x³+3x²+x (LMFDB 450.a.2700.1; 3 such curves in LMFDB) | LMFDB | lmfdb.org/Genus2Curve/Q/450/a/2700/1 |
| [2,12] | 24 | Q | exact [DB] | y²+(x³+1)y = x⁵−4x⁴−9x³+28x²−6x−16 (LMFDB 450.a.36450.1) | LMFDB | lmfdb.org/Genus2Curve/Q/450/a/36450/1 |
| [2,2,6] | 24 | Q | exact [DB] | y²+xy = 10x⁵−18x⁴+8x³+x²−x (LMFDB 600.a.18000.1) | LMFDB | lmfdb.org/Genus2Curve/Q/600/a/18000/1 |
| [21] | 21 | Q | contains | 2-param family (Lemma 11 applied to universal 7-torsion curve) | HLP Thm 1 + §3.6 | arxiv.org/abs/math/9809210 |
| [21] | 21 | Q | **exact** [L]+[M] | C₂₁: y² = 4x⁶−12x⁵+13x⁴−6x³+3x²−2x+1; bielliptic: covers E₁: Y²=16X³−8X²+9X+2 (3-torsion) and E₂: Y²=2X³+9X²−8X+16 (7-torsion). Same G2-invariants as LMFDB 676.a.5408.1, i.e. Jac(C₂₁) ~ J₀(26) (isogeny class 676.a = 26.a×26.b) [M] | Leprévost, JTNB 7 (1995), Thm 1.2.1 + §3.2 | numdam.org/item/JTNB_1995__7_1_283_0/ |
| [21] | 21 | Q | **exact [M]+[DB]** | X₀(26): y²+(x³+1)y = 2x⁵+2x⁴+4x³+2x²+2x = LMFDB 676.a.562432.1; J ~ 26.a (Z/3) × 26.b (Z/7) | Magma model + LMFDB | lmfdb.org/Genus2Curve/Q/676/a/562432/1 |
| [21] | 21 | **Qbar only** | **exact** [L]+[M]+[DB] | X₁(18): y² = x⁶−4x⁵+10x⁴−10x³+5x²−2x+1 = LMFDB 324.a.648.1 (Q-simple, CM by imag. quad. field over Q, geometrically M₂(Q)) | torsion: Ogg 1973 (PSPM 24, p. 226), quoted in Leprévost JTNB §1.1; endomorphisms: LMFDB | lmfdb.org/Genus2Curve/Q/324/a/648/1 |
| [2,10] | 20 | **Qbar only** | **exact** [L]+[M]+[DB] | X₁(16): y² = (x−1)(x+1)(x²+1)(x²+2x−1) = LMFDB 256.a.512.1 (Q-simple, geometrically split) | Ogg 1973 p. 226 via Leprévost JTNB §1.1; LMFDB | lmfdb.org/Genus2Curve/Q/256/a/512/1 |
| [20] | 20 | Q | contains | 2-param family (Δ-matching N=10, N'=4) | HLP Thm 1 + §3.5 | arxiv.org/abs/math/9809210 |
| [20] | 20 | Q | exact [DB] | y²+(x²+x)y = 3x⁵−13x³+16x²+65x+40 (LMFDB 2250.a.324000.1) | LMFDB | lmfdb.org/Genus2Curve/Q/2250/a/324000/1 |
| [19] | 19 | **Qbar only** | **exact** [L]+[M]+[DB] | X₁(13): y² = x⁶−2x⁵+x⁴−2x³+6x²−4x+1 = LMFDB 169.a.169.1 (Q-simple; CM field Q(√−3) acts /Q; geometrically M₂(Q)) | order-19 class: Leprévost JTNB §1.1 (explicit divisor); J₁(13)(Q)≅Z/19: Ogg 1973 p. 226 (also Mazur–Tate 1973) | lmfdb.org/Genus2Curve/Q/169/a/169/1 |

### a.2 Groups of order ≤ 18 (all split/Q unless noted; witnesses with exact torsion)

All of these also occur inside the HLP families above ("contains"); the point of this sub-table is
an *exact-torsion* witness for each group on a geometrically split Jacobian. All rows [DB] (LMFDB
`torsion_subgroup` + `is_simple_geom=false`), with count of such curves in LMFDB and one example.

| Group | Order | LMFDB example | Notes |
|---|---|---|---|
| [3,6] | 18 | 2700.a.81000.1 (6 curves) | y²+(x³+1)y = 5x⁵+26x⁴+12x³+26x²+5x |
| [2,8] | 16 | 1008.a.27216.1 (11 curves) | |
| [4,4] | 16 | 504.a.27216.1 (3 curves) | y²+(x³+x)y = 3x⁴+15x²+21; group not in HLP's Table 1 |
| [2,2,4] | 16 | 4950.a.742500.1 (2 curves) | |
| [15] | 15 | 2500.a.50000.1 = **X₀(50)** (7 curves) | y²+(x³+1)y = x⁵+2x³+x; J ~ 50.a (Z/3) × 50.b (Z/5) [M] |
| [12] | 12 | 1088.b.2176.2 (41 curves) | |
| [2,6] | 12 | 1170.a.10530.1 (21 curves) | |
| [10] | 10 | 363.a.43923.1 (5 curves) | |
| [9] | 9 | 1696.b.434176.1 | **Qbar only** (is_simple_base = true, End_Q = Z, geometrically Q×Q) |
| [3,3] | 9 | 676.b.17576.1 [M], 26244.c.157464.1 (8 curves) | 676.b: y²+(x²+x)y = −(x⁶−3x⁵+6x⁴−6x³+6x²−3x+1), M₂(Q) over Q |
| [8] | 8 | 101400.a.608400.1 (27 curves) | |
| [2,4] | 8 | 10080.c.141120.1 (93 curves) | |
| [2,2,2] | 8 | 12800.c.128000.1 (15 curves) | |
| [7] | 7 | 10816.g.86528.1 (11 curves) | |
| [6] | 6 | (217 curves) | |
| [5] | 5 | 10000.a.160000.1 (75 curves); also 2500.a.400000.1 (isogenous to J₀(50)) [M] | |
| [4], [2,2], [3], [2], [1] | ≤4 | (339 / 268 / 289 / 944 / 526 curves) | |

Also exact [M] this session: **X₀(37)** = LMFDB 1369.a.50653.1 (Magma model y²−x³y =
2x⁵−5x⁴+7x³−6x²+3x−1; G2-invariants match), J₀(37)(Q)_tors = **[3]**, J ~ 37.a × 37.b (split/Q);
the second curve 1369.a.1369.1 in the class also has torsion [3].

### a.3 Records and structural remarks

- Largest **point order** known on a split genus-2 Jacobian /Q: **70** (Howe 2015; split/Q).
  Largest **group order**: **128** = [2,2,4,8] (HLP family). On the absolutely simple side the
  record point order is **40** (Elkies 2001–02; see §(d)), with P&P's order-48 curves of unknown
  split/simple status (§(c)).
- For J Q-isogenous to E₁×E₂ *via the (2,2)- or (3,3)-gluings used in all the above
  constructions*, odd-order rational points push to rational points on E₁×E₂ (up to the small
  kernel), so the possible odd parts are governed by Mazur's theorem on the factors; e.g. a
  rational point of prime order ℓ ∈ {11,13,17,19,23,29,31} forces such a J to be **not**
  Q-split (for general Q-isogenies a caveat remains: quotients along ℓ-congruences with
  inverse isogeny characters could in principle create rational ℓ-torsion, so this is a
  heuristic constraint outside the gluing framework, not a theorem cited from the literature).
  The known Qbar-only-split examples with orders 19, 21 (X₁(13), X₁(18)) show the geometric
  split condition alone does not obstruct ℓ > 12.
- Gaps: no geometrically split example appears to be known (literature or LMFDB) with exact
  cyclic torsion [14], [16], [2,14], [2,16], [4,8], [4,12], [2,20] — several of these are
  plausible targets for (2,2)/(3,3)-gluing constructions.

---

## (b) Per-paper notes

### (b1) Howe, Leprévost, Poonen, "Large torsion subgroups of split Jacobians of curves of genus two or three", Forum Math. 12 (2000) 315–364; arXiv:math/9809210

The central reference. Theorem 1 (genus 2): for every G in their Table 1 there is a family of
genus-2 curves /Q, parameterized by the rational points of a nonempty Zariski-open subset of the
listed variety, with G ⊆ J(Q)_tors ("contained": exactness is proven only for the two single
curves, orders 63 and 49). All constructions produce J **(2,2)-isogenous over Q** to E₁×E₂
(Prop. 3 & 4: gluing along a Galois-module isomorphism ψ: E[2]→F[2] not induced by an
isomorphism of curves), so every entry is split over Q. Complete genus-2 table
(group, |G|, parameterizing variety):

- Z/20 (20, P²); Z/21 (21, P²); Z/3×Z/9 (27, P²); Z/30 (30, P²); Z/35 (35, pos.rank ell. curve);
- Z/6×Z/6 (36, P²); Z/3×Z/12 (36, P²); Z/40 (40, pos.rank ell. surface); Z/45 (45, pos.rank ell. curve);
- Z/2×Z/24 (48, P²); Z/7×Z/7 (49, P⁰ = single curve); Z/5×Z/10 (50, pos.rank ell. surface);
- Z/60 (60, pos.rank ell. curve); Z/63 (63, P⁰); Z/8×Z/8 (64, P²); Z/2×Z/4×Z/8 (64, P²);
- Z/6×Z/12 (72, pos.rank ell. surface); Z/2×Z/6×Z/6 (72, pos.rank ell. surface);
- Z/2×Z/2×Z/24 (96, pos.rank ell. curve); Z/2×Z/2×Z/4×Z/8 (128, pos.rank ell. surface).

Construction sections: §3.4 (both curves with full rational 2-torsion; (2,2N)-structures, N≤4)
gives [8,8]/[2,4,8] (N₁=N₂=4), [2,24] (N₁=3,N₂=4), [6,6] (N₁=N₂=3). §3.5 (one rational
2-torsion point; solve Δ_N(t)y² = Δ_N'(u)) gives [20], [30], [3,12], [40], [5,10], [6,12], [60].
§3.6 (trivial rational 2-torsion, odd N≤9; matching cubic 2-division fields) gives [21], [3,9],
[63], [7,7], [35], [45]. §3.7 ("gaining 2-power torsion" via a Galois-cohomology criterion,
Prop. 12) gives [2,2,4,8], [2,6,6], [2,2,24]. Genus-3 analogue in Theorem 2/Table 2 (largest
group order 864, plane quartic 15625(X⁴+Y⁴+Z⁴)−96914(X²Y²+X²Z²+Y²Z²)=0, exact torsion
Z/6×Z/12×Z/12 — outside our scope).

Also: Corollary 5 identifies y² = −2x⁶−10x⁴+26x²+242 as X₀(22), Jacobian of conductor 121
((2,2)-isogenous to X₁(11)×X₀(11)); Corollary 8 (rank-28 genus-2 Jacobian, Q-isogenous to E×E).
Announcement note: Howe–Leprévost–Poonen, "Sous-groupes de torsion d'ordres élevés de
Jacobiennes décomposables de courbes de genre 2", C. R. Acad. Sci. Paris 323 (1996) 1031–1034
(table of the same groups up to order 128, with the Z/63 example).

### (b2) Howe, "Genus-2 Jacobians with torsion points of large order", Bull. LMS 47 (2015) 127–135; arXiv:1407.2654

- Thm 2.1: five curves with a rational point of order **28**, produced by the
  Bröker–Howe–Lauter–Stevenhagen (3,3)-gluing algorithm applied to pairs (4-torsion curve,
  7-torsion curve) — all split over Q by construction; e.g. C₂₈,₁: y²+(x²+x)y =
  x⁶+3x⁵+5x⁴−4x²−10x+4 (Remark 2.2: exact torsion Z/28 for all five). [M: exact [28] confirmed.]
- Thm 2.4: single curve C₇₀ (order **70**, from 3-gluing of 858k1 and 66c2), the largest order
  known; Remark 2.6: exact torsion Z/70. [M: confirmed.]
- §3 (Thm 3.1, Thm 3.2, Rem 3.3, Cor 3.5): infinite family with a rational point of order
  **48** (2-gluing of a 6-torsion curve with a fixed 8-torsion curve, then doubling a 24-point
  via HLP Prop. 12); parameterized by the rank-2 elliptic curve y²=x³+14x²+196x. Prior to
  this, only P&P's two order-48 curves existed [P&P Dokl. Math. 86 (2012), Thm 3].
  [M: both printed curves have exact torsion [48].]
- §4 Table 2: small-coefficient curves found by search, with attributions: order 27 —
  Leprévost JTNB Thm 1.2.1 curve + three new; 28 — the two P&P curves; 29 — Leprévost;
  33 — P&P; 34 — Elkies; 36 — the two P&P curves; 39 — Elkies + one new. Howe makes **no
  split/simple claims** for the §4 curves. This session's verdicts [M]: new-27 #1
  (y²+(x³−2x+1)y=x³): absolutely simple (D4 @ p=5); new-27 #2 (y²+(2x³+3x²−3x+2)y=6x³+6):
  split/Q ~ 54.b×162.d; new-27 #3 (y²+(6x³+9x²+6x−1)y=−3x²): Q-simple (irreducible Frobenius
  charpolys at several p) but with heavily mixed reducibility — likely split over a quadratic
  field, status open; new-39 (y²+(6x³+6x²−7x−9)y=−2x−2): absolutely simple (D4 @ p=5).
- Intro survey (as of 2014): orders realized 1≤n≤30, 32≤n≤36, n ∈ {39,40,45,48,60,63};
  infinite families for 1≤n≤26 and n ∈ {30,32,35,40,45,60}; for n ∈ {14,16,18,22,26} the
  infinite families are from an *unpublished* Leprévost preprint (specializing his 7-, 9-, 11-,
  13-families to gain a rational 2-torsion point; 16 via his Manuscripta 75 methods).
- Also reports a targeted search for order **31** (with R. Bröker) over all y²=f, f quintic or
  sextic with |coeffs| ≤ 20, and y²+gy=f with small coefficients: **no example found**.

### (b3) Leprévost's corpus

- C. R. Acad. Sci. 313 (1991) 451–454 ("…d'ordre 13"): 1-parameter family with a rational
  divisor class of order 13. (gallica.bnf.fr/ark:/12148/bpt6k57325582/f455.image)
- C. R. Acad. Sci. 313 (1991) 771–774 ("…d'ordre 15, 17, 19 ou 21"): 1-parameter families,
  orders 15, 17, 19, 21. (…/f775.image)
- Manuscripta Math. 75 (1992) 303–326: genus-g families with points of order 2g²+2g+1 and
  2g²+3g+1 (g=2: 13, 15) over Q(t).
- Sém. Th. Nombres Paris 1991–92 (Birkhäuser 1993): genus-g family, order 2g²+4g+1 (g=2: 17).
- C. R. Acad. Sci. 316 (1993) 819–821: announcement of the JTNB results below.
- **JTNB 7 (1995) 283–306** ("Jacobiennes de certaines courbes de genre 2 : torsion et
  simplicité") — fully read this session (numdam):
  * §1.1: models of the genus-2 modular curves and their torsion (citing Ogg, PSPM 24 (1973),
    p. 226): X₁(13): y²=x⁶−2x⁵+x⁴−2x³+6x²−4x+1 with explicit order-19 divisor class
    D₀=(0,1)−(∞⁺); J₁(13)(Q)_tors ≅ Z/19; X₁(18): y²=x⁶−4x⁵+10x⁴−10x³+5x²−2x+1,
    J₁(18)(Q)_tors ≅ Z/21; X₁(16): y²=(x−1)(x+1)(x²+1)(x²+2x−1), J₁(16)(Q)_tors ≅ Z/2×Z/10.
  * Thm 1.2.1: single curves C_l, y²=f_l, whose Jacobians have a rational point of order
    l = 21, 22, 23, 25, 26, 27, 29, and two non-Qbar-isomorphic curves for l = 24:
    f₂₁ = 4x⁶−12x⁵+13x⁴−6x³+3x²−2x+1; f₂₂ = (2x²−2x+1)(2x⁴−2x³+x²−4x+4);
    f₂₃ = x⁶−10x⁵+33x⁴−36x³+28x²−16x+4; f₂₄,₁ = (2x²−2x−1)(2x⁴−10x³+7x²+4x−4);
    f₂₄,₂ = (x²−x+1)(x⁴−3x³+8x²−3x+1); f₂₅ = 36x⁶−156x⁵+241x⁴−192x³+102x²−36x+9;
    f₂₆ = (6x²−6x+1)(6x⁴−30x³+49x²−20x+4); f₂₇ = (2x³−15x²−3)(2x³−15x²+12x−3);
    f₂₉ = (2x−1)(2x⁵−x⁴−4x²+8x−4). [M: all torsion orders confirmed; exact torsion equals
    the cyclic group in every case.]
  * Thm 1.2.2 / Cor 1.2.3: 1-parameter family C_t with Z/3×Z/9 ⊆ J(Q): y² = x⁶+6(t−1)x⁵
    +3(t−1)(3t−5)x⁴−(18t²−34t+18)x³+3(t−1)(5t−3)x²−6t(t−1)x+t².
  * §3.2 split/simple verdicts (his own): **C₂₁, C₂₄,₂, C₂₇ are NOT simple over Q** (explicit
    bielliptic covers; see table rows above). **C₂₂, C₂₃, C₂₄,₁, C₂₆, C₂₉ are absolutely
    simple** (Lemme 3.1.2, the D4 criterion — the same criterion later used by Elkies and by
    this repo). **C₂₅ is simple over Q** (V₄ Galois group argument at p=5) and — per his
    numerics — "seems isogenous over Q(√3) to a product of elliptic curves, but we have no
    proof". The family C_t is generically simple over Q (specialization t=3 proven simple/Q);
    absolute simplicity left open. [M: D4 certificates reproduce all simplicity claims;
    reducibility patterns support the Q(√3) suspicion for C₂₅ and suggest the same phenomenon
    for C_t(3).]
- J. reine angew. Math. 473 (1996) 59–68 ("Sur une conjecture sur les points de torsion
  rationnels des jacobiennes de courbes"): review inaccessible this session (see §(c)).
- Manuscripta Math. 92 (1997) 47–63: genus-g families with a rational point of order 2g(2g+1)
  (g=2: **20**) and, for l = 2g²+5g+5 (g=2: **23**), families with a point of order l, l/2 or
  l/4 (zbMATH review). Split/simple status not addressed in the review.
- Unpublished: "A provisional report on the rational torsion groups of Jacobians of curves of
  genus two" (1996 preprint, cited by HLP as [29]) — per Elkies and Howe it contains a
  1-parameter family with a point of order **30** and the infinite families for 14, 16, 18,
  22, 26.

### (b4) Platonov–Petrunin(–Zhgun) (torsion via S-units / polynomial continued fractions)

- P&P, Dokl. Math. 85 (2012) 286–288 ("New orders of torsion points…"; Russian original Dokl.
  Akad. Nauk 443:6 664–667): contains (Thm 4, p. 288) a curve with a rational point of order
  **28** — the first ever; small model (Howe 2015 Table 2): y²+(3x³+2x²+1)y = −x²−x.
  [M: exact torsion [28]; **absolutely simple** (D4 @ p=7).]
- P&P, Dokl. Math. 86 (2012) 642–643 ("On the torsion problem…"; Russian original DAN 446:3
  263–264): Cor. 1: order **33** curve, small model y²+(3x³+9x²+x+2)y = −8x [M: exact [33];
  **absolutely simple** (D4 @ p=13)]; Thm 2: order **36** curve, small model
  y²+(6x³−3x²−x+2)y = 3x³−4x²+2x [M: exact [36]; **split/Q** ~ 90.c×270.d by trace matching];
  Thm 3 (p. 643): **two curves with points of order 48** (models recovered from the RMS 69:1
  (2014) survey §6 and Magma-verified, split per its §7 — see the resolution block above).
- Platonov–Zhgun–Petrunin, Dokl. Math. 87 (2013) 318–321 ("On the simplicity of Jacobians for
  hyperelliptic curves of genus 2 … with torsion points of high order"; DAN 450:4 385–388):
  proves simplicity results for their high-order curves and exhibits (Thm 3, p. 320) an
  order-28 curve, small model y²+(2x³−3x²+3x+4)y = 4x, and an order-36 curve, small model
  y²+(6x³+3x²−x+2)y = 2x²+2x. [M: the order-28 small model has **the same G2-invariants and
  |disc|** as Howe's glued C₂₈,₁ — i.e. it is the same curve up to isomorphism, and it is
  split/Q; the order-36 curve is split/Q ~ 54.b×90.c (a Z/9-curve times a Z/12-curve — in
  hindsight a (3,3)-gluing with 9·12/3 = 36). The simplicity assertions of the Russian school
  therefore pertain to their order-28 (2012a) and order-33 curves, which ARE absolutely
  simple; their order-36 (and possibly order-48) curves are split.]
- Summary of the school's contribution (per the Russian survey literature, e.g. Platonov's
  Uspekhi/RMS 69:1 (2014) survey and the swsys-web.ru account): before them the realized
  orders were 11,13,15,17,19,20,21,22,23,24,25,26,27,29,30,32,34,35,39,40,45,60,63
  (Flynn, Leprévost, HLP, Ogawa, Elkies); they added **28, 33, 36, 48** and completed
  "Poonen's conjecture" that every order ≤ 30 occurs.

### (b5) Kulesz (Dem'janenko–Manin families)

Kulesz, "Application de la méthode de Dem'janenko–Manin à certaines familles de courbes de
genre 2 et 3", J. Number Theory 76 (1999) 130–146 (doi:10.1006/jnth.1998.2339). Genus-2 family
a·y² = 6(x²−4)(x²+2)(x²+8), a ∈ N such that the CM curve y² = x³−ax has rank 1; the curves
admit two independent morphisms to that elliptic curve, so the Jacobians are **isogenous over Q
to a product of two elliptic curves** (split/Q) — but the paper's aim is the determination of
all rational points, and the only torsion visible on the family is 2-torsion
(J(Q)_tors ⊇ Z/2×Z/2 from the factored sextic). No new torsion groups for our table. His other
papers (Acta Arith. 87 (1999) 103–120; JLMS 2003 with genus-2/3 curves of high rank; C. R. 329
(1999) 503–506) are rank-oriented, again with bielliptic (split) families and only small
2-torsion stated.

### (b6) Modular curves of genus 2

Genus-2 X₀(N): N ∈ {22, 23, 26, 28, 29, 31, 37, 50}. Split cases (all split over Q, since the
weight-2 eigenforms of those levels are rational):

| Curve | Model (Magma SmallModularCurve / LMFDB) | J ~ | J(Q)_tors | Evidence |
|---|---|---|---|---|
| X₀(22) | y²−x³y = −x⁴+5x³−10x²+12x−8 (≅ HLP's y²=−2x⁶−10x⁴+26x²+242, cond. 121) | 11.a × 11.a | **[5,5]** | [M] both models; consistent with C(22)≅(Z/5)² and Ohta's theorem (odd part, squarefree N) |
| X₀(26) | LMFDB 676.a.562432.1 | 26.a × 26.b | **[21]** | [M]+[DB] |
| X₀(28) | LMFDB 196.a.21952.1 (cond. 196) | 14.a × 14.a | **[6,6]** | [M]+[DB] |
| X₀(37) | LMFDB 1369.a.50653.1 | 37.a × 37.b | **[3]** | [M]+[DB]; classical (cuspidal group of J₀(37) has order 3, Mazur–Swinnerton-Dyer) |
| X₀(50) | LMFDB 2500.a.50000.1 | 50.a × 50.b | **[15]** | [M]+[DB] |

Simple cases (RM, geometrically simple; excluded from the main table): X₀(23): torsion [11],
X₀(29): [7], X₀(31): [5] (the classical cuspidal groups, Ogg) [all M; D4 certificates confirm
geometric simplicity for all three].

Genus-2 X₁(N): N ∈ {13, 16, 18} — all three are **Q-simple but geometrically split**
(LMFDB: end_alg = CM (imaginary quadratic acting over Q), geom_end_alg = M₂(Q), is_simple_geom
= false; they are Qbar-isogenous to E² for a Q-curve E without CM):
J₁(13)(Q) ≅ Z/19 (Ogg 1973 p. 226; the order-19 class is the cuspidal divisor, cf. also
Mazur–Tate, "Points of order 13 on elliptic curves", Invent. Math. 22 (1973)); J₁(16)(Q)_tors ≅
Z/2×Z/10; J₁(18)(Q)_tors ≅ Z/21 (both Ogg 1973 p. 226). [All three confirmed exactly by M and
DB.] These are the only three groups currently known on split Jacobians that are NOT split
over Q — see also LMFDB 1696.b.434176.1 ([9], Qbar only).

Note: the further genus-2 curves X₀(N)/⟨w_d⟩ (Hasegawa, Proc. Japan Acad. 71 (1995) 235–239,
lists all genus-2 quotients) include many with split Jacobians; their rational torsion is less
documented and was not systematically compiled here.

### (b7) Elkies (web page, 2001–2004): absolutely simple side

"Curves of genus 2 over Q whose Jacobians are absolutely simple abelian surfaces with torsion
points of high order" (people.math.harvard.edu/~elkies/g2_tors.html). States the split state of
the art is HLP with largest order 63; on the simple side reports (all proven absolutely simple
by Leprévost's D4 criterion): a 1-parameter family with order **32** (sample curve
y²=(15x−1)(1056x⁴+156183x³+26297x²+649x−121)); two curves with order **34**
(y²=(9x²+2x+1)(32x³+81x²−6x+1) and y²=(10−x)(3x+2)(72x⁴+96x³+45x²−38x+5)); one with **39**
(y²=x⁶+4x⁴+10x³+4x²−4x+1); one with **40** (y²=(3x+4)(x⁴+5x³+8x²+(19/4)x+1)) — the simple
record. Prior simple records: 29 (Leprévost JTNB, single curve), 30 (Leprévost unpublished
1996, 1-param family). [M: exact torsion [32],[34],[34],[39],[40] and D4 certificates
all confirmed.]

### (b8) LMFDB (The LMFDB Collaboration, lmfdb.org, genus-2 curves over Q)

Census of exact torsion groups on curves with `is_simple_geom = false` (i.e. geometrically
split), from the g2c database (66,158 curves with |disc| ≤ 10⁶, of which 2,926 are
geometrically split): the 31 groups listed in §a.1/a.2. Largest torsion order among geometrically split curves in the database: 36
(X₀(28)). The database `torsion_subgroup` values are provably computed; endomorphism
certifications are rigorous (Costa–Mascot–Sijsling–Voight).

### (b9) Related but yielding no table entries

- Frengley, "On 12-congruences of elliptic curves" (arXiv:2208.05842, Math. Z. 2023?) and his
  thesis "Explicit moduli spaces for curves of genus 1 and 2" (Cambridge, 2023): produces
  infinitely many pairs of 12-congruent non-isogenous elliptic curves — exactly the input
  needed for (12,12)-gluings — but no genus-2 torsion application is stated.
- Bröker–Howe–Lauter–Stevenhagen, "Genus-2 curves and Jacobians with a given number of points"
  (LMS J. Comput. Math. 18 (2015)); their Algorithm 5.4 (explicit (3,3)-gluing) is the engine
  of Howe 2015; no torsion tables of its own.
- Fisher (7- and 11-congruences), Frengley–Fisher: machinery for (N,N)-split Jacobians with
  larger N; no rational-torsion applications located.
- Flynn, "Large rational torsion on abelian varieties" (J. Number Theory 36 (1990) 257–265) and
  "Sequences of rational torsions on abelian varieties" (Invent. Math. 106 (1991) 433–442):
  1-parameter genus-2 families with rational points of order 11 (and further orders growing
  with g); simplicity/splitness not addressed (generic members expected simple) — side list.

---

## (c) Claims NOT verified from the primary source (flagged)

1. **P&P order-48 curves** — **RESOLVED 2026-08-13**, see the resolution block at the top of
   this file: the models were recovered from Platonov, Russian Math. Surveys 69:1 (2014) §6
   (which supersedes the Doklady access problem), the survey's §7 proves all of J36, J48,1,
   J48,2 non-simple over Qbar, and independent Magma verification confirms exact torsion and
   identifies the curves as LMFDB 4860.f1, 1764.a1, 5292.c2.  (Original flag, retained for
   context: the Doklady paper itself was inaccessible and the equations were then unverified.)
2. **Ogawa, order 23** (Proc. Japan Acad. 70 (1994) 295–298, doi:10.3792/pjaa.70.295): equation
   not obtained (Project Euclid download blocked); listed on the strength of Howe 2015 [18] and
   the Russian survey. Note 23 is prime > 12, so (at least within the gluing framework) the
   curve cannot be split over Q; its geometric status is unknown to me.
3. **Leprévost, JRAM 473 (1996)** ("Sur une conjecture sur les points de torsion rationnels des
   jacobiennes de courbes"): content not accessed (zbMATH review license-blocked); it is cited
   by HLP among the genus-2 torsion family constructions. No entries were taken from it.
4. **Leprévost's unpublished 1996 preprint** (order-30 family; infinite families for 14, 16,
   18, 22, 26): reported by Elkies (g2_tors page) and Howe 2015; the preprint itself was not
   seen. No split claims are attached to it.
5. **C₂₅ split over Q(√3)**: Leprévost's own numerics (JTNB p. 303) + this session's mixed
   Frobenius-reducibility pattern (reducible exactly at a set of primes avoiding p ≢ ±1 mod 12)
   support it, but it is unproven; C₂₅ is proven simple over Q only.
6. **Splitness of the P&P order-36 curves and of Howe's second new order-27 curve**: these are
   THIS session's computational findings (all-primes L-factor reducibility + complementary trace
   matching with LMFDB isogeny classes 90.c×270.d, 54.b×90.c, 54.b×162.d at all p<100). They are
   not stated in the source papers, and no formal isogeny certificate has been written down.
   (Conversely, the D4 simplicity certificates quoted for other curves ARE proofs.)
7. **J₀(22)(Q)_tors = [5,5], J₀(26) = [21], J₀(28) = [6,6], J₀(50) = [15]**: proven here by
   Magma `TorsionSubgroup` on verified models; classical/expected via cuspidal-group theory
   (Ogg; Ohta, J. Math. Soc. Japan 2014; Yoo, "The rational torsion subgroup of J₀(N)", Adv.
   Math. 2023, arXiv:2106.01020; Yoo, "The rational cuspidal subgroup of J₀(N)",
   arXiv:2504.12564), but I did not locate an explicit literature statement of these four values
   to cite verbatim, so the literature citation is "folklore + general theorems", with the
   computation as primary evidence.
8. The identification "PZP 2013 Thm 3 order-28 curve = Howe's C₂₈,₁" rests on equality of
   G2-invariants and |disc| of the two published models [M]; the models themselves were
   transcribed from Howe 2015 (Table 2 and Thm 2.1).

---

## (d) Deliberately excluded: famous SIMPLE (non-split) constructions (side list)

| Order | Source | Status |
|---|---|---|
| 11 family | Flynn, JNT 36 (1990) | simplicity not addressed |
| 13, 15, 17, 19, 21 families | Leprévost C.R. 313 (1991) ×2; Manuscripta 75 (1992); Sém. Paris 1991–92 | simplicity not addressed in the originals (13-family widely treated as generically simple) |
| 20 family, 23-related family | Leprévost, Manuscripta 92 (1997) | not addressed |
| 22, 23, 24, 26, 29 singles | Leprévost, JTNB 7 (1995) Thm 1.2.1 + §3.2 (C₂₂, C₂₃, C₂₄,₁, C₂₆, C₂₉) | **absolutely simple** (D4 criterion; [M] confirmed) |
| 25 single | Leprévost, JTNB 7 (1995) | simple/Q; conjecturally split over Q(√3) |
| 23 | Ogawa, Proc. Japan Acad. 70 (1994) | unverified (see §c2); cannot be split/Q |
| 27 single (new #1), 39 single (new) | Howe 2015 §4 Table 2 | **absolutely simple** [M] (D4 @ 5) — not stated by Howe |
| 28 (2012a Thm 4), 33 | Platonov–Petrunin, Dokl. Math. 85/86 (2012); simplicity: Platonov–Zhgun–Petrunin, Dokl. Math. 87 (2013) | **absolutely simple** ([M] D4 @ 7, @ 13) |
| 30 family | Leprévost, unpublished (4/1996) | reported by Elkies/Howe |
| 32 family; 34 (two); 39; 40 | Elkies, g2_tors web page (2001–02) | **absolutely simple** (D4; [M] confirmed); 40 = simple record |
| (29 single) | Leprévost, JTNB (f₂₉) | absolutely simple — the 1995 simple record |

Also excluded as geometrically simple: X₀(23), X₀(29), X₀(31) (RM Jacobians, torsion [11], [7],
[5]). Excluded as non-genus-2: all of HLP's genus-3 results (Table 2; up to order 864), Flynn's
higher-genus sequences, Leprévost's genus-g families for g ≥ 3.

---

## Primary sources (with URLs)

- E. W. Howe, F. Leprévost, B. Poonen, Forum Math. 12 (2000) 315–364. arXiv:math/9809210;
  https://math.mit.edu/~poonen/papers/large.pdf; https://ewhowe.com/papers/paper11.html
- E. W. Howe, Bull. London Math. Soc. 47 (2015) 127–135. arXiv:1407.2654;
  https://doi.org/10.1112/blms/bdu107
- F. Leprévost, J. Théor. Nombres Bordeaux 7 (1995) 283–306.
  http://www.numdam.org/item/JTNB_1995__7_1_283_0/
- F. Leprévost, C. R. Acad. Sci. Paris 313 (1991) 451–454 and 771–774 (Gallica:
  ark:/12148/bpt6k57325582/f455.image, …/f775.image); Manuscripta Math. 75 (1992) 303–326
  (doi:10.1007/BF02567087); Manuscripta Math. 92 (1997) 47–63; J. reine angew. Math. 473 (1996)
  59–68 (doi:10.1515/crll.1995.473.59); C. R. 316 (1993) 819–821.
- E. W. Howe, F. Leprévost, B. Poonen, C. R. Acad. Sci. Paris 323 (1996) 1031–1034 (announcement).
- N. D. Elkies, https://people.math.harvard.edu/~elkies/g2_tors.html
- V. P. Platonov, M. M. Petrunin, Dokl. Math. 85 (2012) 286–288 (doi:10.1134/S1064562412020330);
  Dokl. Math. 86 (2012) 642–643 (doi:10.1134/S1064562412050304); V. P. Platonov, V. S. Zhgun,
  M. M. Petrunin, Dokl. Math. 87 (2013) 318–321 (doi:10.1134/S1064562413030216); V. P. Platonov,
  Russ. Math. Surv. 69:1 (2014) (survey; mathnet id rm9563).
- H. Ogawa, Proc. Japan Acad. Ser. A 70 (1994) 295–298 (doi:10.3792/pjaa.70.295).
- L. Kulesz, J. Number Theory 76 (1999) 130–146 (doi:10.1006/jnth.1998.2339).
- A. P. Ogg, "Rational points on certain elliptic modular curves", Proc. Symp. Pure Math. 24
  (1973) 221–231. (Torsion of J₁(13), J₁(16), J₁(18) at p. 226, as cited by Leprévost.)
- E. V. Flynn, J. Number Theory 36 (1990) 257–265 (doi:10.1016/0022-314X(90)90089-A);
  Invent. Math. 106 (1991) 433–442.
- Y. Hasegawa, Proc. Japan Acad. Ser. A 71 (1995) 235–239 (genus-2 X₀(N) and quotient models).
- The LMFDB Collaboration, https://www.lmfdb.org (Genus 2 curves over Q), accessed 2026-08-12.
- H. Yoo, arXiv:2106.01020 (Adv. Math. 2023), arXiv:2504.12564; M. Ohta, "Eisenstein ideals and
  the rational torsion subgroups of modular Jacobian varieties II", (squarefree level) — context
  for the J₀(N) torsion values.
- S. Frengley, arXiv:2208.05842 (12-congruences; no torsion application).

*Magma verification scripts and logs for every [M] tag: session scratchpad `lit/verify2.m`,
`lit/verify3.m` (Magma V2.29-9; `TorsionSubgroup`, D4/Galois-group certificates, L-factor trace
matching against LMFDB `ec_classdata.aplist`).*
