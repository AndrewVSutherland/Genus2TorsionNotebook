# Verification scripts for "Rational Torsion on Simple Genus Two Jacobians"

Every claim of the paper that rests on a computation (rather than on an
argument given in the paper or a citation to the literature) is verified by
one of the eight Magma scripts in this directory.  The displayed curves,
torsion groups and sources of the paper's two census tables live in the
machine-readable data files `table1.txt` / `table2.txt`; the LaTeX tables
`table1.tex` / `table2.tex` are typeset from them by `make_tables.py`, and
the verification scripts read the same files, so the typeset tables and the
verified claims cannot drift apart.  All other data a script needs is
embedded in the script.  Every script asserts everything it claims (any
failure aborts with a runtime error) and ends with an
`ALL ... VERIFIED` / `ALL ... PASSED` line.  Total runtime is under a minute
on a 2026 desktop with Magma V2.29-9 (N.B. **Magma V2.28-9 or later is required**).

Run everything:

```
for f in verify_*.m; do magma -b "$f" || break; done
```

## Table data and links

`table1.txt` and `table2.txt` contain one line per table row, in table
order, with fields separated by `|`:

```
group | [[f],[h]] | label | display | source | route | comment
```

`group` is the list of invariant factors of the torsion subgroup, and
`[[f],[h]]` are the integer coefficient lists of f and h **in ascending
degree** (the LMFDB convention) for the displayed model y² + h(x)y = f(x);
`display`, when present, records a factored presentation of f as a list of
factor coefficient lists whose product is asserted (by `make_tables.py`
and by the table scripts) to equal f exactly.  `route` records, per curve, the splitness-certificate
route used by `verify_split_certificates.m` (`R1`, `Q2`, or `COVER:n`
with `n` indexing the explicit cover data in that script; `-` in
`table1.txt`) -- it is curve-specific, not a function of the torsion
group.  See the header of `make_tables.py` for the full field spec.  Running `python3 make_tables.py`
regenerates `table1.tex` and `table2.tex`.

Rows whose curve is in the production LMFDB are linked by label
(production labels are permalinks).  Rows found only in the extended
Booker–Sutherland database are linked as
`https://alpha.lmfdb.org/Genus2Curve/Q/?jump=[[f],[h]]`, which asks the
LMFDB to look the curve up by its equation (this matches any isomorphic
model), because labels on the alpha site are not permanent.  The equations
are the intended stable identifiers; the `label` column is a dated
database-snapshot identifier (2026-08) kept as provenance.

## What each script verifies

| Script | Paper item | Claim checked | Time |
|---|---|---|---|
| `verify_simple_torsion_table.m` | **Table 1** | the displayed curve of each of the 73 rows (including the trivial-torsion row 277.a.3) has exactly the displayed rational torsion subgroup | 2 s |
| `verify_split_torsion_table.m` | **Table 2** | the displayed curve of each of the 77 rows (including the trivial-torsion row 1083.b.390963.1) has exactly the displayed rational torsion subgroup | 3 s |
| `verify_simple_certificates.m` | **Table 1**, Remark 6.1, §"[2,22]-torsion" | every Table 1 Jacobian (73 rows) is geometrically simple (good-prime Frobenius root-power certificate, witness prime and χ printed); End(Jac) = Z over Q̄ for **all 71 non-RM rows**, primarily by Zywina's two-Frobenius-polynomial criterion (ordinary primes q, p with the root groups Φ free of rank 3 — so the Frobenius torus is maximal and Mumford–Tate holds by Larsen–Pink — and [L(W_p):L] = 8, forcing Weyl order 8, hence MT = GSp₄, which also re-proves geometric simplicity), with an odd conductor-exponent-1 prime (semistable toric rank 1) and disjoint-splitting-fields + #tors > 18 as fallbacks; for the RM rows [31] and [2,22], a rigorous upper bound End⁰ ∈ {Q, Q(√2)} resp. {Q, Q(√5)} (two strict primes with non-isomorphic degree-8-closure quartic Frobenius fields bound the center in the unique real quadratic subfield by CLV, excluding CM and imaginary quadratic; #tors > 18 excludes QM), plus Frobenius evidence for the cited RM lower bound; and, for both [2,22] curves: exact torsion, the two rational Weierstrass points, and a rational point P with [P − W_i] of order 22 for both i; the certified curves are read directly from `table1.txt` (each certificate names its curve by its `[[f],[h]]` coefficients) | 45 s |
| `verify_split_certificates.m` | **Table 2** | geometric splitness of **all 77 rows, each by an explicit geometric object**: 59 by a degenerate Richelot codomain (a (2,2)-isogeny to a product over Q, over a quadratic field — including the new [11] row's E×E^σ over Q(√11) — or over Q(√2) at depth 2), and all 18 odd-glued curves by an explicit nonconstant map to a genus-1 curve — 13 degree-3 covers to named Cremona curves over Q (including Howe's [70] → 858k1), covers of degree 5 and 3 over quadratic fields, degree-2 covers over cubic fields, and a degree-2 cover over a sextic field for [19] = J₁(13) — each verified by a single exact polynomial identity N(p,q)·q = g·h² plus nonvanishing of the Wronskian and disc(E); the certified curves are read directly from `table2.txt` (each certificate names its curve by its `[[f],[h]]` coefficients) | 8 s |
| `verify_22212_theorem.m` | **Theorem 3.2** | the three listed points lie on S°; each curve X_P has torsion exactly [2,2,2,12]; each Jacobian is geometrically simple; the three curves are pairwise non-isomorphic | 1 s |
| `verify_22212_surface.m` | §"Geometry of S", §"How we found S°" | S has exactly 36 singular points, all rational ordinary double points with coordinates in {0,1,−1}, forming the G-orbits of [1:0:0:1:0] (12) and [1:1:0:1:1] (24); the displayed birational map sends S into the surface Z | 1 s |
| `verify_2220_family.m` | §"A second order-80 group" | over Q(t) the family has classes of exact order 5 and 20; the quartic-root condition is parameterized by t(z); the cubic-root cover is a genus-2 curve with normalization Y : y² + (x²+1)y = x⁵ − x of conductor 2528; rank Jac(Y) = 1 (RankBounds = [1,1]); Chabauty provably computes Y(Q) (9 points); transporting back, z = −1/7 and z = −7/9 are the only nondegenerate rational fibers, both isomorphic to Eq. (12) with exact torsion [2,2,20] | 3 s |
| `verify_moduli_identities.m` | Lemma 2.3(b)(d), Lemma 4.5, Prop. 4.4, Cor. 4.3, Thm. 4.7, **Thm. 3.3** | Zarhin's halving formula 2D₀ = [(0,0) − ∞] generically; the Lemma 2.3(d) expressions are exactly q(−a_i²); the E_u ≅ V_t isomorphism; all Prop. 4.4 identities over Q(s) incl. the four Lemma 2.3(c) square conditions (plus an exact fiber); rank E₇ = 1; the Disc(q) = K3-quartic identity and the rational curve on the [2,2,2,8] surface (plus two exact fibers), with the density mechanism of that section's theorem: the generic fiber S ∩ {d = tc} is genus 1 over Q(t), the displayed elliptic model 𝒜_t carries the section (1,1) of infinite order, and at sample fibers Jac(S_t) ≅ 𝒜_t; and the Thm. 3.3 infinitude ingredients: on H : u−v−a+b = 0 the quartic splits **modulo the quadric** as (smooth quadric)·(a+v)(b−v), the plane pair cutting double lines in the degenerate locus {c = 0} and the smooth quadric cutting a smooth genus-1 (2,2)-curve through P₁ whose Jacobian is E = 288b1 itself; E : y² = x³−21x−20 has conductor 288 and rank exactly 1 with Q = (−3,4) generating E(Q)/tors; ψ maps E onto that curve identically over Q(E); ψ(2Q) = P₁ exactly, and ψ(4Q) ≠ P₁ lies in S° (ψ nonconstant) | 1 s |

## What rests on citations rather than on these scripts

* The RM *lower bound* for [31] and [2,22] — that the real-multiplication
  endomorphism actually exists ([31]: curve attached to newform 1830.2.a.q,
  Costa et al. dataset; [2,22]: extended-database certified endomorphisms).
  This is the one endomorphism statement that no amount of L-polynomial
  data can certify (finitely many Frobenius polynomials never prove an
  endomorphism exists; cf. Zywina, Remark 1.14(ii)); the script certifies
  everything else — simplicity, End⁰ ∈ {Q, F} with F the claimed real
  quadratic field — and prints the constant discriminant-core evidence
  (√2 resp. √5 at ~40 primes).
* "Smallest known conductor" and "only known examples" statements are
  database-scope statements, dated in the paper.
* Theorem 1.1, Theorem 3.1, Lemmas 2.1/2.2, Prop. 4.2, Cor. 4.6 and the
  adjunction/general-type argument for S are proved in the paper or cited.
