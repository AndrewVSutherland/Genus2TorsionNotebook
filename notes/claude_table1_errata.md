# Errata / update list for NotesAndTodo.tex Table 1 — 2026-07-20

Found while building paper/torsion_realizations.tex (production-LMFDB census
cross-check, results of run_sql on g2c_curves with is_simple_geom):

1. Z/18 row links 388.a.776.1 — that is the Z/21 curve (also linked on the
   Z/21 row). The minimal-conductor geometrically simple [18] curve in
   production is 1180.a.18880.1.
2. Z/39 row has no reference, but Elkies' 2002 page exhibits an order-39 curve
   (y^2 = x^6+4x^4+10x^3+4x^2-4x+1, differences of nonrational Weierstrass
   points) — should cite elkies2002 like the [32]/[34]/[40] rows.
3. Z/9 row links 745.a.745.1; minimal-conductor production example is
   713.b.713.1 (cosmetic).
4. Missing rows for groups realized by this project (no LMFDB entry):
   [2,2,2,12] (three curves), [6,6], [2,4,8], [2,2,20], [2,4,4], [2,2,4,4].
   ([2,2,2,8] already has its alpha row 3942.b.3.)
5. The (4,16) row: confirmed not geometrically simple (bielliptic; see
   m18_416 notes) — should be removed or marked split.
6. Elkies' page also notes a rational 31-SUBGROUP with no rational generator
   (not a rational-point realization; worth a footnote only).
