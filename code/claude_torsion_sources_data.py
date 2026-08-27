# claude_torsion_sources_data.py — the per-group dated & classified source matrix
# (Drew's protocol, 2026-07-21: for each torsion subgroup, EVERY paper that mentions
# it in any way, sorted by date, each classified).
#
# Simplicity categories:
#   GS    = explicitly geometrically/absolutely simple (proof or certificate in source)
#   QS    = explicitly simple over Q only
#   SPLIT = explicitly split / decomposable construction
#   NA    = simplicity not addressed
# CONSTRUCTS describes what the source actually builds (curve/family; pointwise
# rational torsion vs Galois-stable subgroup only).
#
# Rows verified in: notes/claude_torsion_refs_dossier.md (first pass, line-verified
# against nicholls.txt), Sol's review 2026-07-21, and the three-agent sweep of
# 2026-07-21 (see notes/claude_torsion_sources_audit.md).  UNVERIFIED rows are
# marked verified=False and rendered with a * in the companion document.
#
# Each row: (year, refkey, constructs, simplicity, evidence, verified)

SOURCES = {
  # ---- filled/extended by the 2026-07-21 agent sweep; skeleton rows below ----
  "[5]": [
    (2000, "BoxallGrant2000", "single CM curve y^2-y=x^5, exact [5]", "GS",
     "TAMS 352, torsion computed; absolute simplicity stated", True),
    (2001, "BGL2001", "5-torsion families (2-dim)", "NA",
     "JLMS 64, families with a 5-torsion point", True),
    (2024, "Elkies2024", "3-dim 5-torsion moduli + families", "GS",
     "LuCaNT; generic End = Z for the (2)^4+5 family", True),
  ],
  "[11]": [
    (1990, "Flynn1990", "1-param family, rational 11-torsion point", "NA",
     "JNT 36; simplicity not discussed", True),
    (2018, "DS2018", "continued-fraction families, 11-torsion", "NA",
     "JNT 189", True),
  ],
  "[13]": [
    (1990, "Flynn1990", "13-torsion construction", "NA", "JNT 36 (with Fly91)", True),
    (1991, "Leprevost1991a", "1-param family, order-13 divisor class", "NA",
     "C.R. 313, 451-454", True),
  ],
  "[3,3]": [
    (2014, "BFT2014", "rational parametrization of genus-2 curves with (3,3)-structure",
     "NA", "Acta Arith. 165; pointwise-rationality and simplicity per agent verification",
     False),
    (2024, "LSSV2024", "exact (Z/3)^2 on a PQM Jacobian (non-generic)", "GS",
     "Forum Math. Sigma 12, Table 2, D=6", True),
  ],
  # ... remaining groups are merged in from the agent reports by
  # claude_gen_torsion_sources.py's MERGE section before rendering.
}

REFS_BIB = {
  # refkey -> full bibliography entry (LaTeX); merged with the main table's bib.
  "BFT2014": r"N. Bruin, E.~V. Flynn, and D. Testa, \emph{Descent via (3,3)-isogeny on Jacobians of genus 2 curves}, Acta Arith.\ \textbf{165} (2014), 201--223.",
  "DS2018": r"K. Daowsud and T.~A. Schmidt, \emph{Continued fractions of certain Laurent series and torsion on Jacobians}, J.~Number Theory \textbf{189} (2018), 115--130.",
}
