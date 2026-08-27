#!/usr/bin/env python3
"""make_tables.py -- regenerate table1.tex and table2.tex (the two census
tables of the paper) from the machine-readable curve data in table1.txt and
table2.txt.  Run from this directory:

    python3 make_tables.py

The .txt files are the primary source: one row per table row, in table order,
with fields separated by '|' (see the header comments in those files):

    group | [[f],[h]] | label | display | source | route | comment

  * group    -- invariant factors of the torsion subgroup, e.g. [2,2,14]
                ([] for the trivial group).
  * [[f],[h]] -- integer coefficient lists of f and h in ASCENDING degree
                (the LMFDB convention) for the model y^2 + h(x)y = f(x).
                This is the machine truth; the verification scripts read it.
  * label    -- LMFDB label: cond.class.disc.num (production, www.lmfdb.org)
                or cond.class.num (extended database, alpha.lmfdb.org), or
                '-' if the curve is in neither.  Alpha labels are snapshot
                identifiers (2026-08) and are NOT permalinks; they are kept
                as provenance but never used to build links (see below).
  * display  -- '-' for the standard expanded rendering, or a list of factor
                coefficient lists (ascending degree; a length-1 list is a
                constant factor) whose product must equal f exactly -- this
                is asserted below, so the factored presentation can never
                drift from the data.  The factor ORDER is meaningful (it
                reflects the construction) and is preserved.
  * source   -- verbatim LaTeX for the Source(s) column.
  * route    -- splitness-certificate route for this curve (table2 only:
                R1, Q2, or COVER:n with n indexing the cover data in
                verify_split_certificates.m; '-' in table1).  Curve-specific
                verification metadata, ignored here.
  * comment  -- free text (provenance of non-database curves); ignored here.

Link policy: production rows are linked by label (production labels are
permalinks).  Extended-database rows are linked as
    https://alpha.lmfdb.org/Genus2Curve/Q/?jump=[[f],[h]]
which makes the LMFDB look the curve up by equation (matching isomorphic
models), because alpha labels will change; the equation is the stable
identifier.  Rows with label '-' are not linked.
"""

import ast, os, sys

DIR = os.path.dirname(os.path.abspath(__file__)) or "."

# ---------------------------------------------------------------- data files

def read_rows(path):
    rows = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("|")
            assert len(parts) == 7, f"{path}:{ln}: expected 7 fields"
            group = ast.literal_eval(parts[0])
            f, h = ast.literal_eval(parts[1])
            label = parts[2]
            display = None if parts[3] == "-" else ast.literal_eval(parts[3])
            source = parts[4]
            assert all(isinstance(c, int) for c in f + h), f"{path}:{ln}"
            if display is not None:
                assert poly_mul_list(display) == f, \
                    f"{path}:{ln}: display factorization does not multiply out to f"
                assert h == [0], f"{path}:{ln}: factored display requires h = 0"
            rows.append(dict(group=group, f=f, h=h, label=label,
                             display=display, source=source, line=ln))
    return rows

# ------------------------------------------------------- polynomial helpers

def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        for j, cb in enumerate(b):
            out[i + j] += ca * cb
    return out

def poly_mul_list(factors):
    prod = [1]
    for fac in factors:
        prod = poly_mul(prod, fac)
    while len(prod) > 1 and prod[-1] == 0:
        prod.pop()
    return prod

# ------------------------------------------------------------- LaTeX output

def poly_tex(coeffs, var="x"):
    """Render an ascending coefficient list in descending powers."""
    terms = []
    for i in range(len(coeffs) - 1, -1, -1):
        c = coeffs[i]
        if c == 0:
            continue
        if i == 0:
            mono = ""
        elif i == 1:
            mono = var
        else:
            mono = f"{var}^{{{i}}}"
        if c == 1 and i > 0:
            cs = ""
        elif c == -1 and i > 0:
            cs = "-"
        else:
            cs = str(c)
        t = cs + mono
        if terms and not t.startswith("-"):
            t = "+" + t
        terms.append(t)
    return "".join(terms) if terms else "0"

def is_square(n):
    if n < 0:
        return False
    r = int(n ** 0.5)
    while r * r < n:
        r += 1
    return r * r == n

def factor_tex(fac):
    """Render one display factor."""
    if len(fac) == 1:                      # constant factor
        c = fac[0]
        return "-" if c == -1 else str(c)
    if fac == [0, 1]:                      # the factor x, printed bare
        return "x"
    # monic linear factor with (large) square constant: print x +- a^2
    if len(fac) == 2 and fac[1] == 1 and is_square(abs(fac[0])):
        a = int(round(abs(fac[0]) ** 0.5))
        if a >= 30:
            return f"(x{'+' if fac[0] > 0 else '-'}{a}^2)"
    return f"({poly_tex(fac)})"

def eqn_tex(row):
    f, h, display = row["f"], row["h"], row["display"]
    if display is not None:
        return "$y^2 = " + "".join(factor_tex(fac) for fac in display) + "$"
    ft = poly_tex(f)
    ht = poly_tex(h)
    if ht == "0":
        return f"$y^2 = {ft}$"
    if ht == "1":
        return f"$y^2 + y = {ft}$"
    hpar = f"({ht})" if ("+" in ht or "-" in ht[1:]) else ht
    return f"$y^2 + {hpar}y = {ft}$"

def grp_tex(group):
    if not group:
        return r"$[\,]$"
    return r"$\grp{" + ",\\,".join(str(n) for n in group) + "}$"

def coeff_str(f, h):
    return "[[" + ",".join(str(c) for c in f) + "],[" + \
           ",".join(str(c) for c in h) + "]]"

def curve_cell(row):
    eqn = eqn_tex(row)
    label = row["label"]
    if label == "-":
        return eqn
    parts = label.split(".")
    if len(parts) == 4:                    # production label: permalink
        url = "https://www.lmfdb.org/Genus2Curve/Q/" + "/".join(parts)
    else:                                  # extended DB: lookup by equation
        assert len(parts) == 3, f"bad label {label}"
        url = ("https://alpha.lmfdb.org/Genus2Curve/Q/?jump=" +
               coeff_str(row["f"], row["h"]))
    return rf"\href{{{url}}}{{{eqn}}}"

# --------------------------------------------------------------- the tables

HEADER = r"""\begingroup\scriptsize
\renewcommand{\arraystretch}{1.18}
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{1.35cm}
 >{\raggedright\arraybackslash}p{0.47\textwidth}
 >{\raggedright\arraybackslash}p{%(srcwidth)s}@{}}
\caption{All %(count)d finite abelian groups known to arise as $J(\Q)_{\mathrm{tors}}$ for a geometrically %(kind)s genus two Jacobian over \(\Q\), the
smallest-conductor known example, and sources.\label{%(tablabel)s}}\\
\toprule
Group & Curve & Source(s)\\
\midrule
\endfirsthead
\multicolumn{3}{@{}l}{{\normalsize\emph{Table~\ref{%(tablabel)s} (continued)}}}\\
\toprule
Group & Curve & Source(s)\\
\midrule
\endhead
\midrule
\endfoot
\bottomrule
\endlastfoot
"""

FOOTER = "\\end{longtable}\n\\endgroup\n"

TABLES = [
    ("table1.txt", "table1.tex",
     dict(kind="simple", tablabel="tab:census", srcwidth="1.1cm")),
    ("table2.txt", "table2.tex",
     dict(kind="split", tablabel="tab:splitcensus", srcwidth="1.25cm")),
]

def render(rows, meta):
    body = "".join(
        f"{grp_tex(r['group'])} & {curve_cell(r)} & {r['source']}\\\\\n"
        for r in rows)
    return HEADER % dict(meta, count=len(rows)) + body + FOOTER

def main():
    for txt, tex, meta in TABLES:
        rows = read_rows(os.path.join(DIR, txt))
        groups = [tuple(r["group"]) for r in rows]
        assert len(set(groups)) == len(groups), f"{txt}: duplicate group"
        out = render(rows, meta)
        with open(os.path.join(DIR, tex), "w") as fh:
            fh.write(out)
        print(f"wrote {tex}: {len(rows)} rows "
              f"({sum(1 for r in rows if r['display'] is not None)} factored, "
              f"{sum(1 for r in rows if r['label'] != '-')} linked)")

if __name__ == "__main__":
    main()
