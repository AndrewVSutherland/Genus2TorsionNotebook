#!/usr/bin/env python3
"""Generate paper/split_torsion_table.tex — torsion subgroups with a
geometrically split genus-2 Jacobian realization over Q, in the same format
as paper/torsion_realizations.tex (explicit minimal-model equations
hyperlinked to the LMFDB, Source(s) column, bibliography).

Inputs: product/data/split_min_witnesses.csv (minimal-conductor split witness
per torsion group in the extended database, exported from g2c_curves_new on
2026-08-12) plus the non-database witnesses embedded below (HLP/Howe curves
and this paper's new gluing witnesses, verified in
product/logs/verify_witnesses*.log).

Run from the repository root:  python3 code/claude_gen_split_table.py
"""
import csv, re, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "product", "data", "split_min_witnesses.csv")
OUT = os.path.join(ROOT, "paper", "split_torsion_table.tex")

# ----- production-LMFDB membership: (class, abs_disc) -> production label
# (queried from g2c_curves 2026-08-13; everything else links to alpha)
PROD = {
    ("169.a", "169"): "169.a.169.1",
    ("1083.a", "390963"): "1083.b.390963.1",
    ("196.a", "21952"): "196.a.21952.1",
    ("256.a", "512"): "256.a.512.1",
    ("294.a", "294"): "294.a.294.1",
    ("324.a", "648"): "324.a.648.1",
    ("336.a", "172032"): "336.a.172032.1",
    ("360.a", "6480"): "360.a.6480.1",
    ("363.a", "43923"): "363.a.43923.1",
    ("450.a", "36450"): "450.a.36450.1",
    ("484.a", "1936"): "484.a.1936.1",
    ("600.a", "18000"): "600.a.18000.1",
}

# ----- Source(s) column overrides, keyed by invariant tuple.
# Default: BSSVY for production-linked rows, BookerSutherland for alpha rows.
SRC = {
    (19,): r"\cite{MazurTate1973,Ogg1973}",
    (20,): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (21,): r"\cite{Ogg1973}$^\dagger$, \cite{Leprevost1995}",
    (24,): r"\cite{Leprevost1995}",
    (25,): r"\cite{Leprevost1995}$^\dagger$, \cite{PZP2013}",
    (27,): r"\cite{Leprevost1995}",
    (28,): r"\cite{Howe2015}",
    (30,): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (35,): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (36,): r"\cite{PP2012b,Platonov2014}",
    (40,): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (45,): r"\cite{HLP2000}$^\dagger$, new",
    (48,): r"\cite{PP2012b,Platonov2014}, \cite{Howe2015}",
    (60,): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (63,): r"\cite{HLP2000}",
    (70,): r"\cite{Howe2015}",
    (2, 10): r"\cite{Ogg1973}",
    (2, 12): r"\cite{HLP2000}$^\dagger$, \cite{BSSVY}",
    (2, 24): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (3, 3): r"\cite{BSSVY}",
    (3, 6): r"\cite{BSSVY}",
    (3, 9): r"\cite{Leprevost1995}$^\dagger$, \cite{BSSVY}",
    (3, 12): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (5, 5): r"\cite{HLP2000}, \cite{BookerSutherland}",
    (5, 10): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (6, 6): r"\cite{HLP2000}$^\dagger$, \cite{BSSVY}",
    (6, 12): r"\cite{HLP2000}$^\dagger$, new",
    (7, 7): r"\cite{HLP2000}",
    (8, 8): r"\cite{HLP2000}$^\dagger$, new",
    (2, 2, 24): r"\cite{HLP2000}$^\dagger$, new",
    (2, 4, 8): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (2, 6, 6): r"\cite{HLP2000}$^\dagger$, \cite{BookerSutherland}",
    (2, 2, 4, 4): r"new",
    # trivial torsion (2026-08-23): witness 1083.b.390963.1 (BSSVY).  NB the
    # alpha curve 961.a2 claims split/QxQ with trivial torsion at smaller
    # conductor, but that is an ALPHA DATA ERROR: it is isomorphic to
    # 961.a.961.1 (J0(31) class), geometrically simple with RM by sqrt 5 --
    # verified in Magma (root-power certificate at p=7), 2026-08-23.
    (): r"\cite{BSSVY}",
    (2, 2, 4, 8): r"\cite{HLP2000}$^\dagger$, new",
}

# ----- non-database witnesses: read from the shared structured source
# (product/data/split_table_witnesses.json, also consumed by
# code/claude_gen_split_table_check.py) and DERIVE the display equations from
# the coefficients, so table and verifier cannot drift.
import json
WJSON = os.path.join(ROOT, "product", "data", "split_table_witnesses.json")


def poly_tex(coeffs):
    """ascending integer coefficient list -> descending TeX polynomial"""
    terms = []
    for k in range(len(coeffs) - 1, -1, -1):
        c = coeffs[k]
        if c == 0:
            continue
        if k == 0:
            mono = ""
        elif k == 1:
            mono = "x"
        else:
            mono = "x^{%d}" % k
        if mono == "":
            body = str(abs(c))
        elif abs(c) == 1:
            body = mono
        else:
            body = str(abs(c)) + mono
        sign = "-" if c < 0 else "+"
        terms.append((sign, body))
    if not terms:
        return "0"
    s = ("-" if terms[0][0] == "-" else "") + terms[0][1]
    for sign, body in terms[1:]:
        s += sign + body
    return s


def eqn_tex(eqn):
    """LMFDB eqn string '[[f],[h]]' -> minimal-model equation TeX"""
    m = re.match(r"\[\[(.*?)\],\[(.*?)\]\]$", eqn.replace(" ", ""))
    f = [int(t) for t in m.group(1).split(",")] if m.group(1) else []
    h = [int(t) for t in m.group(2).split(",")] if m.group(2) else []
    ftex = poly_tex(f)
    if not h:
        return "y^2 = %s" % ftex
    htex = poly_tex(h)
    if len([c for c in h if c != 0]) > 1 or (len(h) > 1 and h[-1] != 1) or htex not in ("y",):
        # parenthesize unless h is a bare monomial like x^{3} or 1
        if "+" in htex[1:] or "-" in htex[1:]:
            htex = "(%s)" % htex
    if htex == "1":
        return "y^2 + y = %s" % ftex
    return "y^2 + %sy = %s" % (htex, ftex)


def group_tex(invs):
    if len(invs) == 0:
        return r"$[\,]$"
    return "$[%s]$" % ",\\,".join(str(d) for d in invs)


def parse_label(label):
    m = re.match(r"(\d+)\.([a-z]+)(\d+)$", label)
    return m.group(1), m.group(2), m.group(3)


def eq_from_lists(f, h):
    ftex = poly_tex(f)
    if not h:
        return "y^2 = %s" % ftex
    htex = poly_tex(h)
    if "+" in htex[1:] or "-" in htex[1:]:
        htex = "(%s)" % htex
    if htex == "1":
        return "y^2 + y = %s" % ftex
    return "y^2 + %sy = %s" % (htex, ftex)


NEWROWS = {}
for w in json.load(open(WJSON))["witnesses"]:
    invs = tuple(w["invs"])
    if "display" in w:
        eq = w["display"]
    elif "quintic_roots" in w:
        eq = "y^2 = x" + "".join("(x+%d)" % r for r in w["quintic_roots"])
    else:
        eq = eq_from_lists(w["f"], w["h"])
    NEWROWS[invs] = (eq, SRC[invs])


def build_lines():
    """The table body lines, rebuilt from the shared sources (CSV + JSON).
    Also imported by code/claude_gen_split_table_check.py, which compares
    the result against the shipped TeX so the verifier cannot drift."""
    rows = {}
    with open(CSV) as fh:
        for r in csv.DictReader(fh):
            invs = tuple(int(t) for t in r["torsion_subgroup"].strip("[]").split(",")) \
                if r["torsion_subgroup"] != "[]" else ()
            rows[invs] = r

    allgroups = sorted(set(rows) | set(NEWROWS), key=lambda t: (len(t), t))

    lines = []
    for invs in allgroups:
        if invs in NEWROWS and invs not in rows:
            eq, src = NEWROWS[invs]
            lines.append("%s & $%s$ & %s\\\\ \\hline" % (group_tex(invs), eq, src))
            continue
        r = rows[invs]
        cond, letter, num = parse_label(r["label"])
        cls = "%s.%s" % (cond, letter)
        eq = eqn_tex(r["eqn"])
        prod = PROD.get((cls, r["abs_disc"]))
        if prod:
            url = "https://www.lmfdb.org/Genus2Curve/Q/%s/%s/%s/%s" % tuple(prod.split("."))
            default_src = r"\cite{BSSVY}"
        else:
            # alpha labels are not permanent (AVS 2026-08-24): link via the
            # equation-lookup jump, which matches any isomorphic model
            url = ("https://alpha.lmfdb.org/Genus2Curve/Q/?jump=" +
                   r["eqn"].replace(" ", "").replace(",[]]", ",[0]]"))
            default_src = r"\cite{BookerSutherland}"
        src = SRC.get(invs, default_src)
        lines.append("%s & \\href{%s}{$%s$} & %s\\\\ \\hline"
                     % (group_tex(invs), url, eq, src))
    return lines


lines = build_lines()

PREAMBLE = r"""% split_torsion_table.tex — torsion subgroups with a geometrically split
% genus-2 Jacobian realization over Q, with minimal-conductor examples and
% earliest references.  Generated by code/claude_gen_split_table.py.
\documentclass[11pt]{amsart}
\usepackage{array}
\usepackage{longtable}
\usepackage[margin=2cm]{geometry}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
\title{Exact torsion groups realized on geometrically split\\
genus-2 Jacobians over $\mathbb{Q}$}
\date{August 2026}
\begin{document}
\maketitle
\noindent Groups are written by invariant factors: $[n_1,\ldots,n_r]$ denotes
$\mathbb{Z}/n_1\mathbb{Z}\times\cdots\times\mathbb{Z}/n_r\mathbb{Z}$ with
$n_1\mid\cdots\mid n_r$.  Each row records one genus-2 curve $C/\mathbb{Q}$
whose Jacobian is geometrically \emph{split} --- isogenous over
$\overline{\mathbb{Q}}$ to a product of two elliptic curves, i.e.\ not
geometrically simple --- and whose \emph{full} rational torsion group is
\emph{exactly} the displayed group,
$\operatorname{Jac}(C)(\mathbb{Q})_{\rm tors}\cong G$ (the first row is the
trivial group, written $[\,]$).  This is the split companion of
\texttt{torsion\_realizations.tex}; the two tables partition the known
realizations by geometric (non)simplicity.  We display, for each group, the
equation of the curve of smallest conductor in the queried extended-database
snapshot (2026-08-12) with $\operatorname{Jac}$ geometrically split; rows
whose group does not occur in the database show a construction from the
literature or from the present paper.  Curves of small conductor and
discriminant lie in the production LMFDB; the remaining database curves cite
\cite{BookerSutherland}.  For convenience each database equation is
hyperlinked to its LMFDB home page: production curves by their permanent
label, and extended-database curves via an alpha-site
(\url{https://alpha.lmfdb.org}) equation-lookup URL, since alpha labels are
not permanent (the lookup matches any isomorphic model of the curve).  The
equations are the intended stable identifiers.

\emph{Source(s) column.}  It gives the earliest construction of the torsion
group on a split genus-2 Jacobian over $\mathbb{Q}$, together with the source
of an exact-torsion witness when these differ; ``new'' means the exact split
realization first appears in the present paper (witness equations, split
certificates, and elliptic factors are recorded in
\texttt{product/data/new\_split\_witnesses.txt} and the verification logs
cited there).  A dagger ($\dagger$) marks a source that exhibits the torsion
--- typically a positive-dimensional family with $G$ \emph{contained} in the
torsion (Howe--Lepr\'evost--Poonen, Theorem~1) --- without an exact-torsion
witness.  Six rows are marked ``new'': $[45]$, whose exact witness comes from
an earlier stage of this project (verified in
\texttt{results/claude\_ov\_lane8\_verify.log}), and $[8,8]$, $[6,12]$,
$[2,2,24]$, $[2,2,4,8]$, $[2,2,4,4]$, constructed by this paper's gluing
search --- $[2,2,24]$ and $[2,2,4,8]$ by instantiating the square conditions
of \cite[\S3.7]{HLP2000}, and $[2,2,4,4]$ by gluing the elliptic curves
\texttt{210.c5} and \texttt{2310.o4} (its Jacobian has a Richelot isogeny to
their product).  All six are archived in
\texttt{product/data/new\_split\_witnesses.txt}.  Five of the six ($[45]$,
$[8,8]$, $[6,12]$, $[2,2,24]$, $[2,2,4,8]$) lie in HLP's Table~1, where they
were previously containment-only, so every group of that table now has an
exact witness.

The $[8,8]$ witness is one member of an apparently infinite exact family: on
the $X_1(8)\times X_1(8)$ gluing chart the parameter locus $t=a^2/(a^2+b^2)$
satisfies the gluing conditions identically, and all $59$ distinct instances
found to parameter height $250$ have exact torsion $[8,8]$ (HLP's family gives
the containment; the observed exactness throughout the family is new --- see
\texttt{product/logs/lane1212\_H250.log}).

Two groups occur only for $\mathbb{Q}$-simple, geometrically split Jacobians:
$[19]$ (the modular Jacobian $J_1(13)$, \cite{MazurTate1973}) and $[25]$;
every other row is split over $\mathbb{Q}$.  The largest order occurring is
$128 = \#[2,2,4,8]$ (the maximum of HLP's Table~1), and the largest cyclic
order is $70$ \cite{Howe2015}.  The $[36]$ and $[48]$ rows are the curves of
Platonov--Petrunin \cite{PP2012b}: the displayed minimal-conductor database
witnesses coincide with their curves $f_{36}$ and $f_{48,1}$ (and the second
curve of \cite{PP2012b} is the database curve \texttt{5292.c2}), with
geometric non-simplicity proven in \cite[\S7]{Platonov2014}.
\medskip

\begingroup\small
\renewcommand{\arraystretch}{1.25}
\begin{longtable}{|>{\raggedright\arraybackslash}p{1.7cm}|>{\raggedright\arraybackslash}p{0.70\textwidth}|>{\raggedright\arraybackslash}p{2.6cm}|}
\hline
Group & Curve & Source(s)\\ \hline
\endfirsthead
\hline
Group & Curve & Source(s)\\ \hline
\endhead
"""

BIB = r"""\end{longtable}
\endgroup

\begin{thebibliography}{99}
\bibitem{BSSVY} A.~R. Booker, J. Sijsling, A.~V. Sutherland, J. Voight, and D. Yasaki,
\emph{A database of genus-2 curves over the rational numbers},
LMS J.~Comput.\ Math.\ \textbf{19(A)} (2016), 235--254.
\bibitem{BookerSutherland} A.~R. Booker and A.~V. Sutherland,
\emph{Genus 2 curves of small conductor}, in preparation.
\bibitem{HLP2000} E.~W. Howe, F. Lepr\'evost, and B. Poonen, \emph{Large torsion
subgroups of split Jacobians of curves of genus two or three}, Forum Math.\
\textbf{12} (2000), 315--364.
\bibitem{Howe2015} E.~W. Howe, \emph{Genus-2 Jacobians with torsion points of
large order}, Bull.\ London Math.\ Soc.\ \textbf{47} (2015), 127--135.
\bibitem{Leprevost1995} F. Lepr\'evost, \emph{Jacobiennes de certaines courbes de
genre 2: torsion et simplicit\'e}, J.~Th\'eor.\ Nombres Bordeaux \textbf{7} (1995),
283--306.
\bibitem{MazurTate1973} B. Mazur and J. Tate, \emph{Points of order 13 on elliptic
curves}, Invent.\ Math.\ \textbf{22} (1973), 41--49.
\bibitem{Ogg1973} A. Ogg, \emph{Rational points on certain elliptic modular curves},
Proc.\ Sympos.\ Pure Math.\ XXIV, Amer.\ Math.\ Soc., 1973, 221--231.
\bibitem{Platonov2014} V.~P. Platonov, \emph{Number-theoretic properties of
hyperelliptic fields and the torsion problem in Jacobians of hyperelliptic curves
over the rational number field}, Russian Math.\ Surveys \textbf{69} (2014), no.~1,
1--34.
\bibitem{PP2012b} V.~P. Platonov and M.~M. Petrunin, \emph{On the torsion problem
in Jacobians of curves of genus 2 over the rational number field},
Dokl.\ Math.\ \textbf{86} (2012), 642--643.
\bibitem{PZP2013} V.~P. Platonov, V.~S. Zhgun, and M.~M. Petrunin,
\emph{On the simplicity of Jacobians for hyperelliptic curves of genus 2 over the
field of rational numbers with torsion points of high order},
Dokl.\ Math.\ \textbf{87} (2013), 318--321.
\end{thebibliography}
\end{document}
"""

if __name__ == "__main__":
    with open(OUT, "w") as fh:
        fh.write(PREAMBLE)
        fh.write("\n".join(lines) + "\n")
        fh.write(BIB)
    print("wrote %s (%d rows)" % (OUT, len(lines)))
