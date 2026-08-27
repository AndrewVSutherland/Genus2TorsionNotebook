#!/usr/bin/env python3
# claude_gen_torsion_table.py — generate paper/torsion_realizations.tex:
# torsion subgroups with geometrically simple genus-2 realizations /Q,
# minimal-conductor example (LMFDB-linked when possible), earliest reference.
# Reference conventions (AVS 2026-07-20): earliest literature where one exists;
# www.lmfdb.org + no pre-2016 literature -> BSSVY; alpha-only + no literature ->
# Booker--Sutherland (in prep); new in this paper -> "new".

# production LMFDB (www) minimal-conductor geometrically simple representatives
# (from g2c_curves: is_simple_geom, 2026-07-20 query)
WWW = {
 "[2]":  ("295.a.295.2",  [[-608,389,22,-40,0,1],[1,1,1]]),
 "[3]":  ("997.b.997.1",  [[0,0,-1,2,-2,1],[1]]),
 "[4]":  ("1070.a.2140.1",[[0,-1,0,1],[1,0,0,1]]),
 "[5]":  ("277.a.277.2",  [[-6,11,-19,14,-9,1],[1]]),
 "[6]":  ("1038.a.1038.1",[[3,21,46,26,-12,1],[0,1,1]]),
 "[7]":  ("461.a.461.1",  [[-2,3,0,-3,0,1],[0,0,0,1]]),
 "[8]":  ("464.a.464.1",  [[0,0,0,-1,-2,-2,-1],[1,1]]),
 "[9]":  ("713.b.713.1",  [[0,0,0,0,-1],[1,1,0,1]]),
 "[10]": ("389.a.389.1",  [[7,16,0,-8,-2,1],[0,1,0,1]]),
 "[11]": ("353.a.353.1",  [[0,0,1],[1,1,0,1]]),
 "[12]": ("762.a.3048.1", [[1,1,1],[0,1,1,1]]),
 "[13]": ("349.a.349.1",  [[0,0,-1,-1],[1,1,1,1]]),
 "[14]": ("249.a.249.1",  [[0,1,1],[1,0,0,1]]),
 "[15]": ("277.a.277.1",  [[0,-1,-1],[1,1,1,1]]),
 "[16]": ("830.a.6640.1", [[1,1,-2,0,1,-1],[1,0,0,1]]),
 "[17]": ("1996.b.510976.1",[[1,3,1,1],[0,1,1,1]]),
 "[18]": ("1180.a.18880.1",[[0,2,4,0,-2],[1,0,0,1]]),
 "[20]": ("394.a.3152.1", [[0,0,0,0,0,-1],[1,1]]),
 "[21]": ("388.a.776.1",  [[0,1,2,0,-1],[1,1,0,1]]),
 "[22]": ("1192.a.19072.1",[[1,-1,-2,1],[0,1,0,1]]),
 "[24]": ("1908.a.183168.1",[[0,2,4,3,2],[1,0,0,1]]),
 "[27]": ("604.a.9664.2", [[0,-1,1,1,-1],[1,0,0,1]]),
 "[28]": ("249.a.6723.1", [[2,3,1,1,0,-1],[1,0,0,1]]),
 "[29]": ("976.a.999424.1",[[0,0,-1,2,0,-2,1],[1,1]]),
 "[39]": ("1116.a.214272.1",[[0,-1,1,2,1],[1,0,0,1]]),
 "[2,2]": ("464.a.29696.2",[[0,1,16,72,33,4],[0,1]]),
 "[2,4]": ("997.a.997.1", [[0,-1,0,16,-8,1],[0,1]]),
 "[2,6]": ("704.a.45056.1",[[0,0,-2,-1,4,4],[1]]),
 "[2,8]": ("464.a.29696.1",[[0,0,-2,-4,3,8],[1,1]]),
 "[2,10]":("555.a.8325.1",[[0,1,1,-4,-2,3],[1,1]]),
 "[2,12]":("762.a.82296.1",[[0,-1,2,14,-8,1],[0,1,1]]),
 "[2,14]":("1416.b.135936.1",[[1,1,0,-1,-2],[0,1,0,1]]),
 "[2,2,2]":("2600.a.338000.1",[[0,1,-3,-5,8,10],[0,1]]),
 "[2,2,4]":("3978.a.930852.1",[[0,6,-8,-3,3,1],[0,1,1]]),
 "[2,2,6]":("816.a.39168.1",[[0,1,-1,-4,0,3],[1,0,1]]),
}

# alpha-only rows (only in the extended Booker-Sutherland database, not production).
# Minimal-conductor GENERIC (geom_end_alg='Q') geometrically simple witness per group,
# queried from g2c_curves_new (LMFDB mirror, 2026-07-22); label is the alpha 3-part
# URL form cond.class.num; eqn = [[f-coeffs],[h-coeffs]] for y^2+h(x)y=f(x).
# EXCEPTIONS: [2,22] and [2,2,14] have NO generic witness in the extended database --
# their only geometrically simple curves have real multiplication (flagged RM below).
#
# CONDUCTOR PROVENANCE (2026-07-25).  Every conductor in the table body is the
# LMFDB label's first component, i.e. published data; this script computes no
# conductor itself.  The ONE hard-coded conductor is the [2,2,14] value in the
# endomorphism-refinement paragraph of the preamble.  It must NOT be taken from
# PARI genus2red: genus2red is documented for p > 2 and, whenever its returned
# factorisation matrix carries the sentinel row [2,-1], it omits the 2-part of
# N entirely.  That bug put an odd 9550095752925 in this file until 2026-07-25;
# the correct value, from Magma's Conductor on the reduced minimal Weierstrass
# model, is 2^2 * 9550095752925 = 38200383011700.  See
# notes/claude_ov_erratum_2026_07_25.md.
ALPHA = {
 # trivial exact torsion (added 2026-08-23, AVS): min-conductor geometrically
 # simple witness with trivial torsion; production's smallest is 461.a.461.2,
 # but the extended DB has this third member of the class 277.a (whose other
 # members are the [15] and [5] rows), generic, conductor 277.
 "[]":     ("277.a.3",    [[-14580,51520,-63684,31054,-3823,-637,-18],[1]]),
 "[19]":   ("1468.b.1",   [[-2,5,-3,-5,1,1,1],[1,1,1]]),
 "[23]":   ("1696.c.1",   [[1,-2,2,1,-2,-1,1],[0,0,1]]),
 "[25]":   ("168750.e.2", [[2,-23,41,75,25,-9],[1,1]]),
 "[26]":   ("11016.g.1",  [[0,0,-3,-3,-3,-12],[1,1]]),
 "[30]":   ("8694.e.1",   [[0,-1,16,-12,-10,4,2],[0,0,1,1]]),
 "[32]":   ("66852.b.1",  [[-47,-49,49,45,1,-1,1],[0,1,1]]),
 # [2,4,4]: min-conductor generic witness of the four in g2c_curves_new
 # (172260.bj3, 180180.cg7, 727650.qf8, 956340.fc7 -- all geom_end_alg='Q').
 # Verified independently in Magma 2026-07-25: exact TorsionSubgroup [2,4,4],
 # Conductor 172260, strict End=Z certificate at p0=37, q0=47 (compositum
 # 64 = 8*8, linearly disjoint).  See notes/claude_ov_db_histogram_2026_07_25.md.
 "[2,4,4]": ("172260.bj.3", [[0,180,-210,-1850,2204,-116],[0,1,1]]),
 "[33]":   ("13716.a.1",  [[2,13,21,7,9,-7,1],[1,1,1]]),
 "[34]":   ("17856.m.1",  [[15,12,7,3,0,-2],[0,0,1,1]]),
 "[36]":   ("59040.r.3",  [[9,6,-14,13,30,3],[0,0,1]]),
 "[40]":   ("1206.b.1",   [[0,0,-2,-2,4,3],[1,1]]),
 "[2,16]": ("4392.a.1",   [[8,4,-3,1,-1,-1],[0,0,1,1]]),
 "[2,18]": ("4428.a.2",   [[7,2,1,3,-5,-1,1],[0,1,1]]),
 "[2,20]": ("936600.b.1", [[0,-180,-251,185,28,-30,4],[0,0,1]]),
 "[2,26]": ("8136.c.1",   [[0,0,2,3,-12,-3,9],[1,0,1]]),
 "[2,28]": ("28200.e.1",  [[0,-3,4,34,-42,-3,9],[1,0,1]]),
 "[2,2,8]":  ("3942.b.3",   [[0,-1,-6,-3,20,9],[0,1,1]]),
 "[2,2,10]": ("10512.n.1",  [[0,0,-7,12,6,-12],[1,0,1]]),
 "[2,2,12]": ("12300.e.2",  [[0,-5,-5,16,8,-7,1],[0,1,1]]),
 "[2,2,2,2]": ("1575.c.1",  [[0,-1,-4,8,23,-27],[0,1,1]]),
 "[2,2,2,4]": ("9450.b.2",  [[0,1,-6,-16,146,-225],[0,1,1]]),
 "[2,2,2,6]": ("39600.dq.1",[[2,-12,11,18,-14,-6],[1,0,1]]),
 "[3,3]": ("5100.a.1",  [[3,3,7,2,4,0,1],[1]]),
 "[3,6]": ("3780.b.1",  [[21,-2,7,-5,0,-1],[1,1,1,1]]),
 "[3,9]": ("8436.a.2",  [[326,-75,-230,-113,40,75,16],[1,0,0,1]]),
 "[4,4]": ("1575.c.14", [[20493,6534,10832,2601,1823,242,99],[0,1,1]]),
 "[4,8]": ("913500.cf.3",[[164,682,450,-435,272,-81],[0,1,1]]),
}
# non-generic (real multiplication) rows: the ONLY known geom-simple witness has RM
ALPHA_RM = {
 "[2,22]":  ("19044.h.2",   [[0,-6,12,-5,9,-3,1],[0,1,1]]),
}
# RM-only rows NOT in any LMFDB database (conductor exceeds the 2^20 bound of the
# extended database).  [31]: the curve attached to newform 1830.2.a.q (RM by
# Z[sqrt2], conductor 1830^2 = 3348900), from the Costa et al. ModularAbelianSurfaces
# dataset (olddata/label_to_curve_QQ.txt, commit 687c10d); exact torsion Z/31Z,
# D4 simplicity certificate at p=11, verified 2026-07-30
# (notes/claude_z31_rm_witness_2026_07_30.md).  First EXPLICIT order-31
# realization and first on a genus-2 Jacobian (existence on some inexplicit
# class member follows from AC2026 Thm 4.2/4.3); irrational-support generator.
NEW_RM = {
 "[31]": [[-126,816,-2466,4300,-4587,2841,-839],[0,-1,-1]],
}

# curves realized in the present project (no LMFDB entry); equations from repo data
NEW = {
 # 2026-07-24 presentations: same curves as the certificates file except
 # [2,2,2,12], whose witness was REPLACED by Ari Shnidman's split model
 # (different curve, independently certified; notes/ari_22212_split_model_2026_07_23.md).
 # The factored forms are exact identities with the old displays
 # (notes/claude_table_row_rewrites_2026_07_23.md); keep them in sync with
 # code/claude_endz_certificates.m.
 "[2,2,14]":   r"y^2 = (x+1)(x+9)(2x-115)(6x+55)(3x^2-220x+2652)",
 "[2,2,20]":   r"y^2 = -(x-1)(6x+1)(2x+7)(6217x+1008)(21x^2-161x+144)",
 "[2,4,8]":    r"y^2 = x(3x-5)(5x+333)(16x^2+23x+576)",
 "[6,6]":      r"y^2 = x(39x^2-69x+125)(48x^2+8x+39)",
 # [2,4,4] MOVED to ALPHA on 2026-07-25: the extended database contains four
 # geometrically simple curves with exact torsion [2,4,4], all geom_end_alg='Q',
 # so this row is NOT first-realized here.  The project's own M(2,4,4) pilot
 # curve y^2 = x(x+16)(x+(644/799)^2)(x^2+x+16) remains a valid independent
 # realization and keeps its certificate in data/claude_endz_certificates.txt,
 # but it is no longer the displayed witness.
 "[2,2,2,8]":  r"y^2 = x(x+1)(x+55^2)(x+99^2)(x+125^2)",
 "[2,2,4,4]":  r"y^2 = x(x+36^2)(x+57^2)(x+64^2)(x+132^2)",
 "[2,2,2,12]": r"y^2 = x(x-120^2)(x-143^2)(x-266^2)(x-218^2)(x-241^2)",
}
# literature-only realization (not in LMFDB, not this paper)
LIT_ONLY = { "[2,2,2,10]": r"y^2 = x(x+1)(x-1)(3x-7)(8x-13)(24x+25)" }

# Earliest references, reconciled against the 2026-07-20 literature-mining pass
# (nicholls.txt line-verified; see notes/claude_torsion_refs_dossier.md).
# \dag = the cited source exhibited the torsion but geometric simplicity of a
# realizing curve was established later.
DAG = r"$^\dagger$"
REFS = {
 "[]":   r"\cite{BSSVY}",  # realized in the production DB (e.g. 461.a.461.2)
 # cyclic
 "[5]":  r"\cite{Ogg1973}" + DAG + r", \cite{BoxallGrant2000}$^\ddagger$",
 "[6]":  r"\cite{Flynn1991}" + DAG,
 "[7]":  r"\cite{Ogg1973}" + DAG + r", \cite{LPS2004}",
 "[9]":  r"\cite{Flynn1991}" + DAG,
 "[10]": r"\cite{Flynn1991}" + DAG,
 "[11]": r"\cite{Ogg1973}" + DAG + r", \cite{BLP2009}",
 "[13]": r"\cite{Flynn1990,Leprevost1991a}" + DAG + r", \cite{PZP2013}",
 "[14]": r"\cite{PP2012,PZP2013}",
 "[15]": r"\cite{Leprevost1991b}" + DAG,
 "[17]": r"\cite{Leprevost1991b}" + DAG + r", \cite{PZP2013}",
 "[18]": r"\cite{PP2012,PZP2013}",
 "[19]": r"\cite{Leprevost1991b}" + DAG + r", \cite{PZP2013}",
 "[20]": r"\cite{Leprevost1997}" + DAG,
 "[21]": r"\cite{Leprevost1991b}" + DAG + r", \cite{Nicholls2018}",
 "[22]": r"\cite{Leprevost1993}" + DAG + r", \cite{Leprevost1995}",
 "[23]": r"\cite{Ogawa1994}" + DAG + r", \cite{Leprevost1995}",
 "[24]": r"\cite{Leprevost1993}" + DAG + r", \cite{Leprevost1995}",
 "[25]": r"\cite{Leprevost1993}" + DAG + r", \cite{Nicholls2018}",
 "[26]": r"\cite{Leprevost1993}" + DAG + r", \cite{Leprevost1995}",
 "[27]": r"\cite{Howe2015}" + DAG + r", \cite{Nicholls2018}",
 "[28]": r"\cite{PP2012,PZP2013}",
 "[29]": r"\cite{Leprevost1993}" + DAG + r", \cite{Leprevost1995}",
 "[31]": r"\cite{CEHJMPV}" + DAG + r", new",
 "[30]": r"\cite{PP2015}",
 "[36]": r"\cite{PP2015}",
 "[32]": r"\cite{Elkies2002}", "[34]": r"\cite{Elkies2002}",
 "[39]": r"\cite{Elkies2002}", "[40]": r"\cite{Elkies2002}",
 "[33]": r"\cite{PZP2013}",
 # non-cyclic
 "[3,3]": r"\cite{BFT2014}",
 "[3,9]": r"\cite{Leprevost1995}" + DAG + r", \cite{BookerSutherland}",
 "[2,2,14]": r"\cite{BookerSutherland}, new",
 "[2,2,2,10]": r"\cite{Elkies2024}",
}

BSSVY = r"\cite{BSSVY}"
BS    = r"\cite{BookerSutherland}"
NEWREF = "new"

def grp_tex(g):
    if g == "[]":
        return r"$[\,]$"
    return "$" + g.replace(",", ",\\,") + "$"

def poly_tex(coeffs, var="x"):
    terms = []
    for i in range(len(coeffs)-1, -1, -1):
        c = coeffs[i]
        if c == 0: continue
        if i == 0: mono = ""
        elif i == 1: mono = var
        else: mono = f"{var}^{{{i}}}"
        if c == 1 and i > 0: cs = ""
        elif c == -1 and i > 0: cs = "-"
        else: cs = str(c)
        t = cs + mono
        if terms and not t.startswith("-"): t = "+" + t
        terms.append(t)
    return "".join(terms) if terms else "0"

def eqn_tex(eqn):
    f, h = eqn
    ht = poly_tex(h)
    ft = poly_tex(f)
    if ht == "0": return f"$y^2 = {ft}$"
    if ht == "1": return f"$y^2 + y = {ft}$"
    hpar = f"({ht})" if ("+" in ht or "-" in ht[1:]) else ht
    return f"$y^2 + {hpar}y = {ft}$"

import sys, os
# link mode: "all" (default; equation linked to its LMFDB page, www or alpha),
# "www" (link production rows only -- alpha labels are impermanent), or
# "none" (no links -- cleanest for the published version).  Set via env LINKS.
LINKS = os.environ.get("LINKS", "all")

def www_url(label):
    c, cl, d, n = label.split(".")
    return f"https://www.lmfdb.org/Genus2Curve/Q/{c}/{cl}/{d}/{n}"

def alpha_url(eqn):
    # alpha labels are not permanent (AVS 2026-08-24): link via the
    # equation-lookup jump, which matches any isomorphic model
    f, h = eqn
    return ("https://alpha.lmfdb.org/Genus2Curve/Q/?jump=[[" +
            ",".join(str(c) for c in f) + "],[" +
            ",".join(str(c) for c in (h if h else [0])) + "]]")

def linked(eqn_str, url, kind):
    # kind in {"www","alpha",None}; wrap the equation in \href per LINKS mode
    if url is None: return eqn_str
    if LINKS == "none": return eqn_str
    if LINKS == "www" and kind != "www": return eqn_str
    return rf"\href{{{url}}}{{{eqn_str}}}"

def sortkey(g):
    if g == "[]":
        return (0, [])
    v = [int(t) for t in g.strip("[]").split(",")]
    return (len(v), v)

def new_eqn(latex):
    # NEW / LIT_ONLY equations are stored as raw LaTeX bodies (no $...$)
    return f"${latex}$"

rows = []
allg = set(WWW) | set(ALPHA) | set(ALPHA_RM) | set(NEW) | set(NEW_RM) | set(LIT_ONLY)
for g in sorted(allg, key=sortkey):
    if g in WWW:
        label, eqn = WWW[g]
        curve = linked(eqn_tex(eqn), www_url(label), "www")
        ref = REFS.get(g, BSSVY)
    elif g in ALPHA:
        label, eqn = ALPHA[g]
        curve = linked(eqn_tex(eqn), alpha_url(eqn), "alpha")
        ref = REFS.get(g, BS)
    elif g in ALPHA_RM:
        label, eqn = ALPHA_RM[g]
        curve = linked(eqn_tex(eqn), alpha_url(eqn), "alpha")
        ref = REFS.get(g, BS) + r"$^{\rm RM}$"
    elif g in NEW_RM:
        curve = eqn_tex(NEW_RM[g])
        ref = REFS.get(g, NEWREF) + r"$^{\rm RM}$"
    elif g in LIT_ONLY:
        curve = new_eqn(LIT_ONLY[g])
        ref = REFS.get(g, "??")
    else:
        curve = new_eqn(NEW[g])
        ref = REFS.get(g, NEWREF)
    rows.append((grp_tex(g), curve, ref))

header = r"""% torsion_realizations.tex — torsion subgroups with a geometrically simple
% genus-2 Jacobian realization over Q, with minimal-conductor examples and
% earliest references.  Generated by code/claude_gen_torsion_table.py.
\documentclass[11pt]{amsart}
\usepackage{array}
\usepackage{longtable}
\usepackage[margin=2cm]{geometry}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
\title{Exact torsion groups realized on geometrically simple\\
genus-2 Jacobians over $\mathbb{Q}$}
\date{July 2026}
\begin{document}
\maketitle
\noindent Groups are written by invariant factors: $[n_1,\ldots,n_r]$ denotes
$\mathbb{Z}/n_1\mathbb{Z}\times\cdots\times\mathbb{Z}/n_r\mathbb{Z}$ with
$n_1\mid\cdots\mid n_r$.  Each row records one genus-2 curve $C/\mathbb{Q}$ whose
Jacobian is geometrically simple and whose \emph{full} rational torsion group is
\emph{exactly} the displayed group, $\operatorname{Jac}(C)(\mathbb{Q})_{\rm
tors}\cong G$ (the first row is the trivial group, written $[\,]$).  The \emph{Source(s)} column, by
contrast, dates the \emph{group}, and the earliest source there may establish
only a rational point of the given order, a subgroup, or a positive-dimensional
family containing $G$ --- not necessarily an exact-torsion or simple example; the
companion document \texttt{torsion\_sources.tex} classifies every source
individually.  We display, for each group, the equation of the curve of smallest
conductor in the queried database snapshot (2026-07-22) with $\operatorname{Jac}$
geometrically simple and (where such a curve exists; see the endomorphism note)
generic endomorphisms.  Curves of small conductor lie in the production LMFDB;
those found only in the extended database of Booker--Sutherland cite
\cite{BookerSutherland} in the \emph{Source(s)} column.  For convenience each
equation is hyperlinked to its LMFDB home page: production curves by their
permanent label, and curves that are only in the extended database via a
alpha-site (\url{https://alpha.lmfdb.org}) equation-lookup URL, since alpha
labels are not permanent (the lookup matches any isomorphic model of the
curve).  The equations are the intended stable identifiers.

\emph{Source(s) column.} It gives the earliest construction of the torsion,
together with the earliest certification of geometric simplicity when these
differ; ``new'' means the exact simple realization first appears in the present
paper.  A dagger ($\dagger$) marks a source that exhibited the torsion but did
not establish geometric simplicity of a realizing curve; in daggered rows the
certification is supplied by the linked database entry or by the co-cited later
source.  Lepr\'evost's classical order-21, 25 and 27 curves were later shown
non-simple or geometrically split, so those rows cite the first proven-simple
realizations: the certified database curve \texttt{388.a.776.1} realizes $[21]$,
Nicholls \cite{Nicholls2018} realizes $[25]$, and for $[27]$ both Nicholls
\cite{Nicholls2018} and the certified database curve \texttt{604.a.9664.2}
(each distinct from Lepr\'evost's split curve \cite{Leprevost1995}).
Lepr\'evost's $[3,9]$ family is proven
simple over $\mathbb{Q}$ only.  Elkies' page exhibits a rational order-31
\emph{subgroup} with no rational generator (not a realization); the first
\emph{explicit} order-31 realization --- previously the smallest cyclic order
with no explicit realization, and the first on a genus-2 Jacobian --- is the
$[31]$ row: the curve attached to the modular form
\texttt{1830.2.a.q} from the modular-abelian-surfaces dataset \cite{CEHJMPV},
whose exact rational torsion ($\mathbb{Z}/31\mathbb{Z}$, with generator of
irrational Mumford support) and geometric simplicity are certified here.
Order 31 also appears (for $\mathrm{GL}_2$-type abelian varieties, without
curve equations) in the torsion conjectures of Alessandr\`i--Coppola
\cite{AC2026}, whose Theorem 4.2/4.3 proves an order-31 point exists on some
inexplicit member of this isogeny class --- the novelty here is the explicit
curve, its exact group, and its certificates; their $g=2$ lists include 37 as
well, but at the one
admissible level ($2190$, newform \texttt{2190.2.a.v}) no member of the class
is a Jacobian over $\mathbb{Q}$ --- the real multiplication field
$\mathbb{Q}(\sqrt3)$ has narrow class number $2$, so no member carries a
rational principal polarization.

\emph{Endomorphism refinement.} With two exceptions, every group here is realized
by a curve with \emph{generic} geometric endomorphisms,
$\operatorname{End}(\operatorname{Jac}(C)_{\overline{\mathbb{Q}}})=\mathbb{Z}$
(the displayed curve is the smallest-conductor known such witness).  For database
rows this is the LMFDB endomorphism certification (Costa--Mascot--Sijsling--Voight).
The exceptions, marked $^{\rm RM}$, are $[2,22]$ and $[31]$: no generic witness is
known for either.  For $[2,22]$ the extended database's only geometrically simple
curves have real multiplication ($\operatorname{End}_{\overline{\mathbb{Q}}}$ a
real quadratic order), so the displayed curve is geometrically simple but not
generic.  For $[31]$ the only known witness is the displayed curve attached to
the newform \texttt{1830.2.a.q} (real multiplication by $\mathbb{Z}[\sqrt2]$,
conductor $1830^2 = 3348900$, outside the extended database's conductor range);
its $31$-torsion arises from an Eisenstein congruence at the level prime
$61$ ($31 \mid 61+1$).  For
$[2,2,14]$ the extended database likewise contains only RM witnesses; the
displayed curve --- of conductor
$2^2 3^2 5^2 7^2 13^2\cdot 19\cdot 23\cdot 37\cdot 317 = 38200383011700$, the
smallest of ten generic witnesses found here by a three-root refinement of a
contact-7 construction --- is new to this paper.
For each displayed non-database
equation with generic endomorphisms (the rows marked ``new'' without an
$^{\rm RM}$ flag, and the $[2,2,2,10]$ family member) we give a
self-contained certificate in \texttt{data/claude\_endz\_certificates.txt}: a
good prime whose $L$-polynomial is irreducible with all Frobenius-root powers of
degree $4$ (geometric simplicity, by the root-power criterion), a second good
prime whose irreducible $L$-polynomial has splitting field linearly disjoint
from the first (so the centre of the geometric endomorphism algebra is
$\mathbb{Q}$, excluding real and complex multiplication; cf.\
\cite{CLV2021,Zywina2022}), and the observation that $\#\operatorname{Jac}(C)(
\mathbb{Q})_{\rm tors}>18$, which excludes quaternionic multiplication since an
(O-)PQM genus-2 Jacobian over $\mathbb{Q}$ has at most $16$ (resp.\ $18$) rational
torsion points \cite{LSSV2024}.  The RM-flagged non-database row $[31]$ cannot
carry such an $\operatorname{End}=\mathbb{Z}$ certificate (its geometric
endomorphism ring is $\mathbb{Z}[\sqrt2]$); its dedicated certificate --- exact
torsion invariants, the $D_4$/root-power geometric-simplicity witness at
$p=11$, and the real-multiplication identification (constant Frobenius
discriminant core $2$ at all $164$ good primes $p \le 1000$, Hecke-trace match
to \texttt{1830.2.a.q}) --- is in \texttt{data/claude\_z31\_rm\_certificate.txt}.
The one caveat is that the earliest cited simple
realization is occasionally itself non-generic --- $\ddagger$ marks $[5]$, whose
first simple realizations are the RM Jacobian $J_0(31)$ and the CM curve
$y^2-y=x^5$; the displayed witness in every such row is nonetheless generic.
\medskip

\begingroup\small
\renewcommand{\arraystretch}{1.25}
\begin{longtable}{|>{\raggedright\arraybackslash}p{1.7cm}|>{\raggedright\arraybackslash}p{0.70\textwidth}|>{\raggedright\arraybackslash}p{2.1cm}|}
\hline
Group & Curve & Source(s)\\ \hline
\endfirsthead
\hline
Group & Curve & Source(s)\\ \hline
\endhead
"""

body = "".join(f"{g} & {c} & {r}\\\\ \\hline\n" for g, c, r in rows)

footer = r"""\end{longtable}
\endgroup

\begin{thebibliography}{99}
\bibitem{BFT2014} N. Bruin, E.~V. Flynn, and D. Testa, \emph{Descent via (3,3)-isogeny
on Jacobians of genus 2 curves}, Acta Arith.\ \textbf{165} (2014), 201--223.
\bibitem{BLP2009} N. Bernard, F. Lepr\'evost, and M. Pohst, \emph{Jacobians of
genus-2 curves with a rational point of order 11}, Experiment.\ Math.\ \textbf{18}
(2009), 65--70.
\bibitem{BoxallGrant2000} J. Boxall and D. Grant, \emph{Examples of torsion
points on genus two curves}, Trans.\ Amer.\ Math.\ Soc.\ \textbf{352} (2000),
4533--4555.
\bibitem{BSSVY} A.~R. Booker, J. Sijsling, A.~V. Sutherland, J. Voight, and D. Yasaki,
\emph{A database of genus-2 curves over the rational numbers},
LMS J.~Comput.\ Math.\ \textbf{19(A)} (2016), 235--254.
\bibitem{BookerSutherland} A.~R. Booker and A.~V. Sutherland,
\emph{Genus 2 curves of small conductor}, in preparation.
\bibitem{AC2026} J. Alessandr\`i and N. Coppola, \emph{Torsion points on
$\mathrm{GL}_2$-type abelian varieties}, arXiv:2602.21047 (v3, April 2026);
data repository \url{https://github.com/NirvanaC93/Torsion-of-GL2-abelian-varieties}.
\bibitem{CEHJMPV} E. Costa, N.~D. Elkies, S. Hashimoto, R. Jha, K. Martin,
B. Poonen, and J. Voight, \emph{Modular abelian surfaces} (data repository),
\url{https://github.com/edgarcosta/ModularAbelianSurfaces}, 2024.
% TODO(collaborators): confirm the author list of the dataset/paper before
% publication; taken from the repository (commit 687c10d).
\bibitem{Elkies2002} N.~D. Elkies, \emph{Curves of genus 2 over $\mathbb{Q}$ whose
Jacobians are absolutely simple abelian surfaces with torsion points of high order},
\url{https://people.math.harvard.edu/~elkies/g2_tors.html}, 2001--2002, updated 2010.
\bibitem{Elkies2024} N.~D. Elkies, \emph{Families of genus-2 curves with 5-torsion},
LuCaNT: LMFDB, Computation, and Number Theory, Contemp.\ Math.\ \textbf{796},
Amer.\ Math.\ Soc., 2024, 165--185.
\bibitem{Flynn1990} E.~V. Flynn, \emph{Large rational torsion on abelian varieties},
J.~Number Theory \textbf{36} (1990), 257--265.
\bibitem{Flynn1991} E.~V. Flynn, \emph{Sequences of rational torsions on abelian
varieties}, Invent.\ Math.\ \textbf{106} (1991), 433--442.
\bibitem{HLP1996} E.~W. Howe, F. Lepr\'evost, and B. Poonen, \emph{Sous-groupes de
torsion d'ordres \'elev\'es de Jacobiennes d\'ecomposables de courbes de genre 2},
C.~R.\ Acad.\ Sci.\ Paris S\'er.\ I \textbf{323} (1996), 1031--1034.
\bibitem{HLP2000} E.~W. Howe, F. Lepr\'evost, and B. Poonen, \emph{Large torsion
subgroups of split Jacobians of curves of genus two or three}, Forum Math.\
\textbf{12} (2000), 315--364.
\bibitem{Howe2015} E.~W. Howe, \emph{Genus-2 Jacobians with torsion points of
large order}, Bull.\ London Math.\ Soc.\ \textbf{47} (2015), 127--135.
\bibitem{Leprevost1991a} F. Lepr\'evost, \emph{Famille de courbes de genre 2 munies
d'une classe de diviseurs rationnels d'ordre 13},
C.~R.\ Acad.\ Sci.\ Paris S\'er.\ I Math.\ \textbf{313} (1991), 451--454.
\bibitem{Leprevost1991b} F. Lepr\'evost, \emph{Familles de courbes de genre 2 munies
d'une classe de diviseurs rationnels d'ordre 15, 17, 19 ou 21},
C.~R.\ Acad.\ Sci.\ Paris S\'er.\ I Math.\ \textbf{313} (1991), 771--774.
\bibitem{Leprevost1993} F. Lepr\'evost, \emph{Points rationnels de torsion de
jacobiennes de certaines courbes de genre 2}, C.~R.\ Acad.\ Sci.\ Paris S\'er.\ I
\textbf{316} (1993), 819--821.
\bibitem{Leprevost1995} F. Lepr\'evost, \emph{Jacobiennes de certaines courbes de
genre 2: torsion et simplicit\'e}, J.~Th\'eor.\ Nombres Bordeaux \textbf{7} (1995),
283--306.
\bibitem{Leprevost1997} F. Lepr\'evost, \emph{Sur certains sous-groupes de torsion
de jacobiennes de courbes hyperelliptiques de genre $g\ge 1$},
Manuscripta Math.\ \textbf{92} (1997), 47--63.
\bibitem{LPS2004} F. Lepr\'evost, M. Pohst, and A. Sch\"opp, \emph{Rational torsion
of $J_0(N)$ for hyperelliptic modular curves and families of Jacobians of genus 2 and
genus 3 curves with a rational point of order 5, 7 or 10}, Abh.\ Math.\ Sem.\ Univ.\
Hamburg \textbf{74} (2004), 193--203.
\bibitem{Ogg1973} A. Ogg, \emph{Rational points on certain elliptic modular curves},
Proc.\ Sympos.\ Pure Math.\ XXIV, Amer.\ Math.\ Soc., 1973, 221--231.
\bibitem{PP2015} V.~P. Platonov and M.~M. Petrunin, \emph{New curves of genus 2 over
the field of rational numbers whose Jacobians contain torsion points of high order},
Dokl.\ Math.\ \textbf{91} (2015), 220--221.
\bibitem{Nicholls2018} C. Nicholls, \emph{Descent methods and torsion on Jacobians
of higher genus curves}, Ph.D.\ thesis, University of Oxford, 2018.
\bibitem{Ogawa1994} H. Ogawa, \emph{Curves of genus 2 with a rational torsion
divisor of order 23}, Proc.\ Japan Acad.\ Ser.\ A \textbf{70} (1994), 295--298.
\bibitem{PP2012} V.~P. Platonov and M.~M. Petrunin, \emph{New orders of torsion
points in Jacobians of curves of genus 2 over the rational number field},
Dokl.\ Math.\ \textbf{85} (2012), 286--288.
\bibitem{PZP2013} V.~P. Platonov, V.~S. Zhgun, and M.~M. Petrunin,
\emph{On the simplicity of Jacobians for hyperelliptic curves of genus 2 over the
field of rational numbers with torsion points of high order},
Dokl.\ Math.\ \textbf{87} (2013), 318--321.
\bibitem{CLV2021} E. Costa, D. Lombardo, and J. Voight, \emph{Identifying central
endomorphisms of an abelian variety via Frobenius endomorphisms},
Res.\ Number Theory \textbf{7} (2021), no.~46.
\bibitem{Zywina2022} D. Zywina, \emph{Determining monodromy groups of abelian
varieties}, Res.\ Number Theory \textbf{8} (2022), no.~89.
\bibitem{LSSV2024} J. Laga, C. Schembri, A. Shnidman, and J. Voight,
\emph{Rational torsion points on abelian surfaces with quaternionic
multiplication}, Forum Math.\ Sigma \textbf{12} (2024), e92.
\end{thebibliography}
\end{document}
"""

import sys
out = sys.argv[1] if len(sys.argv) > 1 else "paper/torsion_realizations.tex"
with open(out, "w") as fh:
    fh.write(header + body + footer)
print(f"wrote {out}: {len(rows)} rows")
