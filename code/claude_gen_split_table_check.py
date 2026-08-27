#!/usr/bin/env python3
# SUPERSEDED 2026-08-24: paper/scripts_and_data/verify_split_torsion_table.m
# is now a STATIC data-driven script that reads paper/scripts_and_data/table2.txt
# (the machine-readable table data from which make_tables.py typesets the
# paper table).  Do NOT regenerate it from here -- this generator predates
# the data files and would clobber the data-driven version.  Kept for the
# record of how the original embedded-data scripts were produced.

"""Generate code/verify_split_torsion_table.m — a Magma script that verifies
the exact torsion subgroup of every genus-2 curve in
paper/split_torsion_table.tex, via

    CheckTorsion(I, fh):
      assert I eq Invariants(AbelianGroup(TorsionSubgroup(Jacobian(
          SimplifiedModel(HyperellipticCurve(R!fh[1], R!fh[2]))))));

Inputs: product/data/split_min_witnesses.csv (database rows) plus the nine
non-database witnesses embedded below (identical to those in
code/claude_gen_split_table.py / product/data/new_split_witnesses.txt).
Rows are ordered by coefficient size so the cheap checks run first.

Run from the repository root:  python3 code/claude_gen_split_table_check.py
Then:  cd paper/scripts_and_data && magma -b verify_split_torsion_table.m
"""
import csv, re, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "product", "data", "split_min_witnesses.csv")
OUT = os.path.join(ROOT, "paper", "scripts_and_data", "verify_split_torsion_table.m")

def expand_quintic(roots):
    # x * prod(x + r) -> ascending coefficient list
    poly = [0, 1]  # x
    for r in roots:
        # multiply by (x + r)
        new = [0] * (len(poly) + 1)
        for i, c in enumerate(poly):
            new[i] += c * r
            new[i + 1] += c
        poly = new
    return poly

# non-database rows: read from the SAME structured source as the table
# generator (product/data/split_table_witnesses.json) so witness equations
# cannot drift between the published table and this verifier.
import json
WJSON = os.path.join(ROOT, "product", "data", "split_table_witnesses.json")
NEWROWS = []
for w in json.load(open(WJSON))["witnesses"]:
    if "quintic_roots" in w:
        f, h = expand_quintic(w["quintic_roots"]), []
    else:
        f, h = w["f"], w["h"]
    NEWROWS.append((w["invs"], f, h, w["tag"]))

# ---- drift guard: the shipped table must equal, LINE FOR LINE, the table
# rebuilt right now from the shared sources (CSV + witness JSON + the table
# generator's own rendering).  This catches every drift class at once —
# missing groups, edited witness equations, changed links or citations —
# whenever either generator is rerun without regenerating the other artifact.
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import claude_gen_split_table as _tablegen

TEX = os.path.join(ROOT, "paper", "split_torsion_table.tex")
shipped = []
with open(TEX) as fh:
    body = False
    for ln in fh:
        ln = ln.rstrip("\n")
        if ln.startswith(r"\endhead"):
            body = True
            continue
        if ln.startswith(r"\end{longtable}"):
            body = False
        if body and ln.strip():
            shipped.append(ln)
expected = _tablegen.build_lines()
if shipped != expected:
    import difflib
    diff = "\n".join(difflib.unified_diff(shipped, expected,
                                          "shipped tex", "rebuilt from sources",
                                          lineterm="", n=0))
    raise SystemExit("paper/split_torsion_table.tex is out of sync with the "
                     "witness sources — regenerate it first.\n" + diff[:2000])
table_groups = {tuple(int(t) for t in
                      re.match(r"\$\[([0-9,\\ ]*)\]\$", ln).group(1)
                      .replace("\\,", "").split(",") if t)
                for ln in shipped}

rows = []
with open(CSV) as fh:
    for r in csv.DictReader(fh):
        invs = [int(t) for t in r["torsion_subgroup"].strip("[]").split(",") if t]
        m = re.match(r"\[\[(.*?)\],\[(.*?)\]\]$", r["eqn"].replace(" ", ""))
        f = [int(t) for t in m.group(1).split(",")] if m.group(1) else []
        h = [int(t) for t in m.group(2).split(",")] if m.group(2) else []
        # tag with the production label when the table links to production
        m2 = re.match(r"(\d+)\.([a-z]+)\d+$", r["label"])
        cls = "%s.%s" % (m2.group(1), m2.group(2))
        tag = _tablegen.PROD.get((cls, r["abs_disc"]), r["label"])
        rows.append((invs, f, h, tag))
for invs, f, h, tag in NEWROWS:
    rows.append((invs, f, h, tag))

check_groups = {tuple(r[0]) for r in rows}
if check_groups != table_groups:
    missing = sorted(table_groups - check_groups)
    extra = sorted(check_groups - table_groups)
    raise SystemExit("group mismatch vs paper/split_torsion_table.tex: "
                     "missing witnesses for %s; extra %s" % (missing, extra))

# cheap first: sort by size of the largest coefficient
rows.sort(key=lambda r: max(abs(c) for c in r[1] + (r[2] or [0])))

def mseq(lst):
    return "[" + ",".join(str(c) for c in lst) + "]" if lst else "[0]"

def iseq(lst):   # invariant lists: the empty group is [], not [0]
    return "[" + ",".join(str(c) for c in lst) + "]"

lines = [
    "// verify_split_torsion_table.m -- verifies the exact torsion subgroup of",
    "// every curve in Table 2 (tab:splitcensus) of the paper",
    "// (= paper/split_torsion_table.tex).",
    "// Generated by code/claude_gen_split_table_check.py; run from this directory:",
    "//   magma -b verify_split_torsion_table.m",
    "SetColumns(0);",
    "SetMemoryLimit(8*10^9);",
    "R<x> := PolynomialRing(Rationals());",
    "",
    "procedure CheckTorsion(I, fh)",
    "    assert I eq Invariants(AbelianGroup(TorsionSubgroup(Jacobian(",
    "        SimplifiedModel(HyperellipticCurve(R!fh[1], R!fh[2]))))));",
    "end procedure;",
    "",
    "t0 := Cputime();",
    "n := 0;",
]
for invs, f, h, tag in rows:
    lines.append('printf "checking %-14o  (' + tag.replace('"', "")
                 + ')\\n", ' + iseq(invs) + ";")
    lines.append("CheckTorsion(" + iseq(invs) + ", [" + mseq(f) + ", "
                 + mseq(h) + "]);")
    lines.append("n +:= 1;")
lines += [
    'printf "ALL %o CHECKS PASSED (%.1o s)\\n", n, Cputime()-t0;',
    "quit;",
]

with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("wrote %s (%d checks)" % (OUT, len(rows)))


# ---------------------------------------------------------------------------
# Also emit verify_split_certificates.m: geometric splitness certificates for
# ALL table rows.  Route per row (established by the 2026-08-23 exploration
# run; see the computation-audit note in notes/):
#   R1  = degenerate Richelot codomain over Q (depth 1)
#   Q2  = degenerate Richelot over Q(sqrt 2), depth <= 2
#   SIG = no (2,2)-route exists (odd gluing); verify the reducible-chi split
#         signature at every good p < 200 and cite the source
# ---------------------------------------------------------------------------
OUTC = os.path.join(ROOT, "paper", "scripts_and_data",
                    "verify_split_certificates.m")
SIGTAGS = set()
Q2TAGS = {"256.a.512.1", "256.a2"}
COVERS = json.load(open(os.path.join(ROOT, "product", "data",
                                     "split_cover_certificates.json")))

certrows = []
coverdefs = []
for invs, f, h, tag in rows:
    if tag in SIGTAGS:
        route = "SIG"; extra = ', 0'
    elif tag in Q2TAGS:
        route = "Q2"; extra = ', 0'
    elif tag in COVERS:
        c = COVERS[tag]
        k = len(coverdefs) + 1
        if c["type"] == "Q":
            coverdefs.append(
                'covers[%d] := < "Q", "%s", %d, PxQ![%s], PxQ![%s], PxQ![%s], PxQ!0, PxQ!0 >;'
                % (k, c["E"], c["n"], c["p"], c["q"], c["h"]))
        else:
            coverdefs.append(
                'K := NumberField(%s);\nPK := PolynomialRing(K);\ncovers[%d] := < "K", "%s", %d, PK![ K | %s ], PK![ K | %s ], PK![ K | %s ], PK!0, PK!0 >;\ncoverK[%d] := < K, K!(%s), K!(%s) >;'
                % (c["Kpoly"], k, "", c["n"], c["p"], c["q"], c["h"], k, c["g2"], c["g3"]))
        route = "COVER"; extra = ', %d' % k
    else:
        route = "R1"; extra = ', 0'
    certrows.append('  <%s, %s, %s, "%s", "%s"%s>'
                    % (iseq(invs), mseq(f), mseq(h), route, tag.replace('"', ""), extra))
certdata = ("rows := [\n" + ",\n".join(certrows) + "\n];\n"
            + "covers := AssociativeArray();\ncoverK := AssociativeArray();\n"
            + "\n".join(coverdefs) + "\n")

certtpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
               "claude_split_certificates_template.m")).read()
with open(OUTC, "w") as fh:
    fh.write(certtpl.replace("//<<DATA>>", certdata))
print("wrote %s (%d rows)" % (OUTC, len(certrows)))
