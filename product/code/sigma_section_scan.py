#!/usr/bin/env python3
"""sigma_section_scan.py -- Stage 3 section detection on sigma-surface point lists.

Reads product/data/sigma_si<si>_d<drop>_pts.txt (si in 2..6, drop in 1..3;
columns: t u deck, rationals like -7/3, deck 1=j-equal deck pair, 0=genuine,
-1=unknown), merges the drops per sigma, and hunts for structure:
  1. fiber statistics (by u = C2-fibration, by t = swapped fibration),
  2. exact linear sections t = a*u + b (O(n^2) pair hashing over Q),
  3. deg<=2/deg<=2 rational-function sections t = (a0+a1 u+a2 u^2)/(b0+b1 u+b2 u^2)
     via width-7 sliding windows over sorted distinct u, float SVD nullspace
     prefilter + exact Fraction verification,
  4. si5 ([2,3,1]) vs si6 ([3,1,2]) (u,t)-swap cross-check.
Candidates are printed with exact coefficients for symbolic Magma verification.

Usage (from product/code/):
    python3 sigma_section_scan.py [--datadir ../data] [--si 2,3,4,5,6]
                                  [--report ../logs/sigma_sections_report.txt]
    python3 sigma_section_scan.py --selftest
"""

import argparse
import itertools
import sys
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
except ImportError:  # deg2 detector needs numpy; the rest works without
    np = None

SIGMA_NAMES = {2: "[2,1,3]=(12)", 3: "[3,2,1]=(13)", 4: "[1,3,2]=(23)",
               5: "[2,3,1]=(123)", 6: "[3,1,2]=(132)"}
LIN_MIN_SUPPORT = 5       # points and distinct u
DEG2_MIN_DISTINCT_U = 5
COMBO_CAP = 4096
SV_REL_TOL = 1e-9
CHECKER_REL_TOL = 1e-6


def parse_rat(s):
    return Fraction(s)


def ffloat(x):
    try:
        return float(x)
    except OverflowError:
        return float("inf") if x > 0 else float("-inf")


def load_sigma(si, datadir):
    """Merge the three drop files; returns (points dict {(t,u): deck}, per-file counts)."""
    pts = {}
    counts = {}
    for drop in (1, 2, 3):
        fn = Path(datadir) / f"sigma_si{si}_d{drop}_pts.txt"
        if not fn.exists():
            counts[drop] = None
            continue
        n = 0
        with open(fn) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 3:
                    continue
                t, u, deck = parse_rat(parts[0]), parse_rat(parts[1]), int(parts[2])
                n += 1
                key = (t, u)
                old = pts.get(key)
                if old is None or (old == -1 and deck != -1):
                    pts[key] = deck
                elif old != -1 and deck != -1 and old != deck:
                    # should not happen (deck is a function of (t,u)); be conservative
                    pts[key] = 1
        counts[drop] = n
    return pts, counts


# ---------------------------------------------------------------- statistics

def group_by(points, idx):
    g = {}
    for (t, u) in points:
        key = (t, u)[idx]
        val = (t, u)[1 - idx]
        g.setdefault(key, []).append(val)
    return g


def fmt_list(vals, cap=8):
    vals = sorted(vals, key=lambda x: (abs(x.numerator) + x.denominator, x))
    s = ", ".join(str(v) for v in vals[:cap])
    if len(vals) > cap:
        s += f", ... ({len(vals)} total)"
    return s


def stats_block(out, points):
    n = len(points)
    ndeck = sum(1 for d in points.values() if d == 1)
    ngen = sum(1 for d in points.values() if d == 0)
    nunk = n - ndeck - ngen
    out(f"points: {n} unique (t,u); genuine {ngen}, deck {ndeck}, unknown {nunk}")
    for idx, label, fibvar in ((1, "by-u (C2-fibration fibers)", "u"),
                               (0, "by-t (swapped fibration)", "t")):
        g = group_by(points, idx)
        hist = {}
        for k, v in g.items():
            hist[len(v)] = hist.get(len(v), 0) + 1
        out(f"{label}: {len(g)} distinct {fibvar}; points-per-{fibvar} histogram "
            + str(dict(sorted(hist.items()))))
        top = sorted(g.items(), key=lambda kv: (-len(kv[1]), abs(kv[0].numerator) + kv[0].denominator))[:15]
        for k, v in top:
            if len(v) < 2:
                break
            out(f"    {fibvar}={k}  count={len(v)}  partner-vals: {fmt_list(v)}")
    return ngen


# ------------------------------------------------------------ linear sections

def linear_support(points_list, a, b):
    sup = [(t, u) for (t, u) in points_list if t == a * u + b]
    return sup, len({u for (_, u) in sup})


def linear_scan(points_list):
    """Exact t = a*u + b detection. Returns list of (a, b, support_pts, ndistinct_u)."""
    n = len(points_list)
    cands = {}
    if n <= 4000:
        seen = {}
        for i in range(n):
            t1, u1 = points_list[i]
            for j in range(i + 1, n):
                t2, u2 = points_list[j]
                if u1 == u2:
                    continue
                a = (t1 - t2) / (u1 - u2)
                b = t1 - a * u1
                s = seen.setdefault((a, b), set())
                s.add(i)
                s.add(j)
        for (a, b), s in seen.items():
            if len(s) >= LIN_MIN_SUPPORT:
                sup, ndu = linear_support(points_list, a, b)
                if len(sup) >= LIN_MIN_SUPPORT and ndu >= LIN_MIN_SUPPORT:
                    cands[(a, b)] = (sup, ndu)
    else:
        # float bucket prefilter, exact verification of heavy buckets
        buckets = {}
        tf = [ffloat(t) for (t, _) in points_list]
        uf = [ffloat(u) for (_, u) in points_list]
        for i in range(n):
            for j in range(i + 1, n):
                if uf[i] == uf[j]:
                    continue
                af = (tf[i] - tf[j]) / (uf[i] - uf[j])
                bf = tf[i] - af * uf[i]
                key = (round(af, 7), round(bf, 7))
                cnt, reps = buckets.get(key, (0, []))
                if len(reps) < 3:
                    reps.append((i, j))
                buckets[key] = (cnt + 1, reps)
        need = LIN_MIN_SUPPORT * (LIN_MIN_SUPPORT - 1) // 2 - 2  # slack for split buckets
        for (cnt, reps) in buckets.values():
            if cnt < need:
                continue
            for (i, j) in reps:
                t1, u1 = points_list[i]
                t2, u2 = points_list[j]
                if u1 == u2:
                    continue
                a = (t1 - t2) / (u1 - u2)
                b = t1 - a * u1
                if (a, b) in cands:
                    continue
                sup, ndu = linear_support(points_list, a, b)
                if len(sup) >= LIN_MIN_SUPPORT and ndu >= LIN_MIN_SUPPORT:
                    cands[(a, b)] = (sup, ndu)
    return [(a, b, sup, ndu) for (a, b), (sup, ndu) in sorted(cands.items())]


# ------------------------------------------------- deg<=2/deg<=2 sections

def row_of(t, u):
    return [Fraction(1), u, u * u, -t, -t * u, -t * u * u]


def exact_nullspace(rows):
    """Nullspace basis of an m x 6 Fraction matrix."""
    m = [r[:] for r in rows]
    ncols = 6
    pivots = []
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        inv = m[r][c]
        m[r] = [x / inv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == len(m):
            break
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        v = [Fraction(0)] * ncols
        v[fc] = Fraction(1)
        for i, pc in enumerate(pivots):
            v[pc] = -m[i][fc]
        basis.append(v)
    return basis


def normalize_coeffs(vec):
    """Scale a rational 6-vector to coprime integers, first nonzero positive."""
    from math import gcd
    den = 1
    for x in vec:
        den = den * x.denominator // gcd(den, x.denominator)
    ints = [int(x * den) for x in vec]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    if g:
        ints = [x // g for x in ints]
    lead = next((x for x in ints if x != 0), 0)
    if lead < 0:
        ints = [-x for x in ints]
    return tuple(ints)


def deg2_support(coeffs, points_list):
    a0, a1, a2, b0, b1, b2 = [Fraction(c) for c in coeffs]
    sup = []
    for (t, u) in points_list:
        den = b0 + b1 * u + b2 * u * u
        # den == 0 makes the cross-multiplied test vacuous when num and den
        # share the root u (the rational function is undefined there); such
        # points must not count as support
        if den == 0:
            continue
        if a0 + a1 * u + a2 * u * u == t * den:
            sup.append((t, u))
    return sup, len({u for (_, u) in sup})


def deg2_scan(points_list, all_points_list, out):
    """Windowed deg<=2/deg<=2 detector. points_list drives windows (genuine),
    all_points_list is used for the support count. Returns candidate dict."""
    if np is None:
        out("deg2 sections: numpy unavailable -- SKIPPED")
        return {}
    by_u = {}
    for (t, u) in points_list:
        by_u.setdefault(u, []).append(t)
    us = sorted(by_u)
    if len(us) < 7:
        out(f"deg2 sections: only {len(us)} distinct genuine u -- window scan skipped")
        return {}
    cands = {}
    nwin = nskip = nfloat = 0
    for w in range(len(us) - 6):
        base = us[w:w + 6]
        checker = us[w + 6]
        combos = 1
        for u in base:
            combos *= len(by_u[u])
        if combos > COMBO_CAP:
            nskip += 1
            continue
        nwin += 1
        tchoices = [by_u[u] for u in base]
        mats = []
        picks = []
        for combo in itertools.product(*tchoices):
            rowsf = []
            for t, u in zip(combo, base):
                tf, uf = ffloat(t), ffloat(u)
                row = [1.0, uf, uf * uf, -tf, -tf * uf, -tf * uf * uf]
                mx = max(abs(x) for x in row)
                rowsf.append([x / mx for x in row] if mx > 0 else row)
            mats.append(rowsf)
            picks.append(combo)
        A = np.array(mats)
        try:
            _, S, Vt = np.linalg.svd(A)
        except np.linalg.LinAlgError:
            continue
        for k in range(A.shape[0]):
            if S[k][0] <= 0 or S[k][-1] / S[k][0] >= SV_REL_TOL:
                continue
            v = Vt[k][-1]
            a0, a1, a2, b0, b1, b2 = v
            ucf = ffloat(checker)
            den = b0 + b1 * ucf + b2 * ucf * ucf
            if abs(den) < 1e-12:
                continue
            tpred = (a0 + a1 * ucf + a2 * ucf * ucf) / den
            ok = any(abs(tpred - ffloat(tc)) <= CHECKER_REL_TOL * max(1.0, abs(ffloat(tc)))
                     for tc in by_u[checker])
            if not ok:
                continue
            nfloat += 1
            rows = [row_of(t, u) for t, u in zip(picks[k], base)]
            for basisvec in exact_nullspace(rows):
                if all(x == 0 for x in basisvec[3:6]):
                    continue  # denominator identically zero
                key = normalize_coeffs(basisvec)
                if key in cands:
                    continue
                sup, ndu = deg2_support(key, all_points_list)
                supg, ndug = deg2_support(key, points_list)
                if ndug >= DEG2_MIN_DISTINCT_U:
                    cands[key] = (len(supg), ndug, len(sup), ndu)
    out(f"deg2 sections: {nwin} windows scanned, {nskip} skipped (combo cap), "
        f"{nfloat} float-passers, {len(cands)} exact candidates")
    return cands


# ------------------------------------------------------------------ reports

def analyze_sigma(si, pts, counts, out):
    out("")
    out(f"=== sigma si={si} {SIGMA_NAMES.get(si, '?')} ===")
    cstr = ", ".join(f"d{d}: " + ("missing" if counts[d] is None else f"{counts[d]} lines")
                     for d in (1, 2, 3))
    out(f"files: {cstr}")
    if not pts:
        out("no points -- skipped")
        return []
    stats_block(out, pts)
    genuine = sorted([p for p, d in pts.items() if d == 0])
    allpts = sorted(pts)
    found = []

    lin_gen = linear_scan(genuine)
    if lin_gen:
        for (a, b, sup, ndu) in lin_gen:
            out(f"LINEAR SECTION (genuine): t = ({a})*u + ({b})  support {len(sup)} pts / {ndu} distinct u")
            out(f"    pts: {', '.join(f'({t},{u})' for t, u in sup[:6])}" +
                (" ..." if len(sup) > 6 else ""))
            found.append(("linear", si, (a, b), len(sup), ndu, "genuine"))
    else:
        out("linear sections (genuine): NONE")
    lin_all = linear_scan(allpts)
    extra = [c for c in lin_all if (c[0], c[1]) not in {(a, b) for (a, b, _, _) in lin_gen}]
    if extra:
        for (a, b, sup, ndu) in extra:
            ndeck = sum(1 for p in sup if pts[p] == 1)
            out(f"LINEAR SECTION (incl deck/unknown): t = ({a})*u + ({b})  "
                f"support {len(sup)} pts ({ndeck} deck) / {ndu} distinct u  [DECK-DRIVEN?]")
            found.append(("linear", si, (a, b), len(sup), ndu, "with-deck"))
    else:
        out("linear sections incl deck: no additional candidates")

    d2 = deg2_scan(genuine, allpts, out)
    for key, (nsg, ndug, nsa, ndua) in sorted(d2.items()):
        a0, a1, a2, b0, b1, b2 = key
        out(f"DEG2 SECTION: t = ({a0} + {a1}*u + {a2}*u^2)/({b0} + {b1}*u + {b2}*u^2)  "
            f"support genuine {nsg} pts/{ndug} u, all {nsa} pts/{ndua} u")
        found.append(("deg2", si, key, nsa, ndua, "windowed"))
    return found


def swap_check(pts5, pts6, out):
    out("")
    out("=== si5/si6 (u,t)-swap cross-check ===")
    if not pts5 or not pts6:
        out("skipped (missing data)")
        return
    s5 = set(pts5)
    s6sw = {(u, t) for (t, u) in pts6}
    only5 = sorted(s5 - s6sw)[:10]
    only6 = sorted(s6sw - s5)[:10]
    if not only5 and not only6:
        out(f"PASS: si5 point set ({len(s5)}) is exactly the (u,t)-swap of si6 ({len(pts6)})")
    else:
        out(f"FAIL: {len(s5 - s6sw)} points only in si5, {len(s6sw - s5)} only in swapped si6")
        for p in only5:
            out(f"    only-si5: {p}")
        for p in only6:
            out(f"    only-si6-swapped: {p}")
    out("note: a linear section t=a*u+b of si5 appears in si6 as u=a*t+b (swapped fibration)")


# ----------------------------------------------------------------- selftest

def selftest():
    ok = True
    pts = {}
    # planted deg2 section t = (3u^2-1)/(u+2) at u=1..12
    for k in range(1, 13):
        u = Fraction(k)
        t = (3 * u * u - 1) / (u + 2)
        pts[(t, u)] = 0
    # planted linear section t = 2u-3 at u=1..6
    for k in range(1, 7):
        u = Fraction(k)
        pts[(2 * u - 3, u)] = 0
    # noise at planted u's
    for k, tv in [(1, 17), (2, -5), (3, Fraction(7, 2)), (4, 23), (5, -9),
                  (6, Fraction(1, 3)), (7, 40), (8, -21)]:
        pts[(Fraction(tv), Fraction(k))] = 0
    # noise away from planted u's
    for k, tv in [(-5, 2), (-3, 7), (Fraction(-1, 2), 4), (15, 1), (16, -3),
                  (17, 8), (18, Fraction(2, 7)), (19, 5), (21, -2), (22, 9),
                  (23, 11), (25, -6), (26, 13), (27, Fraction(3, 4))]:
        pts[(Fraction(tv), Fraction(k))] = 0
    assert len(pts) == 40, f"selftest fabricated {len(pts)} points, wanted 40"

    lines = []
    out = lines.append
    genuine = sorted(pts)
    lin = linear_scan(genuine)
    want_lin = (Fraction(2), Fraction(-3))
    got_lin = any((a, b) == want_lin and len(sup) >= 5 for (a, b, sup, ndu) in lin)
    d2 = deg2_scan(genuine, genuine, out)
    want_d2 = normalize_coeffs([Fraction(-1), Fraction(0), Fraction(3),
                                Fraction(2), Fraction(1), Fraction(0)])
    got_d2 = want_d2 in d2 and d2[want_d2][1] >= 10
    print("\n".join(lines))
    print(f"selftest linear detector: {'found' if got_lin else 'MISSED'} t=2u-3 "
          f"(candidates: {[(str(a), str(b)) for a, b, _, _ in lin]})")
    if np is None:
        # the deg2 detector is an explicitly optional (numpy-only) mode:
        # its absence must not fail the selftest of the exact-arithmetic core
        print("selftest deg2 detector:   SKIPPED (numpy unavailable)")
        ok = got_lin
    else:
        print(f"selftest deg2 detector:   {'found' if got_d2 else 'MISSED'} t=(3u^2-1)/(u+2) "
              f"(normalized {want_d2}; candidates: {list(d2)})")
        ok = got_lin and got_d2
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# --------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default="../data")
    ap.add_argument("--si", default="2,3,4,5,6")
    ap.add_argument("--report", default="../logs/sigma_sections_report.txt")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(selftest())

    silist = [int(s) for s in args.si.split(",")]
    lines = []

    def out(s):
        lines.append(s)
        print(s)

    out(f"sigma_section_scan: datadir={args.datadir} si={silist}")
    all_found = []
    loaded = {}
    for si in silist:
        pts, counts = load_sigma(si, args.datadir)
        loaded[si] = pts
        all_found += analyze_sigma(si, pts, counts, out)
    if 5 in loaded and 6 in loaded:
        swap_check(loaded[5], loaded[6], out)

    out("")
    out("=== SUMMARY ===")
    if not all_found:
        out("candidate sections: NONE")
    else:
        for kind, si, coeffs, nsup, ndu, note in all_found:
            if kind == "linear":
                a, b = coeffs
                out(f"si={si} LINEAR t = ({a})*u + ({b})  [{nsup} pts, {ndu} u, {note}]")
            else:
                a0, a1, a2, b0, b1, b2 = coeffs
                out(f"si={si} DEG2 t = ({a0} + {a1}*u + {a2}*u^2)/({b0} + {b1}*u + {b2}*u^2)"
                    f"  [{nsup} pts, {ndu} u, {note}]")
    rp = Path(args.report)
    rp.parent.mkdir(parents=True, exist_ok=True)
    rp.write_text("\n".join(lines) + "\n")
    print(f"[report written to {rp}]")


if __name__ == "__main__":
    main()
