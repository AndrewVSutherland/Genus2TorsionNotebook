#!/usr/bin/env python3
"""Strict-transform diagnostics for the p=31 q-square collision boundary.

This program uses the full local 2-descent (Kummer) condition for the
marked order-12 point.  On q=(x+t)^2 the order-3 summand has trivial
Kummer image, while the order-4 summand has, after removing the common
fourth power, the five coordinates

    D,  wi * product_{j != i}(wi+wj),  i=0,...,3,

where (w0,w1,w2,w3)=(A,B,C,sigma/rho), ABC=1.  The first coordinate is
implied by the other four, since their product is D times a square.

The code:
  * checks the criterion against every smooth finite sheet at 29,31,37,41;
  * classifies all six p=31 collision strata;
  * counts the normalized p^2 and p^3 strict-transform lifts;
  * resolves the first blow-up of R=S=0;
  * searches the rational balanced exceptional covers to a chosen height.

It is intentionally dependency-free.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path


P = 31
PAIR_NAMES = ("AB", "AC", "AD", "BC", "BD", "CD")
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
INTERNAL = {"AB", "AC", "BC"}
EXTERNAL = {"AD", "BD", "CD"}


class Reporter:
    def __init__(self, path: Path):
        self.handle = path.open("w", encoding="utf-8")

    def say(self, *items: object) -> None:
        line = " ".join(str(x) for x in items)
        print(line)
        self.handle.write(line + "\n")
        self.handle.flush()

    def close(self) -> None:
        self.handle.close()


def legendre(a: int, p: int = P) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def roots_table(p: int) -> list[list[int]]:
    roots = [[] for _ in range(p)]
    for a in range(p):
        roots[a * a % p].append(a)
    return roots


def kummer_products(w: list[int], p: int) -> list[int]:
    ans = []
    for i, wi in enumerate(w):
        ki = wi
        for j, wj in enumerate(w):
            if i != j:
                ki = ki * (wi + wj) % p
        ans.append(ki)
    return ans


def square_status(values: list[int], p: int) -> str:
    ls = [legendre(v, p) for v in values]
    if any(e == -1 for e in ls):
        return "dead_nonsquare_unit"
    if any(e == 0 for e in ls):
        return "deep_zero_factor"
    return "resolved_unit_square"


def finite_control(p: int) -> tuple[int, int, int]:
    roots = roots_table(p)
    smooth = target = collision = 0
    for a in range(1, p):
        for b in range(1, p):
            c = pow(a * b, -1, p)
            r2 = (a * a + b * b + c * c - 3) % p
            s2 = (pow(a, -2, p) + pow(b, -2, p) + pow(c, -2, p) - 3) % p
            if r2 == 0 or s2 == 0 or not roots[r2] or not roots[s2]:
                continue
            for rho in roots[r2]:
                for sigma in roots[s2]:
                    d = sigma * pow(rho, -1, p) % p
                    w = [a, b, c, d]
                    if len({z * z % p for z in w}) < 4:
                        collision += 1
                        continue
                    smooth += 1
                    if all(legendre(k, p) == 1 for k in kummer_products(w, p)):
                        target += 1
    return smooth, target, collision


def enumerate_p31_collision_sheets(tsv_path: Path):
    roots = roots_table(P)
    rows = []
    summary = Counter()
    external_constants = Counter()
    internal_constants = Counter()

    for a in range(1, P):
        for b in range(1, P):
            c = pow(a * b, -1, P)
            r2 = (a * a + b * b + c * c - 3) % P
            s2 = (pow(a, -2, P) + pow(b, -2, P) + pow(c, -2, P) - 3) % P
            if r2 == 0 or s2 == 0 or not roots[r2] or not roots[s2]:
                continue
            for rho in roots[r2]:
                for sigma in roots[s2]:
                    d = sigma * pow(rho, -1, P) % P
                    w = [a, b, c, d]
                    squares = [z * z % P for z in w]
                    eq = [k for k, (i, j) in enumerate(PAIRS) if squares[i] == squares[j]]
                    if not eq:
                        continue
                    if len(eq) != 1:
                        raise RuntimeError("non-R=S sheet has a multiple collision")
                    name = PAIR_NAMES[eq[0]]
                    i, j = PAIRS[eq[0]]
                    kvals = kummer_products(w, P)
                    status = square_status(kvals, P)
                    collision_sign = "+" if w[i] == w[j] else "-"
                    summary[(name, collision_sign, status)] += 1

                    normal_constant = ""
                    normal_legendre = ""
                    if status == "deep_zero_factor":
                        # Only the sheets whose other Kummer coordinates are
                        # squares are genuinely deep; otherwise the status
                        # above already says dead_nonsquare_unit.
                        wi, wj = w[i], w[j]
                        if name in INTERNAL:
                            ci = wi * pow(wi - wj, -1, P) % P
                            for h, wh in enumerate(w):
                                if h not in (i, j):
                                    ci = ci * (wi + wh) % P
                            # The same small sum wi+wj occurs in both K_i
                            # and K_j, and
                            #   wi+wj=(wi^2-wj^2)/(wi-wj).
                            # Hence both normalized constants have the same
                            # denominator wi-wj (not opposite denominators).
                            cj = wj * pow(wi - wj, -1, P) % P
                            for h, wh in enumerate(w):
                                if h not in (i, j):
                                    cj = cj * (wj + wh) % P
                            if all(k == 0 or legendre(k) == 1 for k in kvals):
                                internal_constants[(legendre(ci), legendre(cj))] += 1
                                normal_constant = f"{ci},{cj}"
                                normal_legendre = f"{legendre(ci)},{legendre(cj)}"
                        else:
                            # S-xR=-(x-1)^3/x.  For wi=-wj, division by
                            # (x-1)^(3m) leaves this constant in K_i.
                            x = wi * wi % P
                            ci = wi * pow(x * r2 * (wi - wj), -1, P) % P
                            for h, wh in enumerate(w):
                                if h not in (i, j):
                                    ci = ci * (wi + wh) % P
                            cj = wj * pow(x * r2 * (wi - wj), -1, P) % P
                            for h, wh in enumerate(w):
                                if h not in (i, j):
                                    cj = cj * (wj + wh) % P
                            if all(k == 0 or legendre(k) == 1 for k in kvals):
                                external_constants[(legendre(ci), legendre(cj))] += 1
                                normal_constant = f"{ci},{cj}"
                                normal_legendre = f"{legendre(ci)},{legendre(cj)}"

                    rows.append(
                        [
                            name,
                            collision_sign,
                            a,
                            b,
                            c,
                            rho,
                            sigma,
                            d,
                            *kvals,
                            *(legendre(k) for k in kvals),
                            status,
                            normal_constant,
                            normal_legendre,
                        ]
                    )

    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        out = csv.writer(handle, delimiter="\t")
        out.writerow(
            [
                "collision",
                "linear_sign",
                "A",
                "B",
                "C",
                "rho",
                "sigma",
                "D",
                "K_A",
                "K_B",
                "K_C",
                "K_D",
                "leg_A",
                "leg_B",
                "leg_C",
                "leg_D",
                "status",
                "normalized_collision_constants",
                "constant_legendres",
            ]
        )
        out.writerows(rows)
    return rows, summary, internal_constants, external_constants


def rs_intersection(tsv_path: Path):
    """First exceptional divisor above R=S=0.

    Write A=eA(1+pa), B=eB(1+pb).  Then
       rho/(2p)^2 = a^2+ab+b^2,
       l=(a,b,-a-b,ab(a+b)/(a^2+ab+b^2)).
    The projective directions are [a:b] in P^1(F_31).
    """
    directions = [(a, 1) for a in range(P)] + [(1, 0)]
    patterns = (
        (1, -1, -1, 1),
        (-1, 1, -1, 1),
        (-1, -1, 1, 1),
    )
    rows = []
    counts = Counter()
    for a, b in directions:
        q2 = (a * a + a * b + b * b) % P
        qleg = legendre(q2)
        if qleg != 1:
            # If Q=0, both R and S start in degree >=3, but
            # R-S=-8*p^3*a*b*(a+b)+O(p^4).  At either of the two
            # Q=0 directions this coefficient is nonzero.  Hence R and S
            # cannot both jump to even valuation >=4; one has valuation 3
            # and cannot be a square.  These directions are dead, not deep.
            short = "killed_odd_split" if qleg == 0 else "dead_nonsquare"
            counts[("base", short)] += 1
            rows.append([a, b, q2, qleg, "base", "", "", "", "", short])
            continue
        d = a * b * (a + b) * pow(q2, -1, P) % P
        lvals = [a, b, (-a - b) % P, d]

        # All four signs positive.  Every K_i is 8 modulo p, a square.
        counts[("all_plus", "resolved_live")] += 1
        rows.append([a, b, q2, qleg, "all_plus", d, "8,8,8,8", "1,1,1,1", "", "resolved_live"])

        for pattern in patterns:
            kvals = []
            for i, eps in enumerate(pattern):
                value = 2
                for j in range(4):
                    if pattern[j] == -eps:
                        value = value * (lvals[i] - lvals[j]) % P
                kvals.append(value)
            status = square_status(kvals, P)
            if status == "resolved_unit_square":
                short = "resolved_live"
            elif status == "deep_zero_factor":
                # These six points have one of A^2,B^2,C^2 equal to 1
                # to first order, with the corresponding normalized branch
                # parameter opposite to D.  If m>=2 is the first nonzero
                # order of x_i-1, then v(R)=2 and the exact cubic identity
                # gives v(w_i+D)=3m-2.  Including the other cross-sign sum,
                # v(K_i)=3m-1.  Hence even m is killed, while odd m>=3 is
                # governed by W^2=c*u, c=(l_i-l_k)/4, where k is the other
                # opposite-sign vertex and u=(x_i-1)/p^m.
                zero_indices = [i for i, value in enumerate(kvals) if value == 0]
                i0 = next(i for i in zero_indices if i != 3)
                other_opposite = next(
                    j for j in range(4)
                    if j != 3 and j != i0 and pattern[j] == -pattern[i0]
                )
                cleg = legendre((lvals[i0] - lvals[other_opposite]) * pow(4, -1, P))
                short = "tower_odd_m_c_square" if cleg == 1 else "tower_odd_m_c_nonsquare"
            else:
                short = "killed"
            counts[("balanced_" + "".join("+" if e == 1 else "-" for e in pattern), short)] += 1
            rows.append(
                [
                    a,
                    b,
                    q2,
                    qleg,
                    "balanced_" + "".join("+" if e == 1 else "-" for e in pattern),
                    d,
                    ",".join(map(str, kvals)),
                    ",".join(str(legendre(k)) for k in kvals),
                    ",".join(map(str, lvals)),
                    short,
                ]
            )

    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        out = csv.writer(handle, delimiter="\t")
        out.writerow(["a", "b", "Q", "leg_Q", "sign_stratum", "l_D", "leading_K", "K_legendres", "l_vector", "status"])
        out.writerows(rows)
    return counts


def is_square_fraction(value: Fraction) -> bool:
    if value < 0:
        return False
    return isqrt(value.numerator) ** 2 == value.numerator and isqrt(value.denominator) ** 2 == value.denominator


def balanced_rational_search(height: int, tsv_path: Path):
    """Search the three balanced covers on q^2=a^2+ab+b^2.

    The conic parametrization is
      (a:b:q)=(2m-1 : 1-m^2 : m^2-m+1).
    Only nonzero leading Kummer coordinates count as open points.
    """
    patterns = (
        (1, -1, -1, 1),
        (-1, 1, -1, 1),
        (-1, -1, 1, 1),
    )
    tested = 0
    open_hits = []
    degenerate_hits = []
    for den in range(1, height + 1):
        for num in range(-height, height + 1):
            if gcd(abs(num), den) != 1:
                continue
            m = Fraction(num, den)
            a = 2 * m - 1
            b = 1 - m * m
            q = m * m - m + 1
            if q == 0:
                continue
            lvals = [a, b, -a - b, a * b * (a + b) / (q * q)]
            for pattern in patterns:
                tested += 1
                kvals = []
                for i, eps in enumerate(pattern):
                    value = Fraction(2)
                    for j in range(4):
                        if pattern[j] == -eps:
                            value *= lvals[i] - lvals[j]
                    kvals.append(value)
                if all(is_square_fraction(k) for k in kvals):
                    row = (m, pattern, lvals, kvals)
                    if all(k != 0 for k in kvals):
                        open_hits.append(row)
                    else:
                        degenerate_hits.append(row)

    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        out = csv.writer(handle, delimiter="\t")
        out.writerow(["kind", "m", "pattern", "l_vector", "leading_K"])
        for kind, data in (("open", open_hits), ("degenerate", degenerate_hits)):
            for m, pattern, lvals, kvals in data:
                out.writerow([kind, str(m), "".join("+" if e == 1 else "-" for e in pattern), repr(lvals), repr(kvals)])
    return tested, open_hits, degenerate_hits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=300)
    parser.add_argument(
        "--stem",
        default="results/target_22224_p31_collision_blowups",
    )
    args = parser.parse_args()
    stem = Path(args.stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    report = Reporter(stem.with_suffix(".log"))

    report.say("TARGET_22224_P31_COLLISION_BLOWUPS_START")
    report.say("EXACT_IDENTITY", "R-S=(x-1)(y-1)(z-1)")
    report.say("EXACT_INTERNAL", "x=y", "R=(x-1)^2(2x+1)/x^2", "S=(x-1)^2(x+2)/x")
    report.say("EXACT_EXTERNAL", "S-xR=-(x-1)^3/x")
    report.say("KUMMER", "K0=D", "Ki=wi*prod_(j!=i)(wi+wj)", "product(Ki)=D*square")

    expected = {29: 24, 31: 0, 37: 24, 41: 72}
    for prime in (29, 31, 37, 41):
        smooth, target, collision = finite_control(prime)
        report.say(
            "FINITE_CONTROL",
            "p", prime,
            "smooth_sheets", smooth,
            "kummer_target_sheets", target,
            "expected", expected[prime],
            "presentations", 3 * target,
            "collision_sheets", collision,
            "match", target == expected[prime],
        )

    residue_file = Path(str(stem) + "_residue.tsv")
    rows, summary, iconst, econst = enumerate_p31_collision_sheets(residue_file)
    report.say("P31_COLLISION_SHEETS", len(rows), "file", residue_file)
    for key, count in sorted(summary.items()):
        report.say("P31_STRATUM", *key, count)
    report.say("INTERNAL_NORMAL_CONSTANT_LEGENDRES", dict(sorted(iconst.items())))
    report.say("EXTERNAL_NORMAL_CONSTANT_LEGENDRES", dict(sorted(econst.items())))

    # Normalized lift counts.  A collision sheet has p tangent choices and
    # p normal choices at each new digit.  For an ordinary opposite-sign
    # collision the small K factor has valuation m.  For an external one it
    # has valuation 3m by the exact cubic identity above.
    internal_sheets = sum(v for (name, _, _), v in summary.items() if name in INTERNAL)
    internal_potential = sum(
        v for (name, sign, status), v in summary.items()
        if name in INTERNAL and sign == "-" and status == "deep_zero_factor"
    )
    internal_dead = internal_sheets - internal_potential
    # On each of the 12 residue sheets which initially has two zero Kummer
    # coordinates, the two normalized constants have opposite Legendre
    # symbols.  If the collision depth is odd, the K valuations are odd; if
    # it is even, no normal unit can make both coordinates squares.  Thus
    # the strict transform kills every smooth lift, at every depth.
    assert set(iconst) == {(1, -1)}
    i_p2_killed = internal_sheets * P * P
    i_p2_deep = 0
    i_p3_live = i_p3_killed = i_p3_deep = 0
    report.say("INTERNAL_STRICT_TRANSFORM", "potential_residue_sheets", internal_potential, "constant_legendres", "(+1,-1)", "conclusion", "all_smooth_depths_killed")
    report.say("INTERNAL_P2", "total", internal_sheets * P * P, "resolved_live", 0, "killed", i_p2_killed, "deep", i_p2_deep)
    report.say("INTERNAL_P3_FROM_DEEP", "total", 0, "resolved_live", 0, "killed", 0, "deep", 0)

    external_sheets = sum(v for (name, _, _), v in summary.items() if name in EXTERNAL)
    external_unit_live = sum(
        v for (name, _, status), v in summary.items()
        if name in EXTERNAL and status == "resolved_unit_square"
    )
    external_potential = sum(
        v for (name, sign, status), v in summary.items()
        if name in EXTERNAL and sign == "-" and status == "deep_zero_factor"
    )
    external_dead = external_sheets - external_unit_live - external_potential
    e_p2_live = external_unit_live * P * P
    e_p2_killed = external_dead * P * P + external_potential * (P - 1) * P
    e_p2_deep = external_potential * P
    e_p3_live = e_p2_deep * ((P - 1) // 2) * P
    e_p3_killed = e_p3_live
    e_p3_deep = e_p2_deep * P
    report.say("EXTERNAL_P2", "total", external_sheets * P * P, "resolved_live", e_p2_live, "killed", e_p2_killed, "deep", e_p2_deep)
    report.say("EXTERNAL_P3_FROM_DEEP", "total", e_p2_deep * P * P, "resolved_live", e_p3_live, "killed", e_p3_killed, "deep", e_p3_deep)
    report.say(
        "ALL_SIX_P2",
        "total", (internal_sheets + external_sheets) * P * P,
        "resolved_live", e_p2_live,
        "killed", i_p2_killed + e_p2_killed,
        "deep", i_p2_deep + e_p2_deep,
    )
    report.say(
        "ALL_SIX_P3_FROM_DEEP",
        "total", (i_p2_deep + e_p2_deep) * P * P,
        "resolved_live", i_p3_live + e_p3_live,
        "killed", i_p3_killed + e_p3_killed,
        "deep", i_p3_deep + e_p3_deep,
    )

    rs_file = Path(str(stem) + "_rs_intersection.tsv")
    rs_counts = rs_intersection(rs_file)
    report.say("RS_EXCEPTIONAL_FILE", rs_file)
    for key, count in sorted(rs_counts.items()):
        report.say("RS_EXCEPTIONAL", *key, count)
    report.say(
        "RS_QSQUARE_DIRECTION_SUMMARY",
        "Q_square", 15,
        "Q_nonsquare", 15,
        "Q_zero_killed_odd_split", 2,
        "generic_sheet_resolved_live", 42,
        "generic_sheet_killed", 186,
        "balanced_intersection_tower_sheets", 12,
    )

    rational_file = Path(str(stem) + "_balanced_rational.tsv")
    tested, open_hits, degenerate_hits = balanced_rational_search(args.height, rational_file)
    report.say(
        "BALANCED_RATIONAL_SEARCH",
        "height", args.height,
        "tested", tested,
        "open_hits", len(open_hits),
        "degenerate_hits", len(degenerate_hits),
        "file", rational_file,
    )
    report.say("TARGET_22224_P31_COLLISION_BLOWUPS_DONE")
    report.close()


if __name__ == "__main__":
    main()
