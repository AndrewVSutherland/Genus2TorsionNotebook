#!/usr/bin/env sage --python
"""Finite-field diagnostic for contact-5 order-40 plus extra 2-torsion.

The order-40 cover is Q(Y,s)=0 from contact5_order40_cover_analysis.py,
where rational points require s+1 and Y to be squares.  Since

    s = r^2 - 1,    t = (2-s)/s,

an independent rational 2-torsion point would also force the residual quartic
f/(x-1) to be reducible modulo every good prime away from boundary.

This script enumerates the finite-field order-40 cover and records how often
the residual quartic is reducible.  If some prime has no such open point, the
[2,40] route is forced to the boundary at that prime.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

from sage.all import GF, PolynomialRing, QQ, prime_range


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "order40_cover", HERE / "contact5_order40_cover_analysis.py"
)
order40_cover = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
# contact5_order40_cover_analysis.py only runs main under __name__ == "__main__".
SPEC.loader.exec_module(order40_cover)


def eval_poly_mod(poly, yv, sv, F):
    total = F(0)
    for mon, coeff in poly.dict().items():
        total += F(coeff) * (yv ** mon[0]) * (sv ** mon[1])
    return total


def residual_quartic(t, F):
    P = PolynomialRing(F, "x")
    x = P.gen()
    inv2 = F(2) ** -1
    inv4 = F(4) ** -1
    b = (t * t - 1) * inv2
    h = 1 + t * x + b * x**2
    f = h**2 - (t + 1) ** 4 * inv4 * x**5
    if f.degree() != 5 or f.discriminant() == 0:
        return None
    if f(1) != 0:
        return None
    q, rem = f.quo_rem(x - 1)
    if rem != 0 or q.degree() != 4 or q.discriminant() == 0:
        return None
    return q


def factor_type(poly):
    return "+".join(str(f.degree()) for f, _ in sorted(poly.factor(), key=lambda fe: fe[0].degree()))


def analyze_prime(poly, p, max_samples=12):
    F = GF(p)
    squares = {a * a for a in F}
    stats = {
        "cover": 0,
        "open": 0,
        "extra2": 0,
        "fullsplit": 0,
        "bad": 0,
    }
    factor_counts = {}
    s_extra = set()
    samples = []

    for sv in F:
        if sv in {F(0), F(-1)}:
            continue
        if sv + 1 not in squares:
            continue
        t = (F(2) - sv) / sv
        for yv in F:
            if yv == 0 or yv not in squares:
                continue
            if eval_poly_mod(poly, yv, sv, F) != 0:
                continue
            stats["cover"] += 1
            q = residual_quartic(t, F)
            if q is None:
                stats["bad"] += 1
                continue
            stats["open"] += 1
            ftype = factor_type(q)
            factor_counts[ftype] = factor_counts.get(ftype, 0) + 1
            if ftype != "4":
                stats["extra2"] += 1
                s_extra.add(int(sv))
                if ftype == "1+1+1+1":
                    stats["fullsplit"] += 1
                if len(samples) < max_samples:
                    samples.append((int(sv), int(yv), int(t), ftype))

    return stats, factor_counts, sorted(s_extra), samples


def main() -> None:
    prime_bound = 101
    if len(sys.argv) >= 2:
        prime_bound = int(sys.argv[1])
    primes = [p for p in prime_range(3, prime_bound + 1) if p not in {5}]

    poly = order40_cover.quotient_polynomial()
    print("Contact-5 order40 + extra2 finite intersection diagnostic")
    print("prime_bound", prime_bound, "primes", len(primes))
    print("columns: p cover open extra2 fullsplit bad factor_counts s_extra samples")

    obstructing = []
    for p in primes:
        stats, factor_counts, s_extra, samples = analyze_prime(poly, int(p))
        if stats["extra2"] == 0:
            obstructing.append(int(p))
        print(
            "p", int(p),
            "cover", stats["cover"],
            "open", stats["open"],
            "extra2", stats["extra2"],
            "fullsplit", stats["fullsplit"],
            "bad", stats["bad"],
            "factor_counts", sorted(factor_counts.items()),
            "s_extra", s_extra[:20],
            "samples", samples,
        )

    print("DONE")
    print("obstructing_primes", obstructing)


if __name__ == "__main__":
    main()
