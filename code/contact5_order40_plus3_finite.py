#!/usr/bin/env sage --python
"""Finite-field diagnostic for contact-5 order-40 plus rational 3-torsion.

The order-40 cover is Q(Y,s)=0 from contact5_order40_cover_analysis.py.
Rational order-40 points require

    s = r^2 - 1,    Y = D^2,

so modulo a good prime we enumerate open cover points with s+1 and Y
nonzero squares.  A rational 3-torsion class would force

    3 | #Jac(C)(F_p)

at every good prime p != 3.  This script records whether the open order-40
cover has any residues satisfying that necessary condition.
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
SPEC.loader.exec_module(order40_cover)


def eval_poly_mod(poly, yv, sv, F):
    total = F(0)
    for mon, coeff in poly.dict().items():
        total += F(coeff) * (yv ** mon[0]) * (sv ** mon[1])
    return total


def family_polynomial_finite(t, F):
    P = PolynomialRing(F, "x")
    x = P.gen()
    inv2 = F(2) ** -1
    inv4 = F(4) ** -1
    b = (t * t - 1) * inv2
    h = 1 + t * x + b * x**2
    f = h**2 - (t + 1) ** 4 * inv4 * x**5
    return f


def is_good_family_parameter(t, F):
    f = family_polynomial_finite(t, F)
    return f.degree() == 5 and f.discriminant() != 0


def quadratic_character(a):
    if a == 0:
        return 0
    return 1 if a.is_square() else -1


def curve_points_odd_degree(f, F):
    # One rational point at infinity for odd degree.
    return 1 + sum(1 + quadratic_character(f(a)) for a in F)


def jacobian_order_mod3(t, F):
    p = F.characteristic()
    f = family_polynomial_finite(t, F)
    if f.degree() != 5 or f.discriminant() == 0:
        return None

    n1 = curve_points_odd_degree(f, F)

    F2 = GF(p**2, "u")
    P2 = PolynomialRing(F2, "x")
    x2 = P2.gen()
    coeffs = [F2(c) for c in f.list()]
    f2 = sum(coeffs[i] * x2**i for i in range(len(coeffs)))
    n2 = curve_points_odd_degree(f2, F2)

    c1 = n1 - p - 1
    c2 = (n2 - p**2 - 1 + c1**2) // 2
    jac_order = 1 + c1 + c2 + p * c1 + p**2
    return jac_order % 3


def analyze_prime(poly, p, max_samples=12):
    F = GF(p)
    squares = {a * a for a in F}
    stats = {
        "cover": 0,
        "open": 0,
        "plus3": 0,
        "bad": 0,
    }
    t_cache = {}
    s_plus3 = set()
    samples = []

    for sv in F:
        if sv in {F(0), F(-1)}:
            continue
        if sv + 1 not in squares:
            continue
        t = (F(2) - sv) / sv
        t_key = int(t)
        for yv in F:
            if yv == 0 or yv not in squares:
                continue
            if eval_poly_mod(poly, yv, sv, F) != 0:
                continue

            stats["cover"] += 1
            if not is_good_family_parameter(t, F):
                stats["bad"] += 1
                continue

            stats["open"] += 1
            if t_key not in t_cache:
                t_cache[t_key] = jacobian_order_mod3(t, F)
            if t_cache[t_key] == 0:
                stats["plus3"] += 1
                s_plus3.add(int(sv))
                if len(samples) < max_samples:
                    samples.append((int(sv), int(yv), int(t)))

    return stats, sorted(s_plus3), samples


def main() -> None:
    prime_bound = int(sys.argv[1]) if len(sys.argv) >= 2 else 101
    primes = [int(p) for p in prime_range(7, prime_bound + 1) if p != 3]

    poly = order40_cover.quotient_polynomial()
    print("Contact-5 order40 + 3 finite intersection diagnostic")
    print("prime_bound", prime_bound, "primes", len(primes))
    print("columns: p cover open plus3 bad s_plus3 samples")

    obstructing = []
    for p in primes:
        stats, s_plus3, samples = analyze_prime(poly, p)
        if stats["plus3"] == 0:
            obstructing.append(p)
        print(
            "p", p,
            "cover", stats["cover"],
            "open", stats["open"],
            "plus3", stats["plus3"],
            "bad", stats["bad"],
            "s_plus3", s_plus3[:30],
            "samples", samples,
            flush=True,
        )

    print("DONE")
    print("obstructing_primes", obstructing)


if __name__ == "__main__":
    main()
