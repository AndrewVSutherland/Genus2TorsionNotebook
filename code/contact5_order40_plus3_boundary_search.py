#!/usr/bin/env sage --python
"""Boundary-filtered rational search for contact-5 order-40 plus 3-torsion.

The finite diagnostic contact5_order40_plus3_finite.py shows that the open
order-40 cover has no residues with 3 | #Jac(F_p) at p = 7, 11, 17.
Therefore a rational order-40 plus rational 3-torsion candidate must reduce
to a boundary/bad point at each of these primes.

This script filters r by those forced boundary residues, solves the rational
fibers of Q(Y, r^2 - 1), keeps square Y, and then applies the good-prime
necessary condition 3 | #Jac(F_p) at the remaining small primes.
"""

from __future__ import annotations

import importlib.util
from math import gcd as py_gcd
from pathlib import Path
import sys

from sage.all import GF, PolynomialRing, QQ, ZZ, prime_range


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "order40_cover", HERE / "contact5_order40_cover_analysis.py"
)
order40_cover = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(order40_cover)


OBSTRUCTION_PRIMES = [7, 11, 17]


def qstr(q):
    q = QQ(q)
    if q.denominator() == 1:
        return str(q.numerator())
    return f"{q.numerator()}/{q.denominator()}"


def is_square_q(q):
    q = QQ(q)
    if q < 0:
        return False
    return ZZ(q.numerator()).is_square() and ZZ(q.denominator()).is_square()


def sqrt_q(q):
    q = QQ(q)
    return QQ(ZZ(q.numerator()).sqrt()) / QQ(ZZ(q.denominator()).sqrt())


def rational_parameters(bound):
    for den in range(1, bound + 1):
        for num in range(-bound, bound + 1):
            if py_gcd(num, den) == 1:
                yield QQ(num) / QQ(den)


def eval_poly_mod(poly, yv, sv, F):
    total = F(0)
    for mon, coeff in poly.dict().items():
        total += F(coeff) * (yv ** mon[0]) * (sv ** mon[1])
    return total


def family_polynomial_q(t):
    P = PolynomialRing(QQ, "x")
    x = P.gen()
    b = (t * t - 1) / 2
    h = 1 + t * x + b * x**2
    f = h**2 - ((t + 1) ** 4 / 4) * x**5
    return P(f), b


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


def boundary_s_sets(poly):
    Y, s = poly.parent().gens()
    leading_y = poly.coefficient({Y: poly.degree(Y)})
    out = {}
    for p in OBSTRUCTION_PRIMES:
        F = GF(p)
        vals = {F(0), F(-1)}
        for sv in F:
            # Y=0 boundary.
            if eval_poly_mod(poly, F(0), sv, F) == 0:
                vals.add(sv)

            # Y=infinity boundary, detected by the leading Y coefficient.
            lv = F(0)
            for mon, coeff in leading_y.dict().items():
                lv += F(coeff) * (sv ** mon[0])
            if lv == 0:
                vals.add(sv)

            # Bad reduction of the family at this finite s-value.
            if sv != 0:
                t = (F(2) - sv) / sv
                if not is_good_family_parameter(t, F):
                    vals.add(sv)
        out[p] = sorted(int(v) for v in vals)
    return out


def s_boundary_at_p(r, p, allowed_s):
    num = ZZ(r.numerator())
    den = ZZ(r.denominator())
    if den % p == 0:
        return True
    F = GF(p)
    rv = F(num) / F(den)
    sv = rv * rv - F(1)
    return int(sv) in allowed_s[p]


def residue_of_rational(q, p):
    num = ZZ(q.numerator())
    den = ZZ(q.denominator())
    if den % p == 0:
        return False, None
    F = GF(p)
    return True, F(num) / F(den)


def passes_plus3_good_prime_filter(t, primes, cache):
    bad_primes = []
    good_checked = 0
    for p in primes:
        ok, tv = residue_of_rational(t, p)
        if not ok:
            bad_primes.append(p)
            continue
        F = GF(p)
        f = family_polynomial_finite(tv, F)
        if f.degree() != 5 or f.discriminant() == 0:
            bad_primes.append(p)
            continue
        good_checked += 1
        key = (p, int(tv))
        if key not in cache:
            cache[key] = jacobian_order_mod3(tv, F)
        if cache[key] != 0:
            return False, good_checked, bad_primes, p
    return True, good_checked, bad_primes, None


def main() -> None:
    height = int(sys.argv[1]) if len(sys.argv) >= 2 else 1000
    prime_bound = int(sys.argv[2]) if len(sys.argv) >= 3 else 101
    out_path = Path(sys.argv[3]) if len(sys.argv) >= 4 else None
    progress = 200000

    filter_primes = [
        int(p) for p in prime_range(7, prime_bound + 1)
        if int(p) not in {3} and int(p) not in OBSTRUCTION_PRIMES
    ]

    poly = order40_cover.quotient_polynomial()
    Y, s = poly.parent().gens()
    Uy = PolynomialRing(QQ, "Y")
    allowed_s = boundary_s_sets(poly)
    cache = {}

    lines = []
    checked = 0
    boundary_survivors = 0
    fibers_with_roots = 0
    square_cover = 0
    plus3_survivors = 0
    first_kill_counts = {}

    print("Contact-5 order40+3 boundary-filtered search")
    print("height", height, "prime_bound", prime_bound)
    print("obstruction_primes", OBSTRUCTION_PRIMES)
    print("boundary_s", allowed_s)
    print("filter_primes", filter_primes)

    for r in rational_parameters(height):
        if r in {QQ(-1), QQ(0), QQ(1)}:
            continue
        checked += 1
        if progress and checked % progress == 0:
            print(
                "progress", checked,
                "boundary_survivors", boundary_survivors,
                "fibers_with_roots", fibers_with_roots,
                "square_cover", square_cover,
                "plus3_survivors", plus3_survivors,
                flush=True,
            )

        if not all(s_boundary_at_p(r, p, allowed_s) for p in OBSTRUCTION_PRIMES):
            continue
        boundary_survivors += 1

        sval = r * r - 1
        if sval in {QQ(0), QQ(-1)}:
            continue
        spec = Uy(poly.subs({s: sval}))
        roots = spec.roots(QQ)
        if not roots:
            continue
        fibers_with_roots += 1

        for yval, mult in roots:
            if not is_square_q(yval):
                continue
            square_cover += 1
            dval = sqrt_q(yval)
            tval = (2 - sval) / sval
            ok, good_checked, bad_primes, first_kill = passes_plus3_good_prime_filter(
                tval, filter_primes, cache
            )
            if not ok:
                first_kill_counts[first_kill] = first_kill_counts.get(first_kill, 0) + 1
                continue

            plus3_survivors += 1
            f, b = family_polynomial_q(tval)
            print(
                "SURVIVOR", "r", r, "s", sval, "Y", yval, "D", dval,
                "t", tval, "b", b, "good_checked", good_checked,
                "bad_primes", bad_primes, "mult", mult,
                flush=True,
            )
            lines.append(
                " ".join(
                    [
                        qstr(r), qstr(sval), qstr(yval), qstr(dval),
                        qstr(tval), qstr(b), str(good_checked),
                        ",".join(str(p) for p in bad_primes), str(mult),
                    ]
                )
            )

    header = (
        "# contact5 order40+3 boundary-filtered search\n"
        f"# height={height} prime_bound={prime_bound}\n"
        f"# obstruction_primes={OBSTRUCTION_PRIMES}\n"
        f"# boundary_s={allowed_s}\n"
        f"# checked={checked} boundary_survivors={boundary_survivors} "
        f"fibers_with_roots={fibers_with_roots} square_cover={square_cover} "
        f"plus3_survivors={plus3_survivors} "
        f"first_kill_counts={sorted(first_kill_counts.items())}\n"
        "# columns: r s Y D t b good_checked bad_primes mult\n"
    )
    text = header + "\n".join(lines) + ("\n" if lines else "")
    if out_path is not None:
        out_path.write_text(text, encoding="ascii")
    else:
        print(text, end="")

    print("DONE height", height)
    print("checked", checked)
    print("boundary_survivors", boundary_survivors)
    print("fibers_with_roots", fibers_with_roots)
    print("square_cover", square_cover)
    print("plus3_survivors", plus3_survivors)
    print("first_kill_counts", sorted(first_kill_counts.items()))


if __name__ == "__main__":
    main()
