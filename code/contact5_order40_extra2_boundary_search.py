#!/usr/bin/env sage --python
"""Boundary-filtered rational search for contact-5 [2,40] candidates.

Finite-field diagnostics show that the open order-40 cover has no extra
2-torsion intersection at p=7,11,13.  Therefore a rational [2,40] candidate
must reduce to the boundary at each of these primes.

For the order-40 cover Q(Y,s)=0 with s=r^2-1, this script keeps only r whose
s-reduction is boundary at p=7,11,13.  Boundary means s=0, s=-1, s=infinity,
or Q(0,s)=0 so that the cover point can have Y=0 modulo p.  It then computes
rational Y-roots of Q(Y,r^2-1), keeps square Y, and finally tests whether the
residual quartic f/(x-1) is reducible over Q.
"""

from __future__ import annotations

import importlib.util
from math import gcd as py_gcd
from pathlib import Path
import sys

from sage.all import GF, PolynomialRing, QQ, ZZ


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "order40_cover", HERE / "contact5_order40_cover_analysis.py"
)
order40_cover = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(order40_cover)


OBSTRUCTION_PRIMES = [7, 11, 13]


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
            # Y=infinity boundary, detected by vanishing leading Y coefficient.
            lv = F(0)
            for mon, coeff in leading_y.dict().items():
                lv += F(coeff) * (sv ** mon[0])
            if lv == 0:
                vals.add(sv)
        out[p] = sorted(int(v) for v in vals)
    return out


def s_boundary_at_p(r, p, allowed_s):
    num = ZZ(r.numerator())
    den = ZZ(r.denominator())
    # r has pole, hence s=r^2-1 has pole.
    if den % p == 0:
        return True
    F = GF(p)
    rv = F(num) / F(den)
    sv = rv * rv - F(1)
    return int(sv) in allowed_s[p]


def family_polynomial(t):
    P = PolynomialRing(QQ, "x")
    x = P.gen()
    b = (t * t - 1) / 2
    h = 1 + t * x + b * x**2
    f = h**2 - ((t + 1) ** 4 / 4) * x**5
    return P(f), P(h), b


def residual_factor_type(t):
    P = PolynomialRing(QQ, "x")
    x = P.gen()
    f, h, b = family_polynomial(t)
    if f.degree() != 5 or f.discriminant() == 0 or f(1) != 0:
        return "bad", [], b
    q, rem = f.quo_rem(x - 1)
    if rem != 0 or q.degree() != 4:
        return "bad", [], b
    fac = q.factor()
    degs = sorted([ff.degree() for ff, ee in fac for _ in range(ee)])
    return "+".join(str(d) for d in degs), fac, b


def main() -> None:
    height = int(sys.argv[1]) if len(sys.argv) >= 2 else 1000
    out_path = Path(sys.argv[2]) if len(sys.argv) >= 3 else None
    progress = 200000

    poly = order40_cover.quotient_polynomial()
    Y, s = poly.parent().gens()
    Uy = PolynomialRing(QQ, "Y")
    allowed_s = boundary_s_sets(poly)

    lines = []
    checked = 0
    boundary_survivors = 0
    fibers_with_roots = 0
    square_cover = 0
    extra2 = 0
    factor_counts = {}

    print("Contact-5 order40+extra2 boundary-filtered search")
    print("height", height)
    print("boundary_s", allowed_s)

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
                "extra2", extra2,
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
            ftype, fac, b = residual_factor_type(tval)
            factor_counts[ftype] = factor_counts.get(ftype, 0) + 1
            if ftype != "4" and ftype != "bad":
                extra2 += 1
                print(
                    "HIT", "r", r, "s", sval, "Y", yval, "D", dval,
                    "t", tval, "b", b, "factor_type", ftype,
                    "mult", mult,
                    flush=True,
                )
            lines.append(
                " ".join(
                    [
                        qstr(r), qstr(sval), qstr(yval), qstr(dval),
                        qstr(tval), ftype, str(mult),
                    ]
                )
            )

    header = (
        "# contact5 order40+extra2 boundary-filtered search\n"
        f"# height={height}\n"
        f"# boundary_s={allowed_s}\n"
        f"# checked={checked} boundary_survivors={boundary_survivors} "
        f"fibers_with_roots={fibers_with_roots} square_cover={square_cover} "
        f"extra2={extra2} factor_counts={sorted(factor_counts.items())}\n"
        "# columns: r s Y D t residual_factor_type mult\n"
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
    print("extra2", extra2)
    print("factor_counts", sorted(factor_counts.items()))


if __name__ == "__main__":
    main()
