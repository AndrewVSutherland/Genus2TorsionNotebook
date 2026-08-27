#!/usr/bin/env python3
"""Root-based enumerator for the contact-7 first-halving surface.

The older enumerator scans triples (s,u,z).  For fixed (s,z,eps), the
remaining surface equation is a quartic in T = u^2:

    ((T-c3)^2 + 8*z - 4*c2)^2 - 64*T*(z*(T-c3) + c1) = 0.

This script solves that quartic over Q, keeps rational square roots T=u^2,
and writes candidate rows for exact Magma verification.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
import importlib.util

import sympy as sp


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "contact7_enum", HERE / "contact7_halving_surface_enum.py"
)
contact7_enum = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(contact7_enum)


T = sp.Symbol("T")


def qstr(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def to_sympy(q: Fraction) -> sp.Rational:
    return sp.Rational(q.numerator, q.denominator)


def from_sympy(q: sp.Rational) -> Fraction:
    q = sp.Rational(q)
    return Fraction(int(q.p), int(q.q))


def is_square_fraction(q: Fraction) -> Fraction | None:
    if q < 0:
        return None
    n = isqrt(q.numerator)
    d = isqrt(q.denominator)
    if n * n == q.numerator and d * d == q.denominator:
        return Fraction(n, d)
    return None


def reduce_mod_p(q: Fraction, p: int) -> int | None:
    if q.denominator % p == 0:
        return None
    return (q.numerator % p) * pow(q.denominator % p, -1, p) % p


def is_boundary_at_p(s: Fraction, u: Fraction, z: Fraction, eps: int, p: int) -> bool:
    rs = reduce_mod_p(s, p)
    ru = reduce_mod_p(u, p)
    rz = reduce_mod_p(z, p)
    if rs is None or ru is None or rz is None:
        return True
    if rs == 0 or (rs * rs - 1) % p == 0:
        return True
    if ru == 0:
        return True
    S = (eps * rs) % p
    if S == 0 or (S + 1) % p == 0:
        return True
    return False


def rational_roots_in_u2(c3: Fraction, c2: Fraction, c1: Fraction, z: Fraction) -> list[Fraction]:
    C3 = to_sympy(c3)
    C2 = to_sympy(c2)
    C1 = to_sympy(c1)
    Z = to_sympy(z)
    A = C3**2 + 8 * Z - 4 * C2
    coeffs = [
        sp.QQ.one,
        -4 * C3,
        4 * C3**2 + 2 * A - 64 * Z,
        -4 * C3 * A + 64 * Z * C3 - 64 * C1,
        A**2,
    ]
    poly = sp.Poly.from_list(coeffs, gens=T, domain=sp.QQ)
    roots = poly.ground_roots()
    return sorted({from_sympy(root) for root in roots})


def height(q: Fraction) -> int:
    return max(abs(q.numerator), q.denominator)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--s-height", type=int, default=12)
    ap.add_argument("--z-height", type=int, default=40)
    ap.add_argument("--u-height", type=int, default=0, help="optional max height for u")
    ap.add_argument("--boundary-prime", type=int, default=5)
    ap.add_argument("--require-boundary", action="store_true")
    ap.add_argument("--out", default="")
    ap.add_argument("--progress", type=int, default=10000)
    args = ap.parse_args()

    s_vals = contact7_enum.params(args.s_height)
    z_vals = contact7_enum.params(args.z_height)

    checked = 0
    quartics = 0
    root_t = 0
    square_t = 0
    hits = 0
    boundary_hits = 0
    lines: list[str] = []

    for s in s_vals:
        if s == 0 or s * s == 1:
            continue
        for z in z_vals:
            for eps in (-1, 1):
                checked += 1
                if args.progress and checked % args.progress == 0:
                    print(
                        "progress",
                        f"checked={checked}",
                        f"quartics={quartics}",
                        f"root_t={root_t}",
                        f"square_t={square_t}",
                        f"hits={hits}",
                        flush=True,
                    )

                b = contact7_enum.b_from_s_z(s, z, eps)
                if b is None:
                    continue
                c3, c2, c1, c0 = contact7_enum.coeffs_from_s_b(s, b, eps)
                if c0 != z * z:
                    continue
                quartics += 1

                try:
                    roots = rational_roots_in_u2(c3, c2, c1, z)
                except Exception as exc:  # keep long searches moving
                    print(
                        "ROOT_FAIL",
                        "s",
                        qstr(s),
                        "z",
                        qstr(z),
                        "eps",
                        eps,
                        "error",
                        repr(exc),
                        flush=True,
                    )
                    continue

                for t0 in roots:
                    root_t += 1
                    u0 = is_square_fraction(t0)
                    if u0 is None or u0 == 0:
                        continue
                    square_t += 1
                    for u in sorted({u0, -u0}):
                        if args.u_height and height(u) > args.u_height:
                            continue
                        c3, c2, c1, c0 = contact7_enum.coeffs_from_s_b(s, b, eps)
                        v = (u * u - c3) / 2
                        w = (v * v + 2 * z - c2) / (2 * u)
                        if c1 != w * w - 2 * v * z:
                            continue
                        ar = contact7_enum.root_a_r(s, b, eps)
                        if ar is None:
                            continue
                        a, r = ar
                        boundary = is_boundary_at_p(s, u, z, eps, args.boundary_prime)
                        if args.require_boundary and not boundary:
                            continue
                        hits += 1
                        if boundary:
                            boundary_hits += 1
                        lines.append(
                            " ".join(
                                [
                                    qstr(s),
                                    qstr(u),
                                    qstr(z),
                                    str(eps),
                                    qstr(a),
                                    qstr(b),
                                    qstr(r),
                                    qstr(v),
                                    qstr(w),
                                    "boundary" if boundary else "open",
                                ]
                            )
                        )

    header = (
        "# contact7 halving surface root candidates\n"
        f"# s_height={args.s_height} z_height={args.z_height} u_height={args.u_height}\n"
        f"# boundary_prime={args.boundary_prime} require_boundary={args.require_boundary}\n"
        f"# checked={checked} quartics={quartics} root_t={root_t} square_t={square_t} "
        f"hits={hits} boundary_hits={boundary_hits}\n"
        "# columns: s u z eps a b r v w p_boundary\n"
    )
    text = header + "\n".join(lines) + ("\n" if lines else "")
    if args.out:
        with open(args.out, "w", encoding="ascii") as fh:
            fh.write(text)
    else:
        print(text, end="")
    print(
        "DONE",
        f"checked={checked}",
        f"quartics={quartics}",
        f"root_t={root_t}",
        f"square_t={square_t}",
        f"hits={hits}",
        f"boundary_hits={boundary_hits}",
    )


if __name__ == "__main__":
    main()
