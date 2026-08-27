#!/usr/bin/env python3
"""5-adic valuation scan for the contact-7 first-halving surface.

This is a diagnostic for the remaining [28] -> [56] route.  It uses the same
root-eliminated first-halving equation as contact7_halving_surface_roots.py,
but records the 5-adic type of each tested (s,z,eps) and each rational
first-halving candidate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from math import isqrt
from pathlib import Path
import importlib.util

import sympy as sp


HERE = Path(__file__).resolve().parent
SPEC_ENUM = importlib.util.spec_from_file_location(
    "contact7_enum", HERE / "contact7_halving_surface_enum.py"
)
contact7_enum = importlib.util.module_from_spec(SPEC_ENUM)
assert SPEC_ENUM.loader is not None
SPEC_ENUM.loader.exec_module(contact7_enum)


T = sp.Symbol("T")


def qstr(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def to_sympy(q: Fraction) -> sp.Rational:
    return sp.Rational(q.numerator, q.denominator)


def from_sympy(q: sp.Rational) -> Fraction:
    q = sp.Rational(q)
    return Fraction(int(q.p), int(q.q))


def rational_roots_in_u2(c3: Fraction, c2: Fraction, c1: Fraction, z: Fraction) -> list[Fraction]:
    c3s = to_sympy(c3)
    c2s = to_sympy(c2)
    c1s = to_sympy(c1)
    zs = to_sympy(z)
    a = c3s**2 + 8 * zs - 4 * c2s
    coeffs = [
        sp.QQ.one,
        -4 * c3s,
        4 * c3s**2 + 2 * a - 64 * zs,
        -4 * c3s * a + 64 * zs * c3s - 64 * c1s,
        a**2,
    ]
    poly = sp.Poly.from_list(coeffs, gens=T, domain=sp.QQ)
    roots = poly.ground_roots()
    return sorted({from_sympy(root) for root in roots})


def is_square_fraction(q: Fraction) -> Fraction | None:
    if q < 0:
        return None
    n = isqrt(q.numerator)
    d = isqrt(q.denominator)
    if n * n == q.numerator and d * d == q.denominator:
        return Fraction(n, d)
    return None


def vp_int(n: int, p: int = 5) -> int:
    if n == 0:
        return 10**9
    n = abs(n)
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def vp(q: Fraction, p: int = 5) -> int:
    if q == 0:
        return 10**9
    return vp_int(q.numerator, p) - vp_int(q.denominator, p)


def bucket(v: int) -> str:
    if v == 10**9:
        return "zero"
    if v <= -2:
        return "<=-2"
    if v == -1:
        return "-1"
    if v == 0:
        return "0"
    if v == 1:
        return "1"
    return ">=2"


def residue(q: Fraction, p: int = 5) -> str:
    if q.denominator % p == 0:
        return "inf"
    return str((q.numerator % p) * pow(q.denominator % p, -1, p) % p)


def input_signature(s: Fraction, z: Fraction, eps: int) -> str:
    s_eps = Fraction(eps) * s
    r = 1 - s * s
    fields = [
        f"vs={bucket(vp(s))}",
        f"vz={bucket(vp(z))}",
        f"v_r={bucket(vp(r))}",
        f"v_epss_plus1={bucket(vp(s_eps + 1))}",
        f"sres={residue(s)}",
        f"zres={residue(z)}",
        f"eps={eps}",
    ]
    return ",".join(fields)


def candidate_signature(s: Fraction, u: Fraction, z: Fraction, eps: int, a: Fraction, b: Fraction, r: Fraction) -> str:
    s_eps = Fraction(eps) * s
    fields = [
        f"vs={bucket(vp(s))}",
        f"vu={bucket(vp(u))}",
        f"vz={bucket(vp(z))}",
        f"va={bucket(vp(a))}",
        f"vb={bucket(vp(b))}",
        f"v_r={bucket(vp(r))}",
        f"v_epss_plus1={bucket(vp(s_eps + 1))}",
        f"sres={residue(s)}",
        f"ures={residue(u)}",
        f"zres={residue(z)}",
        f"eps={eps}",
    ]
    return ",".join(fields)


def boundary_reasons(s: Fraction, u: Fraction, z: Fraction, eps: int, a: Fraction, b: Fraction, r: Fraction) -> tuple[str, ...]:
    s_eps = Fraction(eps) * s
    reasons: list[str] = []
    if vp(s) < 0:
        reasons.append("s_nonintegral")
    if vp(u) < 0:
        reasons.append("u_nonintegral")
    if vp(z) < 0:
        reasons.append("z_nonintegral")
    if vp(a) < 0:
        reasons.append("a_nonintegral")
    if vp(b) < 0:
        reasons.append("b_nonintegral")
    if vp(s) > 0:
        reasons.append("s_zero_mod5")
    if vp(u) > 0:
        reasons.append("u_zero_mod5")
    if vp(z) > 0:
        reasons.append("z_zero_mod5")
    if vp(r) > 0:
        reasons.append("r_zero_mod5")
    if vp(s_eps + 1) > 0:
        reasons.append("epss_plus1_zero_mod5")
    return tuple(reasons) if reasons else ("open_integral",)


def height(q: Fraction) -> int:
    return max(abs(q.numerator), q.denominator)


def run(args: argparse.Namespace) -> None:
    s_vals = contact7_enum.params(args.s_height)
    z_vals = contact7_enum.params(args.z_height)

    checked = quartics = root_t = zero_t = square_t = hits = boundary_hits = 0
    input_counts: Counter[str] = Counter()
    input_root_counts: Counter[str] = Counter()
    input_square_counts: Counter[str] = Counter()
    cand_counts: Counter[str] = Counter()
    reason_counts: Counter[tuple[str, ...]] = Counter()
    rows: list[str] = []

    for s in s_vals:
        if s == 0 or s * s == 1:
            continue
        for z in z_vals:
            for eps in (-1, 1):
                sig_in = input_signature(s, z, eps)
                input_counts[sig_in] += 1
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
                except Exception as exc:
                    print("ROOT_FAIL", qstr(s), qstr(z), eps, repr(exc), flush=True)
                    continue

                for t0 in roots:
                    root_t += 1
                    input_root_counts[sig_in] += 1
                    if t0 == 0:
                        zero_t += 1
                        continue
                    u0 = is_square_fraction(t0)
                    if u0 is None:
                        continue
                    square_t += 1
                    input_square_counts[sig_in] += 1

                    for u in sorted({u0, -u0}):
                        if args.u_height and height(u) > args.u_height:
                            continue
                        v = (u * u - c3) / 2
                        w = (v * v + 2 * z - c2) / (2 * u)
                        if c1 != w * w - 2 * v * z:
                            continue
                        ar = contact7_enum.root_a_r(s, b, eps)
                        if ar is None:
                            continue
                        a, r = ar
                        reasons = boundary_reasons(s, u, z, eps, a, b, r)
                        boundary = reasons != ("open_integral",)
                        if args.only_boundary and not boundary:
                            continue
                        hits += 1
                        if boundary:
                            boundary_hits += 1
                        sig_c = candidate_signature(s, u, z, eps, a, b, r)
                        cand_counts[sig_c] += 1
                        reason_counts[reasons] += 1
                        if len(rows) < args.max_rows:
                            rows.append(
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
                                        "|".join(reasons),
                                    ]
                                )
                            )

    lines: list[str] = []
    lines.append("# contact7 first-halving 5-adic valuation scan")
    lines.append(f"# s_height={args.s_height} z_height={args.z_height} u_height={args.u_height}")
    lines.append(f"# only_boundary={args.only_boundary}")
    lines.append(
        "# "
        + " ".join(
            [
                f"checked={checked}",
                f"quartics={quartics}",
                f"root_t={root_t}",
                f"zero_t={zero_t}",
                f"square_t={square_t}",
                f"hits={hits}",
                f"boundary_hits={boundary_hits}",
            ]
        )
    )
    lines.append("# columns: s u z eps a b r v w boundary_reasons")
    lines.extend(rows)
    lines.append("")
    lines.append("REASON_COUNTS")
    for key, count in sorted(reason_counts.items(), key=lambda item: (str(item[0]), item[1])):
        lines.append(f"{'|'.join(key)} {count}")
    lines.append("")
    lines.append("CANDIDATE_SIGNATURE_COUNTS")
    for key, count in sorted(cand_counts.items()):
        lines.append(f"{key} {count}")
    lines.append("")
    lines.append("INPUT_SIGNATURES_WITH_SQUARE_ROOTS")
    for key, count in sorted(input_square_counts.items()):
        lines.append(f"{key} square_roots={count} roots={input_root_counts[key]} checked={input_counts[key]}")
    lines.append("")
    lines.append("TOP_INPUT_SIGNATURES")
    for key, count in input_counts.most_common(args.top_inputs):
        lines.append(f"{key} checked={count} roots={input_root_counts[key]} square_roots={input_square_counts[key]}")
    text = "\n".join(lines) + "\n"
    if args.out:
        Path(args.out).write_text(text, encoding="ascii")
    else:
        print(text, end="")

    print(
        "DONE",
        f"checked={checked}",
        f"quartics={quartics}",
        f"root_t={root_t}",
        f"zero_t={zero_t}",
        f"square_t={square_t}",
        f"hits={hits}",
        f"boundary_hits={boundary_hits}",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s-height", type=int, default=12)
    parser.add_argument("--z-height", type=int, default=40)
    parser.add_argument("--u-height", type=int, default=0)
    parser.add_argument("--only-boundary", action="store_true")
    parser.add_argument("--max-rows", type=int, default=50)
    parser.add_argument("--top-inputs", type=int, default=40)
    parser.add_argument("--progress", type=int, default=100000)
    parser.add_argument("--out", default="")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
