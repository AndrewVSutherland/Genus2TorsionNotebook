#!/usr/bin/env python3
"""Normalized 13-adic boundary analysis for M(2,2,2,8) plus 3-torsion.

We analyze the local system

    K(a,b,c,d) = s2(a,b,c,d)^2 - 4abcd = 0

together with the eliminated triple-contact equations in (L,U,v), where
L=1/m is forced to be a 13-adic unit by the first contact equation.

The good-reduction triple-contact cover has no points over good K3 points
modulo 13.  This script therefore allows bad cover reduction modulo 13:
v may be 0, q may have repeated reduction, and q may meet f modulo 13.  It
then lifts the resulting boundary solutions to 13^2 by linearization and
records which boundary coordinates can acquire a nonzero first coefficient.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

P = 13
PAIR_LABELS = {
    (0, 1): "E12",
    (0, 2): "E13",
    (0, 3): "E14",
    (1, 2): "E23",
    (1, 3): "E24",
    (2, 3): "E34",
}


def inv(x: int, p: int = P) -> int:
    return pow(x % p, -1, p)


def k3_value(z: tuple[int, int, int, int, int, int, int]) -> int:
    a, b, c, d, _L, _U, _v = z
    s2 = a * b + a * c + a * d + b * c + b * d + c * d
    return s2 * s2 - 4 * a * b * c * d


def k3_mod(a: int, b: int, c: int, d: int, p: int = P) -> int:
    s2 = a * b + a * c + a * d + b * c + b * d + c * d
    return (s2 * s2 - 4 * a * b * c * d) % p


def coeffs_mod(a: int, b: int, c: int, d: int, p: int = P) -> tuple[int, int, int, int]:
    aa, bb, cc, dd = (a * a) % p, (b * b) % p, (c * c) % p, (d * d) % p
    return (
        (aa + bb + cc + dd) % p,
        (aa * bb + aa * cc + aa * dd + bb * cc + bb * dd + cc * dd) % p,
        (aa * bb * cc + aa * bb * dd + aa * cc * dd + bb * cc * dd) % p,
        (aa * bb * cc * dd) % p,
    )


def coeffs_int(z: tuple[int, int, int, int, int, int, int]) -> tuple[int, int, int, int]:
    a, b, c, d, _L, _U, _v = z
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    return (
        aa + bb + cc + dd,
        aa * bb + aa * cc + aa * dd + bb * cc + bb * dd + cc * dd,
        aa * bb * cc + aa * bb * dd + aa * cc * dd + bb * cc * dd,
        aa * bb * cc * dd,
    )


def contact_mod_from_coeffs(
    coeffs: tuple[int, int, int, int], L: int, U: int, v: int, p: int = P
) -> tuple[int, int, int]:
    e1, e2, e3, e4 = coeffs
    m2 = (L * L) % p
    a_aux = (2 * m2 * e1 + 6 * (U * U + v * v) - (m2 + 3 * U) ** 2) % p
    return (
        ((m2 + 3 * U) * a_aux + 8 * v**3 - 4 * e2 * m2 - 4 * U**3 - 24 * U * v * v)
        % p,
        (
            a_aux * a_aux
            + 16 * (m2 + 3 * U) * v**3
            - 16 * e3 * m2
            - 48 * (U * U * v * v + v**4)
        )
        % p,
        (a_aux * v**3 - 2 * e4 * m2 - 6 * U * v**4) % p,
    )


def equations_int(z: tuple[int, int, int, int, int, int, int]) -> list[int]:
    a, b, c, d, L, U, v = z
    e1, e2, e3, e4 = coeffs_int(z)
    m2 = L * L
    a_aux = 2 * m2 * e1 + 6 * (U * U + v * v) - (m2 + 3 * U) ** 2
    return [
        k3_value(z),
        (m2 + 3 * U) * a_aux + 8 * v**3 - 4 * e2 * m2 - 4 * U**3 - 24 * U * v * v,
        a_aux * a_aux
        + 16 * (m2 + 3 * U) * v**3
        - 16 * e3 * m2
        - 48 * (U * U * v * v + v**4),
        a_aux * v**3 - 2 * e4 * m2 - 6 * U * v**4,
    ]


def boundary_signature(z: tuple[int, int, int, int, int, int, int], p: int = P) -> tuple[str, ...]:
    xs = [z[i] % p for i in range(4)]
    squares = [(x * x) % p for x in xs]
    labels: list[str] = []
    for i, x in enumerate(xs):
        if x == 0:
            labels.append(f"Z{i + 1}")
    for i in range(4):
        for j in range(i + 1, 4):
            if squares[i] == squares[j]:
                labels.append(PAIR_LABELS[(i, j)])
    return tuple(labels)


def cover_boundary_flags(
    z: tuple[int, int, int, int, int, int, int], p: int = P
) -> tuple[str, ...]:
    a, b, c, d, _L, U, v = z
    e1, e2, e3, e4 = coeffs_mod(a, b, c, d, p)
    flags: list[str] = []
    if v % p == 0:
        flags.append("v=0")
    if (U * U - 4 * v * v) % p == 0:
        flags.append("discq=0")
    if e4 % p == 0:
        flags.append("e4=0")

    for x in range(p):
        qx = (x * x + U * x + v * v) % p
        fx = (x**5 + e1 * x**4 + e2 * x**3 + e3 * x * x + e4 * x) % p
        if qx == 0 and fx == 0:
            flags.append("gcd(q,f)")
            break
    return tuple(flags)


def mod13_solutions(p: int = P) -> list[tuple[int, int, int, int, int, int, int]]:
    sols: list[tuple[int, int, int, int, int, int, int]] = []
    for a, b, c, d in product(range(p), repeat=4):
        if k3_mod(a, b, c, d, p) != 0:
            continue
        coeffs = coeffs_mod(a, b, c, d, p)
        for L in range(1, p):
            for U in range(p):
                for v in range(p):
                    if contact_mod_from_coeffs(coeffs, L, U, v, p) == (0, 0, 0):
                        sols.append((a, b, c, d, L, U, v))
    return sols


def affine_lift_space(
    z: tuple[int, int, int, int, int, int, int], p: int = P
) -> tuple[tuple[int, ...], list[tuple[int, ...]]] | None:
    f0 = equations_int(z)
    rhs = [(-(f // p)) % p for f in f0]
    matrix = [[0] * 7 for _ in range(4)]

    for j in range(7):
        zp = list(z)
        zp[j] += p
        fj = equations_int(tuple(zp))
        for i in range(4):
            matrix[i][j] = ((fj[i] - f0[i]) // p) % p

    aug = [matrix[i][:] + [rhs[i]] for i in range(4)]
    row = 0
    pivots: list[int] = []
    for col in range(7):
        pivot = next((r for r in range(row, 4) if aug[r][col] % p), None)
        if pivot is None:
            continue
        aug[row], aug[pivot] = aug[pivot], aug[row]
        scale = inv(aug[row][col], p)
        aug[row] = [(x * scale) % p for x in aug[row]]
        for r in range(4):
            if r != row and aug[r][col] % p:
                factor = aug[r][col] % p
                aug[r] = [(aug[r][c] - factor * aug[row][c]) % p for c in range(8)]
        pivots.append(col)
        row += 1
        if row == 4:
            break

    for r in range(row, 4):
        if all(aug[r][c] % p == 0 for c in range(7)) and aug[r][7] % p:
            return None

    free = [c for c in range(7) if c not in pivots]
    y0 = [0] * 7
    basis: list[tuple[int, ...]] = []
    for r, col in enumerate(pivots):
        y0[col] = aug[r][7] % p
    for fc in free:
        vec = [0] * 7
        vec[fc] = 1
        for r, col in enumerate(pivots):
            vec[col] = (-aug[r][fc]) % p
        basis.append(tuple(vec))
    return tuple(y0), basis


def coord_form(index: int) -> tuple[int, ...]:
    coeffs = [0] * 7
    coeffs[index] = 1
    return tuple(coeffs)


def disc_form(U: int, v: int, p: int = P) -> tuple[int, ...]:
    coeffs = [0] * 7
    coeffs[5] = (2 * U) % p
    coeffs[6] = (-8 * v) % p
    return tuple(coeffs)


def form_identically_zero(
    y0: tuple[int, ...], basis: list[tuple[int, ...]], coeffs: tuple[int, ...], p: int = P
) -> bool:
    def dot(y: tuple[int, ...]) -> int:
        return sum(c * yy for c, yy in zip(coeffs, y)) % p

    return dot(y0) == 0 and all(dot(vec) == 0 for vec in basis)


def forms_can_be_simultaneously_nonzero(
    y0: tuple[int, ...], basis: list[tuple[int, ...]], forms: list[tuple[int, ...]], p: int = P
) -> bool:
    # Since len(forms) < 13, proper affine hyperplanes cannot cover the
    # affine lift space over F_13.  So it is enough to check that no
    # requested linear form is identically zero on the lift space.
    return all(not form_identically_zero(y0, basis, form, p) for form in forms)


def analyze() -> dict:
    sols = mod13_solutions(P)
    base_by_signature = Counter(boundary_signature(z) for z in sols)
    cover_flags = Counter((boundary_signature(z), cover_boundary_flags(z)) for z in sols)

    summary = Counter()
    by_signature: dict[tuple[str, ...], Counter] = defaultdict(Counter)
    dimensions = Counter()
    examples: dict[tuple[str, ...], tuple[int, int, int, int, int, int, int]] = {}

    for z in sols:
        sig = boundary_signature(z)
        lift = affine_lift_space(z, P)
        if lift is None:
            summary["no_p2_lift"] += 1
            by_signature[sig]["no_p2_lift"] += 1
            continue

        y0, basis = lift
        summary["has_p2_lift"] += 1
        by_signature[sig]["has_p2_lift"] += 1
        dimensions[len(basis)] += 1

        branch_forms = [coord_form(i) for i in range(4) if z[i] == 0]
        if branch_forms and forms_can_be_simultaneously_nonzero(y0, basis, branch_forms, P):
            summary["can_make_zero_coords_exact1"] += 1
            by_signature[sig]["can_make_zero_coords_exact1"] += 1

        if z[6] == 0 and forms_can_be_simultaneously_nonzero(y0, basis, [coord_form(6)], P):
            summary["can_make_v_exact1"] += 1
            by_signature[sig]["can_make_v_exact1"] += 1

        disc_zero = (z[5] * z[5] - 4 * z[6] * z[6]) % P == 0
        if disc_zero and forms_can_be_simultaneously_nonzero(y0, basis, [disc_form(z[5], z[6])], P):
            summary["can_make_disc_exact1"] += 1
            by_signature[sig]["can_make_disc_exact1"] += 1

        all_forms = list(branch_forms)
        if z[6] == 0:
            all_forms.append(coord_form(6))
        if disc_zero:
            all_forms.append(disc_form(z[5], z[6]))
        if all_forms and forms_can_be_simultaneously_nonzero(y0, basis, all_forms, P):
            summary["can_satisfy_all_first_order"] += 1
            by_signature[sig]["can_satisfy_all_first_order"] += 1
            examples.setdefault(sig, z)

    return {
        "solutions": sols,
        "base_by_signature": base_by_signature,
        "cover_flags": cover_flags,
        "summary": summary,
        "by_signature": by_signature,
        "dimensions": dimensions,
        "examples": examples,
    }


def write_report(path: Path) -> None:
    result = analyze()
    lines: list[str] = []
    sols = result["solutions"]
    lines.append("Normalized 13-adic boundary analysis for M(2,2,2,8) plus 3-torsion")
    lines.append(f"mod13_solutions {len(sols)}")
    lines.append("")
    lines.append("Base boundary signatures among mod-13 cover solutions:")
    for sig, count in result["base_by_signature"].most_common():
        lines.append(f"  {sig}: {count}")
    lines.append("")
    lines.append("Cover boundary flags by base signature:")
    for (sig, flags), count in result["cover_flags"].most_common():
        lines.append(f"  {sig} | {flags}: {count}")
    lines.append("")
    lines.append(f"p2_lift_summary {dict(result['summary'])}")
    lines.append(f"p2_lift_dimensions {dict(result['dimensions'])}")
    lines.append("")
    lines.append("First-order behavior by base signature:")
    for sig, counts in sorted(
        result["by_signature"].items(),
        key=lambda item: (-item[1]["has_p2_lift"], item[0]),
    ):
        example = result["examples"].get(sig)
        lines.append(f"  {sig}: {dict(counts)} example={example}")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("  Pure collision signatures Eij have no mod-13 cover solutions.")
    lines.append("  The cover only appears above zero-boundary signatures.")
    lines.append("  Simple-zero charts Zi and mixed Zi+Ejk charts lift to 13^2,")
    lines.append("  but their zero coordinate is forced to remain zero modulo 13^2.")
    lines.append("  The only first-order live charts are the four triple-zero/one-unit")
    lines.append("  charts: Z1Z2Z3, Z1Z2Z4, Z1Z3Z4, Z2Z3Z4.")
    lines.append("  These charts have representatives where all zero coordinates acquire")
    lines.append("  nonzero first coefficients modulo 13^2.")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="data/m2228_three_torsion_padic13_report.txt",
        help="report output path",
    )
    args = parser.parse_args()
    out = Path(args.output)
    write_report(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
