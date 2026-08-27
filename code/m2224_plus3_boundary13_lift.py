#!/usr/bin/env python3
"""First-order 13-adic lift check for M(2,2,2,4)+3 boundary contacts.

This complements m2224_plus3_boundary13_analysis.m.  It enumerates the
nondegenerate mod-13 cubic-contact solutions on the M(2,2,2,4) boundary and
linearizes the three eliminated contact equations modulo 13^2 in variables

    (a,b,c,d,L,U,v).

For each solution it records whether a p^2 lift exists and whether the boundary
labels containing the base point can be made exact to first order, e.g. a zero
coordinate can become 13 times a unit or a collision a_i = +/- a_j can split at
order 13.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product

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


def trim(poly: list[int]) -> list[int]:
    while poly and poly[-1] % P == 0:
        poly.pop()
    return [c % P for c in poly]


def poly_mod(a: list[int], b: list[int], p: int = P) -> list[int]:
    a = trim(a[:])
    b = trim(b[:])
    if not b:
        raise ZeroDivisionError
    db = len(b) - 1
    ib = inv(b[-1], p)
    while len(a) >= len(b) and a:
        coeff = a[-1] * ib % p
        shift = len(a) - len(b)
        for i in range(db + 1):
            a[shift + i] = (a[shift + i] - coeff * b[i]) % p
        trim(a)
    return a


def poly_gcd(a: list[int], b: list[int], p: int = P) -> list[int]:
    a = trim(a[:])
    b = trim(b[:])
    while b:
        a, b = b, poly_mod(a, b, p)
    if not a:
        return []
    ia = inv(a[-1], p)
    return [(c * ia) % p for c in a]


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


def contact_eqs_int(z: tuple[int, int, int, int, int, int, int]) -> list[int]:
    _a, _b, _c, _d, L, U, v = z
    e1, e2, e3, e4 = coeffs_int(z)
    M = L * L
    A = 2 * M * e1 + 6 * (U * U + v * v) - (M + 3 * U) ** 2
    return [
        (M + 3 * U) * A + 8 * v**3 - 4 * e2 * M - 4 * U**3 - 24 * U * v * v,
        A * A + 16 * (M + 3 * U) * v**3 - 16 * e3 * M - 48 * (U * U * v * v + v**4),
        A * v**3 - 2 * e4 * M - 6 * U * v**4,
    ]


def raw_contact_witnesses_by_key(p: int = P) -> dict[tuple[int, int, int, int], list[tuple[int, int, int]]]:
    out: dict[tuple[int, int, int, int], list[tuple[int, int, int]]] = defaultdict(list)
    for L, U, v, e1 in product(range(1, p), range(p), range(p), range(p)):
        M = L * L % p
        A = (2 * M * e1 + 6 * (U * U + v * v) - (M + 3 * U) ** 2) % p
        e2 = ((M + 3 * U) * A + 8 * v**3 - 4 * U**3 - 24 * U * v * v) * inv(4 * M, p)
        e3 = (A * A + 16 * (M + 3 * U) * v**3 - 48 * (U * U * v * v + v**4)) * inv(16 * M, p)
        e4 = (A * v**3 - 6 * U * v**4) * inv(2 * M, p)
        key = (e1 % p, e2 % p, e3 % p, e4 % p)
        out[key].append((L, U, v))
    return out


def boundary_labels(vals: tuple[int, int, int, int], p: int = P) -> tuple[str, ...]:
    labels: list[str] = []
    for i, x in enumerate(vals):
        if x % p == 0:
            labels.append(f"Z{i + 1}")
    for i in range(4):
        for j in range(i + 1, 4):
            tag = PAIR_LABELS[(i, j)]
            if (vals[i] - vals[j]) % p == 0:
                labels.append(tag + "+")
            if (vals[i] + vals[j]) % p == 0:
                labels.append(tag + "-")
    return tuple(sorted(labels or ["open"]))


def label_key(labels: tuple[str, ...]) -> str:
    return "+".join(labels)


def cover_open(key: tuple[int, int, int, int], wit: tuple[int, int, int], p: int = P) -> bool:
    _L, U, v = wit
    if v % p == 0:
        return False
    if (U * U - 4 * v * v) % p == 0:
        return False
    e1, e2, e3, e4 = key
    f = [0, e4, e3, e2, e1, 1]
    q = [v * v % p, U % p, 1]
    return len(poly_gcd(f, q, p)) <= 1


def linear_lift_space(z: tuple[int, int, int, int, int, int, int], p: int = P):
    f0 = contact_eqs_int(z)
    rhs = [-(f // p) % p for f in f0]
    matrix = [[0] * 7 for _ in range(3)]
    for j in range(7):
        zp = list(z)
        zp[j] += p
        fj = contact_eqs_int(tuple(zp))
        for i in range(3):
            matrix[i][j] = ((fj[i] - f0[i]) // p) % p

    aug = [matrix[i][:] + [rhs[i]] for i in range(3)]
    row = 0
    pivots: list[int] = []
    for col in range(7):
        pivot = next((r for r in range(row, 3) if aug[r][col] % p), None)
        if pivot is None:
            continue
        aug[row], aug[pivot] = aug[pivot], aug[row]
        scale = inv(aug[row][col], p)
        aug[row] = [(x * scale) % p for x in aug[row]]
        for r in range(3):
            if r != row and aug[r][col] % p:
                factor = aug[r][col] % p
                aug[r] = [(aug[r][c] - factor * aug[row][c]) % p for c in range(8)]
        pivots.append(col)
        row += 1
        if row == 3:
            break
    for r in range(row, 3):
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


def label_form(label: str) -> tuple[int, ...]:
    coeffs = [0] * 7
    if label.startswith("Z"):
        coeffs[int(label[1]) - 1] = 1
    else:
        pair = label[1:3]
        i = int(pair[0]) - 1
        j = int(pair[1]) - 1
        coeffs[i] = 1
        coeffs[j] = -1 if label.endswith("+") else 1
    return tuple(c % P for c in coeffs)


def form_identically_zero(y0: tuple[int, ...], basis: list[tuple[int, ...]], form: tuple[int, ...], p: int = P) -> bool:
    def dot(y: tuple[int, ...]) -> int:
        return sum(a * b for a, b in zip(form, y)) % p
    return dot(y0) == 0 and all(dot(v) == 0 for v in basis)


def can_make_all_labels_first_order(labels: tuple[str, ...], y0: tuple[int, ...], basis: list[tuple[int, ...]]) -> bool:
    if labels == ("open",):
        return True
    forms = [label_form(label) for label in labels]
    # F_13 is not covered by fewer than 13 proper affine hyperplanes.
    return all(not form_identically_zero(y0, basis, form) for form in forms)


def main() -> None:
    by_key = raw_contact_witnesses_by_key()
    total = 0
    liftable = 0
    all_first = 0
    by_sig = Counter()
    by_sig_lift = Counter()
    by_sig_first = Counter()
    by_component = Counter()
    by_component_lift = Counter()
    by_component_resolve = Counter()
    examples: list[tuple] = []

    for vals in product(range(P), repeat=4):
        labels = boundary_labels(vals)
        if labels == ("open",):
            continue
        key = coeffs_mod(*vals)
        for wit in by_key.get(key, []):
            if not cover_open(key, wit):
                continue
            z = (*vals, *wit)
            total += 1
            sig = label_key(labels)
            by_sig[sig] += 1
            for lab in labels:
                by_component[lab] += 1
            space = linear_lift_space(z)
            if space is None:
                if len(examples) < 20:
                    examples.append((vals, sig, wit, "no_p2_lift"))
                continue
            liftable += 1
            by_sig_lift[sig] += 1
            y0, basis = space
            for lab in labels:
                if not form_identically_zero(y0, basis, label_form(lab)):
                    by_component_resolve[lab] += 1
                by_component_lift[lab] += 1
            if can_make_all_labels_first_order(labels, y0, basis):
                all_first += 1
                by_sig_first[sig] += 1
                if len(examples) < 20:
                    examples.append((vals, sig, wit, "all_first_order"))
            elif len(examples) < 20:
                examples.append((vals, sig, wit, "some_label_forced_deeper"))

    print("M(2,2,2,4)+3 p=13 cover-open first-order lift analysis")
    print("cover_open_boundary_solutions", total)
    print("p2_liftable", liftable)
    print("all_boundary_labels_first_order", all_first)
    print()
    print("component_summary")
    for lab in sorted(by_component):
        print(lab,
              "cover_open", by_component[lab],
              "p2_lift", by_component_lift[lab],
              "label_can_be_first_order", by_component_resolve[lab])
    print()
    print("signature_summary")
    for sig in sorted(by_sig):
        print(sig,
              "cover_open", by_sig[sig],
              "p2_lift", by_sig_lift[sig],
              "all_labels_first_order", by_sig_first[sig])
    print()
    print("examples")
    for ex in examples:
        print(ex)


if __name__ == "__main__":
    main()
