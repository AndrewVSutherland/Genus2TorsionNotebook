#!/usr/bin/env python3
"""First-order 13-adic lifts of the full A(2,2,2,8)+3 boundary.

This imports the finite raw-boundary calculation and linearizes both:

* the three eliminated cubic-contact equations; and
* every full-cover radicand that vanishes modulo 13.

For a vanishing radicand, being a 13-adic square forces its value to vanish
modulo 13^2.  Nonzero square radicands lift automatically by Hensel.  The
output distinguishes mere affine lifts from lifts that make every zero base
coordinate have valuation exactly one and resolve every collision between
nonzero branch parameters at the first order.
"""

from __future__ import annotations

from collections import Counter
from itertools import product

import target_22224_full_boundary13 as base


P = base.P


def elementary(vals: list[int]) -> tuple[int, int, int, int]:
    a, b, c, d = vals
    A, B, C, D = a * a, b * b, c * c, d * d
    return (
        A + B + C + D,
        A * B + A * C + A * D + B * C + B * D + C * D,
        A * B * C + A * B * D + A * C * D + B * C * D,
        A * B * C * D,
    )


def contact_equations(state: list[int]) -> list[int]:
    vals = state[:4]
    L, U, v = state[4:]
    e1, e2, e3, e4 = elementary(vals)
    M = L * L
    A = 2 * M * e1 + 6 * (U * U + v * v) - (M + 3 * U) ** 2
    return [
        (M + 3 * U) * A + 8 * v**3 - 4 * e2 * M - 4 * U**3 - 24 * U * v * v,
        A * A + 16 * (M + 3 * U) * v**3 - 16 * e3 * M
        - 48 * (U * U * v * v + v**4),
        A * v**3 - 2 * e4 * M - 6 * U * v**4,
    ]


def integer_radicands(state: list[int]) -> list[int]:
    a, b, c, d = state[:4]
    return [
        a * b * c * d,
        a * (a + b) * (a + c) * (a + d),
        b * (b + a) * (b + c) * (b + d),
        c * (c + a) * (c + b) * (c + d),
    ]


def linearized_row(func, state: list[int]) -> tuple[list[int], int]:
    value = func(state)
    assert value % P == 0
    row = []
    for i in range(len(state)):
        lifted = state[:]
        lifted[i] += P
        row.append(((func(lifted) - value) // P) % P)
    return row, (-(value // P)) % P


def rref_solve(rows: list[list[int]], rhs: list[int], nvars: int):
    mat = [[x % P for x in row] + [b % P] for row, b in zip(rows, rhs)]
    pivot_cols: list[int] = []
    r = 0
    for c in range(nvars):
        pivot = next((i for i in range(r, len(mat)) if mat[i][c]), None)
        if pivot is None:
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        z = pow(mat[r][c], P - 2, P)
        mat[r] = [(z * x) % P for x in mat[r]]
        for i in range(len(mat)):
            if i != r and mat[i][c]:
                z = mat[i][c]
                mat[i] = [(x - z * y) % P for x, y in zip(mat[i], mat[r])]
        pivot_cols.append(c)
        r += 1
        if r == len(mat):
            break
    for row in mat:
        if all(x == 0 for x in row[:nvars]) and row[nvars] != 0:
            return None
    free = [c for c in range(nvars) if c not in pivot_cols]
    particular = [0] * nvars
    for i, c in enumerate(pivot_cols):
        particular[c] = mat[i][nvars]
    basis = []
    for f in free:
        vec = [0] * nvars
        vec[f] = 1
        for i, c in enumerate(pivot_cols):
            vec[c] = (-mat[i][f]) % P
        basis.append(vec)
    return particular, basis


def affine_vectors(particular: list[int], basis: list[list[int]]):
    for coeffs in product(range(P), repeat=len(basis)):
        yield [
            (particular[i] + sum(c * vec[i] for c, vec in zip(coeffs, basis))) % P
            for i in range(len(particular))
        ]


def first_nonzero(vals: tuple[int, int, int, int]) -> int:
    return next(i for i, x in enumerate(vals) if x)


def resolution_flags(vals: tuple[int, int, int, int], delta: list[int]) -> tuple[bool, bool, bool]:
    # Every zero parameter should have valuation exactly one.
    zero_exact = not any(vals[i] == 0 and delta[i] == 0 for i in range(4))
    # Every collision between two nonzero square branches should resolve at p^2.
    nonzero_collision_resolved = True
    for i in range(4):
        for j in range(i + 1, 4):
            if vals[i] and vals[j] and (vals[i] * vals[i] - vals[j] * vals[j]) % P == 0:
                # Representatives lie in [0,p-1], so a relation such as
                # a=-b mod p can already contribute one p to a^2-b^2.
                deriv = (
                    (vals[i] * vals[i] - vals[j] * vals[j]) // P
                    + 2 * (vals[i] * delta[i] - vals[j] * delta[j])
                )
                if deriv % P == 0:
                    nonzero_collision_resolved = False
    # If two coordinates vanish, divide them by 13 and ask that their
    # leading square branches already separate on the exceptional divisor.
    scaled_zero_separated = True
    for i in range(4):
        for j in range(i + 1, 4):
            if vals[i] == vals[j] == 0 and (delta[i] * delta[i] - delta[j] * delta[j]) % P == 0:
                scaled_zero_separated = False
    return zero_exact, nonzero_collision_resolved, scaled_zero_separated


def main() -> None:
    by_key = base.contact_witnesses()
    totals = Counter()
    signature = Counter()
    dimensions = Counter()
    examples = []

    for vals in base.projective_representatives():
        rads = base.radicands(vals)
        if any(r not in base.SQUARES for r in rads):
            continue
        key = base.curve_key(vals)
        for witness in by_key.get(key, ()):
            if base.cover_flags(key, witness) != ("contact_open",):
                continue
            totals["presentations"] += 1
            state = list(vals) + list(witness)
            rows: list[list[int]] = []
            rhs: list[int] = []
            for eq_index in range(3):
                row, b = linearized_row(lambda s, i=eq_index: contact_equations(s)[i], state)
                rows.append(row)
                rhs.append(b)
            for rad_index, rad in enumerate(rads):
                if rad == 0:
                    row, b = linearized_row(lambda s, i=rad_index: integer_radicands(s)[i], state)
                    rows.append(row)
                    rhs.append(b)
            # Fix the projective scale by keeping the first nonzero base entry equal to 1.
            scale_row = [0] * 7
            scale_row[first_nonzero(vals)] = 1
            rows.append(scale_row)
            rhs.append(0)
            solution = rref_solve(rows, rhs, 7)
            sig = "+".join(base.labels(vals))
            if solution is None:
                totals["inconsistent"] += 1
                signature[(sig, "inconsistent")] += 1
                continue
            particular, basis = solution
            dim = len(basis)
            totals["liftable"] += 1
            dimensions[dim] += 1
            signature[(sig, "liftable")] += 1
            can_zero_exact = False
            can_collision_resolve = False
            can_scaled_separate = False
            can_progress = False
            witness_delta = None
            # Dimensions here are small; exhaustive affine-space testing also
            # catches unions of forbidden hyperplanes without heuristics.
            for delta in affine_vectors(particular, basis):
                zero_exact, collision_resolved, scaled_separated = resolution_flags(vals, delta)
                can_zero_exact |= zero_exact
                can_collision_resolve |= collision_resolved
                can_scaled_separate |= scaled_separated
                if zero_exact and collision_resolved and scaled_separated:
                    can_progress = True
                    witness_delta = delta
                    break
            if can_zero_exact:
                totals["zero_exact1"] += 1
                signature[(sig, "zero_exact1")] += 1
            if can_collision_resolve:
                totals["nonzero_collision_resolved"] += 1
                signature[(sig, "collision_resolved")] += 1
            if can_scaled_separate:
                totals["scaled_zero_separated"] += 1
                signature[(sig, "scaled_separated")] += 1
            if can_progress:
                totals["first_order_progressing"] += 1
                signature[(sig, "progressing")] += 1
                if len(examples) < 30:
                    examples.append((vals, witness, dim, witness_delta, rads))

    print("TARGET_22224_FULL_BOUNDARY13_LIFT")
    for k in (
        "presentations", "inconsistent", "liftable", "zero_exact1",
        "nonzero_collision_resolved", "scaled_zero_separated", "first_order_progressing",
    ):
        print(k, totals[k])
    print("solution_dimensions", sorted(dimensions.items()))
    print("signature_summary")
    sigs = sorted({s for s, _ in signature})
    for sig in sigs:
        print(
            " ", sig,
            "liftable", signature[(sig, "liftable")],
            "zero_exact1", signature[(sig, "zero_exact1")],
            "collision_resolved", signature[(sig, "collision_resolved")],
            "scaled_separated", signature[(sig, "scaled_separated")],
            "progressing", signature[(sig, "progressing")],
            "inconsistent", signature[(sig, "inconsistent")],
        )
    print("examples")
    for row in examples:
        print(" ", row)


if __name__ == "__main__":
    main()
