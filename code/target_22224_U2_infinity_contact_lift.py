#!/usr/bin/env python3
"""U2 p=13 infinity-chart contact lift and full-cover blow-up test.

For U2=[A,B,C,D*T^2], put z=1/T and use projective invariance to
write the infinity chart as

    [A*z^2, B*z^2, C*z^2, D].

The corrected cleared cubic-contact equations have two smooth unit-L
incidences at z=0.  We exhaustively lift those contact-only branches through
13^5.  Separately, and decisively, we divide the four full-cover radicands by
z^6.  Two resulting units are nonsquares modulo 13, so no nonzero z in
13 Z_13 can lie on the full cover.  Contact-only lift rows must therefore not
be promoted to actual [2,2,2,24] local branches.
"""

from __future__ import annotations

import argparse
import itertools
from collections import Counter
from pathlib import Path

P = 13
A, B, C, D = -1071, -1054, 1116, 1134


def elementary(base):
    a, b, c, d = (x * x for x in base)
    return (
        a + b + c + d,
        a * b + a * c + a * d + b * c + b * d + c * d,
        a * b * c + a * b * d + a * c * d + b * c * d,
        a * b * c * d,
    )


def contact_equations(x):
    z, L, U, v = x
    base = (A * z * z, B * z * z, C * z * z, D)
    e1, e2, e3, e4 = elementary(base)
    M = L * L
    PP = 4 * M * e1 + 12 * (U * U + v * v) - (M + 3 * U) ** 2
    return (
        (M + 3 * U) * PP + 16 * v**3 - 8 * U**3
        - 48 * U * v * v - 8 * M * e2,
        PP * PP + 64 * (M + 3 * U) * v**3
        - 192 * (U * U * v * v + v**4) - 64 * M * e3,
        PP * v**3 - 12 * U * v**4 - 4 * M * e4,
    )


def radicand_units(z):
    """R_i/z^6 on the U2 infinity chart."""
    return (
        A * B * C * D,
        A * (A + B) * (A + C) * (D + A * z * z),
        B * (B + A) * (B + C) * (D + B * z * z),
        C * (C + A) * (C + B) * (D + C * z * z),
    )


def rank_and_space(mat, rhs, p=P):
    rows = [[q % p for q in row] + [rhs[i] % p]
            for i, row in enumerate(mat)]
    nr, nv = len(rows), len(mat[0])
    pivots = []
    r = 0
    for c in range(nv):
        pivot = next((i for i in range(r, nr) if rows[i][c]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inv = pow(rows[r][c], -1, p)
        rows[r] = [inv * q % p for q in rows[r]]
        for i in range(nr):
            if i != r and rows[i][c]:
                q = rows[i][c]
                rows[i] = [(rows[i][j] - q * rows[r][j]) % p
                           for j in range(nv + 1)]
        pivots.append(c)
        r += 1
    if any(all(rows[i][j] == 0 for j in range(nv)) and rows[i][-1]
           for i in range(r, nr)):
        return r, None
    free = [j for j in range(nv) if j not in pivots]
    particular = [0] * nv
    for i, c in enumerate(pivots):
        particular[c] = rows[i][-1]
    basis = []
    for fc in free:
        q = [0] * nv
        q[fc] = 1
        for i, c in enumerate(pivots):
            q[c] = -rows[i][fc] % p
        basis.append(tuple(q))
    return r, (tuple(particular), tuple(basis))


def jacobian_mod(x, p=P):
    x = tuple(q % p for q in x)
    f0 = contact_equations(x)
    out = [[0] * 4 for _ in range(3)]
    # The p-step quotient is congruent to the exact polynomial derivative.
    for j in range(4):
        y = list(x)
        y[j] += p
        fy = contact_equations(y)
        for i in range(3):
            diff = fy[i] - f0[i]
            assert diff % p == 0
            out[i][j] = (diff // p) % p
    return out


def witnesses_mod_p():
    ans = []
    for L, U, v in itertools.product(range(P), repeat=3):
        x = (0, L, U, v)
        if all(q % P == 0 for q in contact_equations(x)):
            rank, _ = rank_and_space(jacobian_mod(x), [0, 0, 0])
            ans.append((x, rank))
    return ans


def lift(nodes, level):
    """Exhaustively lift smooth rank-three contact nodes one digit."""
    mod = P**level
    nxt = set()
    ranks = Counter()
    for x in nodes:
        vals = contact_equations(x)
        assert all(q % mod == 0 for q in vals)
        rhs = [-(q // mod) % P for q in vals]
        rank, space = rank_and_space(jacobian_mod(x), rhs)
        ranks[rank] += 1
        if rank != 3 or space is None:
            continue
        particular, basis = space
        for coeff in itertools.product(range(P), repeat=len(basis)):
            delta = [
                (particular[j] + sum(c * b[j]
                                     for c, b in zip(coeff, basis))) % P
                for j in range(4)
            ]
            y = tuple((x[j] + mod * delta[j]) % (mod * P)
                      for j in range(4))
            if all(q % (mod * P) == 0 for q in contact_equations(y)):
                nxt.add(y)
    return nxt, ranks


def valuation(n, p=P, cap=99):
    if n == 0:
        return cap
    v = 0
    while v < cap and n % p == 0:
        n //= p
        v += 1
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", type=int, default=5)
    ap.add_argument(
        "--output",
        default="results/target_22224_U2_infinity_contact_lift.tsv",
    )
    ap.add_argument(
        "--log",
        default="results/target_22224_U2_infinity_contact_lift.log",
    )
    args = ap.parse_args()

    unit_residues = tuple(q % P for q in radicand_units(0))
    square_residues = {x * x % P for x in range(P)}
    unit_square = tuple(q in square_residues for q in unit_residues)
    assert unit_residues == (9, 7, 11, 4)
    assert unit_square == (True, False, False, True)

    witnesses = witnesses_mod_p()
    smooth = [x for x, rank in witnesses if x[1] and rank == 3]
    singular_unit = [(x, rank) for x, rank in witnesses
                     if x[1] and rank != 3]
    lzero = [(x, rank) for x, rank in witnesses if not x[1]]
    assert smooth == [(0, 4, 1, 3), (0, 9, 1, 3)]

    lines = [
        f"U2_INFINITY_CONTACT_START depth={args.depth}",
        f"chart=[{A}*z^2,{B}*z^2,{C}*z^2,{D}] z_in_13Z13",
        f"radicand_units_mod13={unit_residues} square={unit_square}",
        "FULL_COVER_BLOWUP_OBSTRUCTION="
        "R1/z^6=7 and R2/z^6=11 are nonsquares mod13",
        f"contact_mod13_witnesses={len(witnesses)} smooth_unit={smooth}",
        f"contact_singular_unit={singular_unit}",
        f"contact_Lzero_count={len(lzero)} contact_Lzero={lzero}",
    ]

    nodes = set(smooth)
    rows = []
    for level in range(1, args.depth + 1):
        mod = P**level
        zvals = Counter(x[0] % mod for x in nodes)
        hist = Counter(min(level, valuation(z, cap=level)) for z in zvals)
        lines.append(
            f"contact_level={level} modulus={mod} nodes={len(nodes)} "
            f"zclasses={len(zvals)} zvaluation_hist={dict(hist)}"
        )
        for z, count in sorted(zvals.items()):
            rows.append((level, mod, z, count, 0 if z else -1))
        if level < args.depth:
            nodes, ranks = lift(nodes, level)
            lines.append(
                f"contact_lift_to={level+1} source_rank_hist={dict(ranks)}"
            )

    out = Path(args.output)
    with out.open("w") as f:
        f.write("level\tmodulus\tz\tcontact_incidences\tfull_cover_possible\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")
    lines.extend([
        "SCOPE contact branches are exhaustive only for the two smooth unit-L "
        "components; the full-cover nonsquare obstruction is complete for every "
        "nonzero z in 13 Z_13 and does not depend on L,U,v",
        "ACTUAL_U2_INFINITY_FULL_COVER_BRANCHES=0",
        f"output={out}",
    ])
    text = "\n".join(lines) + "\n"
    Path(args.log).write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
