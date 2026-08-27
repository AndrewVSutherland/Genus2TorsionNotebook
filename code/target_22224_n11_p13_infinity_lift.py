#!/usr/bin/env python3
"""Hensel analysis of the primitive p=13 infinity chart of fiber N11.

N11 is (-960,1800,2535,-1352*T^2), with -1352=-8*13^2.  Write
T=r/s, z=s/(13*r), and divide the signed branch tuple by 13^2*r^2.
Up to a common square it becomes

    (-960*z^2, 1800*z^2, 2535*z^2, -8).

The direct cubic-contact equations are the corrected cleared equations in
(L,U,v).  This program exhaustively Hensel-lifts the *unit-L smooth* special
fiber incidences from z=0.  It records the possible z residues at every
depth.  Singular unit-L and L=0 charts are reported separately and are not
silently promoted to genuine branches.
"""

from __future__ import annotations

import argparse
import itertools
from collections import Counter, defaultdict
from pathlib import Path

P = 13
ALLOWED_Z = (0, 2, 3, 4, 6, 7, 9, 10, 11)


def elementary(base):
    A, B, C, D = (x * x for x in base)
    return (
        A + B + C + D,
        A * B + A * C + A * D + B * C + B * D + C * D,
        A * B * C + A * B * D + A * C * D + B * C * D,
        A * B * C * D,
    )


def equations(x):
    z, L, U, v = x
    base = (-960 * z * z, 1800 * z * z, 2535 * z * z, -8)
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
    f0 = equations(x)
    out = [[0] * 4 for _ in range(3)]
    # The p-step quotient is the polynomial derivative modulo p.
    for j in range(4):
        y = list(x)
        y[j] += p
        fy = equations(y)
        for i in range(3):
            assert (fy[i] - f0[i]) % p == 0
            out[i][j] = ((fy[i] - f0[i]) // p) % p
    return out


def witnesses_mod_p():
    ans = []
    for z in ALLOWED_Z:
        for L, U, v in itertools.product(range(P), repeat=3):
            x = (z, L, U, v)
            if all(q % P == 0 for q in equations(x)):
                rank, _ = rank_and_space(jacobian_mod(x), [0, 0, 0])
                ans.append((x, rank))
    return ans


def lift(nodes, level, keep_ranks=(3,)):
    """Lift solutions modulo p^level to p^(level+1)."""
    mod = P**level
    nxt = set()
    ranks = Counter()
    for x in nodes:
        f = equations(x)
        assert all(q % mod == 0 for q in f)
        rhs = [-(q // mod) % P for q in f]
        rank, space = rank_and_space(jacobian_mod(x), rhs)
        ranks[rank] += 1
        if rank not in keep_ranks or space is None:
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
            if all(q % (mod * P) == 0 for q in equations(y)):
                nxt.add(y)
    return nxt, ranks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", type=int, default=5)
    ap.add_argument("--output", default="results/target_22224_n11_p13_infinity_lift.tsv")
    ap.add_argument("--log", default="results/target_22224_n11_p13_infinity_lift.log")
    args = ap.parse_args()

    witnesses = witnesses_mod_p()
    smooth = [x for x, rank in witnesses if x[1] and rank == 3]
    singular_unit = [(x, rank) for x, rank in witnesses if x[1] and rank != 3]
    lzero = [(x, rank) for x, rank in witnesses if not x[1]]
    lines = [
        f"N11_P13_INFINITY_START depth={args.depth}",
        f"allowed_z_mod13={ALLOWED_Z}",
        f"mod13_witnesses={len(witnesses)} smooth_unit={smooth}",
        f"singular_unit={singular_unit}",
        f"Lzero_count={len(lzero)} Lzero={lzero}",
    ]
    nodes = set(smooth)
    rows = []
    for level in range(1, args.depth + 1):
        mod = P**level
        zvals = Counter(x[0] % mod for x in nodes)
        valuation = Counter()
        for z in zvals:
            q = z
            vv = 0
            while vv < level and q % P == 0:
                vv += 1
                q //= P
            valuation[vv] += 1
        lines.append(
            f"level={level} modulus={mod} nodes={len(nodes)} "
            f"zclasses={len(zvals)} zvaluation_hist={dict(valuation)}"
        )
        for z, count in sorted(zvals.items()):
            rows.append((level, mod, z, count))
        if level < args.depth:
            nodes, ranks = lift(nodes, level)
            lines.append(f"lift_to={level+1} source_rank_hist={dict(ranks)}")

    out = Path(args.output)
    with out.open("w") as f:
        f.write("level\tmodulus\tz\tincidences\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")
    lines.append(f"output={out}")
    text = "\n".join(lines) + "\n"
    Path(args.log).write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
