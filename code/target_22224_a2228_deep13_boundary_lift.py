#!/usr/bin/env python3
"""Deep corrected p-adic lifts of off-rectangle A2228+3 boundary charts.

This is an incidence lift, not the obsolete factor-of-two elimination.  The
eleven variables are

    a,b,c,d, r0,r1,r2,r3, L,U,v,

with four explicit square-root equations for the full A(2,2,2,8) cover and
the three corrected cleared cubic-contact equations.  We also fix the first
projective branch coordinate to 1.  Starting at the corrected off-rectangle
11- and 13-boundary charts, a beam Hensel search follows compatible sheets
through p^5.  Every retained node is checked against all eight exact integer
equations modulo the current p-power.

The output is evidence for genuine deep local branches.  It does not claim
that a finite p^5 lift alone proves an infinite Q_p branch; the rank and
strict-boundary valuation data are printed so that a later blow-up can supply
that final formal certificate.
"""

from __future__ import annotations

import argparse
import itertools
import random
from collections import Counter
from math import prod


SEEDS = (
    (18, 686, 4932, 8631),
    (833, 2529, 4496, 4913),
    (17, 4208, 4471, 5329),
    (25, 264, 936, 6864),
)


def radicands(base: tuple[int, int, int, int] | list[int]) -> tuple[int, ...]:
    a, b, c, d = base
    return (
        a * b * c * d,
        a * (a + b) * (a + c) * (a + d),
        b * (b + a) * (b + c) * (b + d),
        c * (c + a) * (c + b) * (c + d),
    )


def elementary(base: tuple[int, int, int, int] | list[int]) -> tuple[int, ...]:
    a, b, c, d = (z * z for z in base)
    return (
        a + b + c + d,
        a * b + a * c + a * d + b * c + b * d + c * d,
        a * b * c + a * b * d + a * c * d + b * c * d,
        a * b * c * d,
    )


def equations(z: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    base = z[:4]
    roots = z[4:8]
    L, U, v = z[8:11]
    e1, e2, e3, e4 = elementary(base)
    M = L * L
    PP = 4 * M * e1 + 12 * (U * U + v * v) - (M + 3 * U) ** 2
    contacts = (
        (M + 3 * U) * PP + 16 * v**3 - 8 * U**3 - 48 * U * v * v - 8 * M * e2,
        PP * PP + 64 * (M + 3 * U) * v**3
        - 192 * (U * U * v * v + v**4) - 64 * M * e3,
        PP * v**3 - 12 * U * v**4 - 4 * M * e4,
    )
    return tuple(roots[i] * roots[i] - radicands(base)[i] for i in range(4)) + contacts + (base[0] - 1,)


def roots_mod(x: int, p: int) -> tuple[int, ...]:
    return tuple(r for r in range(p) if (r * r - x) % p == 0)


def normalized_seed(seed: tuple[int, ...], p: int) -> tuple[int, ...]:
    scale = pow(seed[0] % p, -1, p)
    return tuple(z * scale % p for z in seed)


def contact_witnesses(base: tuple[int, ...], p: int):
    # L=0 is essential at p=13: a rational p-adic contact can have
    # v_p(L)>0 even though L itself is nonzero over Q_p.  The depth lift
    # below tests whether such a special-fibre point leaves the cleared
    # boundary or is merely an extraneous L=0 component.
    for L in range(p):
        for U in range(p):
            for v in range(p):
                z = base + (0, 0, 0, 0, L, U, v)
                if all(e % p == 0 for e in equations(z)[4:7]):
                    yield (L, U, v)


def incidence_points(seed: tuple[int, ...], p: int):
    base = normalized_seed(seed, p)
    rr = radicands(base)
    sheets = [roots_mod(x, p) for x in rr]
    if any(not s for s in sheets):
        return base, []
    contacts = list(contact_witnesses(base, p))
    points = [base + roots + contact for roots in itertools.product(*sheets) for contact in contacts]
    return base, points


def jacobian_mod(z: tuple[int, ...], p: int) -> list[list[int]]:
    """Return the exact polynomial Jacobian modulo p.

    A unit finite difference is not a derivative for these high-degree
    equations.  The p-step quotient is congruent to the derivative modulo p:

        (F(z + p e_j) - F(z))/p = dF/dz_j (z)  (mod p).
    """
    f0 = equations(z)
    out = [[0] * len(z) for _ in f0]
    for j in range(len(z)):
        zz = list(z)
        zz[j] += p
        fj = equations(zz)
        for i in range(len(f0)):
            diff = fj[i] - f0[i]
            assert diff % p == 0
            out[i][j] = (diff // p) % p
    return out


def affine_solve(mat: list[list[int]], rhs: list[int], p: int):
    nrows = len(mat)
    nvars = len(mat[0])
    aug = [[x % p for x in mat[i]] + [rhs[i] % p] for i in range(nrows)]
    pivots: list[int] = []
    row = 0
    for col in range(nvars):
        pivot = next((r for r in range(row, nrows) if aug[r][col]), None)
        if pivot is None:
            continue
        aug[row], aug[pivot] = aug[pivot], aug[row]
        q = pow(aug[row][col], -1, p)
        aug[row] = [q * x % p for x in aug[row]]
        for r in range(nrows):
            if r != row and aug[r][col]:
                q = aug[r][col]
                aug[r] = [(aug[r][c] - q * aug[row][c]) % p for c in range(nvars + 1)]
        pivots.append(col)
        row += 1
    for r in range(row, nrows):
        if all(aug[r][c] == 0 for c in range(nvars)) and aug[r][-1]:
            return row, None
    free = [c for c in range(nvars) if c not in pivots]
    particular = [0] * nvars
    for r, c in enumerate(pivots):
        particular[c] = aug[r][-1]
    basis = []
    for fc in free:
        v = [0] * nvars
        v[fc] = 1
        for r, c in enumerate(pivots):
            v[c] = -aug[r][fc] % p
        basis.append(tuple(v))
    return row, (tuple(particular), tuple(basis))


def coefficient_vectors(p: int, dim: int, cap: int, rng: random.Random):
    total = p**dim
    if total <= cap:
        yield from itertools.product(range(p), repeat=dim)
        return
    yielded = set()
    deterministic = [tuple([0] * dim)]
    for i in range(dim):
        for q in (1, p - 1):
            z = [0] * dim
            z[i] = q
            deterministic.append(tuple(z))
    for z in deterministic:
        if z not in yielded:
            yielded.add(z)
            yield z
    while len(yielded) < cap:
        z = tuple(rng.randrange(p) for _ in range(dim))
        if z not in yielded:
            yielded.add(z)
            yield z


def valuation(n: int, p: int, cap: int) -> int:
    if n == 0:
        return cap
    n = abs(n)
    v = 0
    while v < cap and n % p == 0:
        n //= p
        v += 1
    return v


def boundary_profile(z: tuple[int, ...], p: int, depth: int) -> tuple[int, ...]:
    base = z[:4]
    factors = list(base)
    for i in range(4):
        for j in range(i + 1, 4):
            factors.extend((base[i] - base[j], base[i] + base[j]))
    return tuple(sorted(valuation(x, p, depth) for x in factors if valuation(x, p, depth)))


def resolution_profile(z: tuple[int, ...], p: int, depth: int) -> tuple[int, ...]:
    """Boundary valuations, including L when the contact starts at L=0."""
    return tuple(sorted(boundary_profile(z, p, depth) + ((valuation(z[8], p, depth),) if z[8] % p == 0 else ())))


def lift_one_level(z: tuple[int, ...], p: int, level: int, cap: int, rng: random.Random):
    mod = p**level
    values = equations(z)
    if any(x % mod for x in values):
        raise AssertionError((p, level, z, values))
    rhs = [-(x // mod) % p for x in values]
    rank, space = affine_solve(jacobian_mod(tuple(x % p for x in z), p), rhs, p)
    if space is None:
        return rank, (), 0
    particular, basis = space
    out = []
    next_mod = mod * p
    for coeffs in coefficient_vectors(p, len(basis), cap, rng):
        delta = [
            (particular[i] + sum(c * b[i] for c, b in zip(coeffs, basis))) % p
            for i in range(len(z))
        ]
        zz = tuple((z[i] + mod * delta[i]) % next_mod for i in range(len(z)))
        if all(x % next_mod == 0 for x in equations(zz)):
            out.append(zz)
    return rank, tuple(out), len(basis)


def run_chart(
    seed: tuple[int, ...],
    p: int,
    depth: int,
    beam: int,
    children_by_level: tuple[int, ...],
    rng: random.Random,
    incidence_offset: int,
    incidence_limit: int,
):
    base, points = incidence_points(seed, p)
    total_points = len(points)
    if incidence_limit:
        points = points[incidence_offset : incidence_offset + incidence_limit]
    print(
        "CHART_START",
        "seed",
        seed,
        "p",
        p,
        "base",
        base,
        "incidence_modp_total",
        total_points,
        "incidence_modp_used",
        len(points),
        "incidence_offset",
        incidence_offset,
    )
    if not points:
        return {"modp": 0, "final": 0}
    rank_hist = Counter()
    for z in points:
        rank, _ = affine_solve(jacobian_mod(z, p), [0] * 8, p)
        rank_hist[rank] += 1
    print("MODP_RANK_HIST", sorted(rank_hist.items()))
    nodes = tuple(points)
    level_counts = [len(nodes)]
    dims = Counter()
    for level in range(1, depth):
        child_cap = children_by_level[min(level - 1, len(children_by_level) - 1)]
        candidates = {}
        inconsistent = 0
        for z in nodes:
            rank, lifts, dim = lift_one_level(z, p, level, child_cap, rng)
            dims[(level, rank, dim)] += 1
            if not lifts:
                inconsistent += 1
            for zz in lifts:
                candidates[zz] = resolution_profile(zz, p, level + 1)
        # Prefer branches which resolve as many boundary factors as possible
        # and at the shallowest possible valuation.
        # Randomize ties.  Lexicographic truncation badly biases singular
        # charts: second-order consistency can depend on high tangent digits
        # even when every candidate has the same boundary profile.
        tie_break = {z: rng.random() for z in candidates}
        ordered = sorted(
            candidates,
            key=lambda z: (len(candidates[z]), sum(candidates[z]), tie_break[z]),
        )
        nodes = tuple(ordered[:beam])
        level_counts.append(len(nodes))
        print(
            "LIFT_LEVEL",
            level + 1,
            "input",
            level_counts[-2],
            "candidates",
            len(candidates),
            "kept",
            len(nodes),
            "inconsistent",
            inconsistent,
            "child_cap",
            child_cap,
            "best_boundary",
            candidates[nodes[0]] if nodes else None,
        )
        if not nodes:
            break
    sample = nodes[0] if nodes else None
    if sample:
        modulus = p**depth
        assert all(x % modulus == 0 for x in equations(sample))
        print(
            "FINAL_SAMPLE",
            "modulus",
            modulus,
            "z",
            sample,
            "boundary",
            boundary_profile(sample, p, depth),
            "L_valuation",
            valuation(sample[8], p, depth),
        )
    print("LIFT_DIM_HIST", sorted(dims.items()))
    return {"modp": len(points), "levels": level_counts, "final": len(nodes), "sample": sample}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="11,13")
    parser.add_argument("--depth", type=int, default=5)
    parser.add_argument("--beam", type=int, default=3000)
    parser.add_argument("--children", type=int, default=800)
    parser.add_argument(
        "--children-by-level",
        default="",
        help="Comma-separated caps for lifts p^n -> p^(n+1); overrides --children",
    )
    parser.add_argument("--incidence-offset", type=int, default=0)
    parser.add_argument(
        "--incidence-limit",
        type=int,
        default=0,
        help="Use this many mod-p incidence sheets (0 means all)",
    )
    parser.add_argument("--seed-indices", default="1,2,3,4")
    parser.add_argument("--random-seed", type=int, default=22224)
    args = parser.parse_args()
    rng = random.Random(args.random_seed)
    children_by_level = (
        tuple(int(x) for x in args.children_by_level.split(","))
        if args.children_by_level
        else (args.children,)
    )
    primes = [int(x) for x in args.primes.split(",")]
    indices = [int(x) - 1 for x in args.seed_indices.split(",")]
    summary = {}
    for i in indices:
        for p in primes:
            summary[(i + 1, p)] = run_chart(
                SEEDS[i],
                p,
                args.depth,
                args.beam,
                children_by_level,
                rng,
                args.incidence_offset,
                args.incidence_limit,
            )
    print("DEEP_BOUNDARY_DONE", summary)


if __name__ == "__main__":
    main()
