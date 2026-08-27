#!/usr/bin/env python3
"""Bounded short-vector search in one orbit-12 CK CRT disk.

The three certified local orbit-12 seeds in
``elkies22210_source_halving_local_and_search_2026_07_11.md`` give ordinary
congruences on the complete rational CK chart

    t = (r3+r4)/(r1+r2),       m = (r1+r2+r3)/r1.

For the displayed seeds their denominators are p-adic units.  Combining the
three congruences gives t=t0 and m=m0 modulo

    Q = 11^3 * 19^2 * 23^2.

A rational parameter a/b in this disk is therefore a primitive vector in
the determinant-Q lattice a-t0*b = 0 (mod Q), with b>0; likewise for m.
This script Gauss-reduces the two lattices, enumerates all their primitive
vectors in a Euclidean ball, constructs the corresponding exact CK points,
and tests the four *exact* Stoll--Zarhin orbit-12 radicands for squares.

The search concerns one selected product of three Hensel disks.  A zero
count is not a global obstruction to the orbit-12 cover.

Run, with an explicit memory and wall-time cap, for example as

    timeout 300s prlimit --as=8589934592 -- \
      python3 code/elkies22210_orbit12_ck_crt_lattice.py --radius 1000000
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import gcd, isqrt, prod


MODULI = (11**3, 19**2, 23**2)
SEEDS = (
    (1, 242, 959, 9, 120),
    (1, 248, 98, 3, 11),
    (1, 392, 118, 8, 10),
)
Q = prod(MODULI)

# The original order (0,1,2,3,4) is a poor integral chart at the 19^2
# seed: both the polynomial coordinate R1=1+t(t+2)m and the chart factor
# 1+t*m^2 have valuation one at 19, so congruences modulo 19^2 lose one
# digit after primitive normalization.  Swapping r4 and r5 keeps the marked
# pair (r1,r2) fixed and makes both quantities units at all three seeds.
CHART_PERMUTATION = (0, 1, 2, 4, 3)

# Cheap exact-square prefilters.  The CRT primes are intentionally omitted:
# every vector in the selected disk already passes those local tests.
SQUARE_FILTER_MODULI = (64, 63, 65, 7, 13, 17, 29, 31, 37, 41, 43)
SQUARE_FILTERS = {
    q: {x * x % q for x in range(q)} for q in SQUARE_FILTER_MODULI
}


def crt(residues: tuple[int, ...]) -> int:
    """Chinese remainder of residues against the fixed coprime moduli."""
    value = 0
    for residue, q in zip(residues, MODULI):
        cofactor = Q // q
        value += residue * cofactor * pow(cofactor, -1, q)
    return value % Q


def chart_residues(seed: tuple[int, ...], q: int) -> tuple[int, int]:
    """Return the (t,m) chart coordinates of a projective CK seed mod q."""
    r1, r2, r3, r4, _ = tuple(seed[index] for index in CHART_PERMUTATION)
    t = (r3 + r4) * pow((r1 + r2) % q, -1, q) % q
    m = (r1 + r2 + r3) * pow(r1 % q, -1, q) % q
    return t, m


def centered(value: int) -> int:
    value %= Q
    return value if value <= Q // 2 else value - Q


def nearest_integer(numerator: int, denominator: int) -> int:
    """Nearest integer to numerator/denominator, with deterministic ties."""
    quotient, remainder = divmod(numerator, denominator)
    return quotient + (2 * remainder > denominator)


def gauss_reduce(residue: int) -> tuple[tuple[int, int], tuple[int, int]]:
    """Gauss-reduce [(Q,0),(residue,1)], a basis of a-residue*b=0 mod Q."""
    first = (Q, 0)
    second = (residue, 1)
    while True:
        n1 = first[0] * first[0] + first[1] * first[1]
        n2 = second[0] * second[0] + second[1] * second[1]
        if n2 < n1:
            first, second = second, first
            n1 = n2
        mu = nearest_integer(
            first[0] * second[0] + first[1] * second[1], n1
        )
        if mu == 0:
            break
        second = (second[0] - mu * first[0], second[1] - mu * first[1])
    assert first[0] * second[1] - first[1] * second[0] in (Q, -Q)
    return first, second


def short_rationals(residue: int, radius: int) -> list[tuple[int, int]]:
    """Enumerate primitive (numerator,positive denominator) vectors in a disk.

    Since radius<Q/2, a denominator determines at most one centered
    numerator.  This enumerates the disk without allocating a large lattice
    box or invoking a general short-vector package.
    """
    assert 0 < radius < Q // 2
    radius_squared = radius * radius
    result = []
    for denominator in range(1, radius + 1):
        numerator = centered(residue * denominator)
        if numerator * numerator + denominator * denominator > radius_squared:
            continue
        if gcd(abs(numerator), denominator) != 1:
            continue
        # Primitivity plus the congruence imply a p-adic unit denominator,
        # but retain the assertion because that is essential to this chart.
        assert gcd(denominator, Q) == 1
        assert (numerator - residue * denominator) % Q == 0
        result.append((numerator, denominator))
    result.sort(key=lambda pair: pair[0] * pair[0] + pair[1] * pair[1])
    return result


def ck_tuple(a: int, b: int, c: int, d: int) -> tuple[int, ...]:
    """Primitive integral CK tuple for t=a/b and m=c/d."""
    bb = b * b
    dd = d * d
    values = (
        bb * dd + a * (a + 2 * b) * c * d,
        a * c * (b * c - d * (a + 2 * b)),
        -bb * dd + bb * c * d + a * (a + b) * c * c,
        b * (a + b) * dd - bb * c * d - a * b * c * c,
        -(a + b) * (b * dd + a * c * c),
    )
    common = 0
    for value in values:
        common = gcd(common, abs(value))
    assert common
    values = tuple(value // common for value in values)
    if values[0] < 0:
        values = tuple(-value for value in values)
    assert sum(values) == 0
    assert sum(value**3 for value in values) == 0
    return values


def smooth_ck(values: tuple[int, ...]) -> bool:
    if any(value == 0 for value in values):
        return False
    squares = [value * value for value in values]
    return len(set(squares)) == 5


def radicands(values: tuple[int, ...]) -> tuple[int, ...]:
    """Exact orbit-12 radicands for the fixed marked pair {r1^2,r2^2}."""
    aa = [value * value for value in values]
    return (
        -(aa[0] - aa[2]) * (aa[0] - aa[3]) * (aa[0] - aa[4]),
        (aa[2] - aa[1]) * (aa[0] - aa[3]) * (aa[0] - aa[4]),
        (aa[3] - aa[1]) * (aa[0] - aa[2]) * (aa[0] - aa[4]),
        (aa[4] - aa[1]) * (aa[0] - aa[2]) * (aa[0] - aa[3]),
    )


def positive_square(value: int) -> bool:
    if value <= 0:
        return False
    for q, residues in SQUARE_FILTERS.items():
        if value % q not in residues:
            return False
    root = isqrt(value)
    return root * root == value


def verify_seed_chart(
    t_residues: tuple[int, ...], m_residues: tuple[int, ...]
) -> None:
    """Check that the chart residues really recover the displayed seeds."""
    for q, seed, t, m in zip(MODULI, SEEDS, t_residues, m_residues):
        values = ck_tuple(t, 1, m, 1)
        r1_inverse = pow(values[0] % q, -1, q)
        normalized = tuple(value * r1_inverse % q for value in values)
        permuted_seed = tuple(seed[index] % q for index in CHART_PERMUTATION)
        assert normalized == permuted_seed
        p = 11 if q == 11**3 else 19 if q == 19**2 else 23
        assert (1 + t * (t + 2) * m) % p
        assert (1 + t * m * m) % p
        square_residues = {root * root % q for root in range(q)}
        local_radicands = tuple(value % q for value in radicands(permuted_seed))
        assert all(value != 0 and value in square_residues for value in local_radicands)
        valuations = []
        for value in local_radicands:
            valuation = 0
            while value % p == 0:
                value //= p
                valuation += 1
            valuations.append(valuation)
        if p == 11:
            # Swapping r4,r5 permutes the documented (2,2,2,0) signature.
            assert sorted(valuations) == [0, 2, 2, 2]
        else:
            assert valuations == [0, 0, 0, 0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--radius",
        type=int,
        default=250_000,
        help="Euclidean radius for each numerator/denominator lattice",
    )
    parser.add_argument(
        "--sample-limit",
        type=int,
        default=8,
        help="maximum number of partial or exact samples to print",
    )
    args = parser.parse_args()
    if not (0 < args.radius < Q // 2):
        parser.error(f"require 1 <= radius < {Q // 2}")

    local = tuple(chart_residues(seed, q) for seed, q in zip(SEEDS, MODULI))
    t_residues = tuple(pair[0] for pair in local)
    m_residues = tuple(pair[1] for pair in local)
    t0 = crt(t_residues)
    m0 = crt(m_residues)
    verify_seed_chart(t_residues, m_residues)

    t_values = short_rationals(t0, args.radius)
    m_values = short_rationals(m0, args.radius)

    pair_count = 0
    smooth_count = 0
    sign_feasible_count = 0
    masks: Counter[int] = Counter()
    partial_samples = []
    exact_samples = []

    for a, b in t_values:
        for c, d in m_values:
            pair_count += 1
            values = ck_tuple(a, b, c, d)
            if not smooth_ck(values):
                continue
            smooth_count += 1
            gs = radicands(values)
            if all(value > 0 for value in gs):
                sign_feasible_count += 1
            mask = 0
            for index, value in enumerate(gs):
                if positive_square(value):
                    mask |= 1 << index
            masks[mask] += 1
            sample = (a, b, c, d, values)
            if mask == 15 and len(exact_samples) < args.sample_limit:
                exact_samples.append(sample)
            elif mask and len(partial_samples) < args.sample_limit:
                partial_samples.append((mask, sample))

    print("ELKIES22210_ORBIT12_CK_CRT_LATTICE")
    print("moduli", MODULI)
    print("modulus_product", Q)
    print("chart_permutation_zero_based", CHART_PERMUTATION)
    for q, seed, t, m in zip(MODULI, SEEDS, t_residues, m_residues):
        print("local_seed", q, seed, "t", t, "m", m)
    print("chart_congruence_t", t0)
    print("chart_congruence_m", m0)
    print("t_gauss_basis", gauss_reduce(t0))
    print("m_gauss_basis", gauss_reduce(m0))
    print("radius", args.radius)
    print("t_short_vectors", len(t_values))
    print("m_short_vectors", len(m_values))
    print("t_shortest", t_values[: min(5, len(t_values))])
    print("m_shortest", m_values[: min(5, len(m_values))])
    print("parameter_pairs", pair_count)
    print("smooth_ck_points", smooth_count)
    print("sign_feasible_points", sign_feasible_count)
    print(
        "exact_square_mask_counts",
        " ".join(f"{mask}:{masks[mask]}" for mask in sorted(masks)),
    )
    for mask, sample in partial_samples:
        a, b, c, d, values = sample
        print(
            "PARTIAL", "mask", mask, "t", (a, b), "m", (c, d),
            "r_chart_order", values,
        )
    for sample in exact_samples:
        a, b, c, d, values = sample
        print(
            "EXACT_COVER_POINT", "t", (a, b), "m", (c, d),
            "r_chart_order", values,
        )
    print("exact_cover_hits", masks[15])
    print("selected_disk_only", True)
    print("DONE")


if __name__ == "__main__":
    main()
