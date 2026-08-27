#!/usr/bin/env python3
"""Local mod-5 lift diagnostic for the contact-7 [28] -> [56] route.

The open first-halving surface over F_5 has four points.  At each of them the
order-4 class H4 is not divisible by 2 in J(F_5).  Since multiplication by 2 is
finite etale on the good-reduction abelian scheme over Z_5, no 5-adic lift in
these open neighborhoods can give a second halving.

This script still enumerates first-halving surface lifts modulo 5^k around the
four F_5 points, so the local picture is explicit.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable


BASE_POINTS = [
    # s, u, z, eps with eps represented as +/-1 over the integers.
    (2, 1, 2, 1),
    (2, 4, 2, 1),
    (3, 1, 2, -1),
    (3, 4, 2, -1),
]


def inv(a: int, modulus: int) -> int:
    return pow(a % modulus, -1, modulus)


def div_mod(a: int, b: int, modulus: int) -> int:
    return (a % modulus) * inv(b, modulus) % modulus


def coeffs_from_s_b_mod(s: int, b: int, eps: int, modulus: int) -> tuple[int, int, int, int]:
    """Return c3,c2,c1,c0 for f(X+r)/X modulo modulus.

    These formulas match code/contact7_halving_surface_enum.py and are valid
    on the open chart where S+1 is a 5-adic unit.
    """
    s %= modulus
    b %= modulus
    sgn = eps % modulus
    S = sgn * s % modulus
    den1 = (S + 1) % modulus

    c3 = (b * b - 5 * s * s - 2) % modulus

    c2_num = (
        -2 * b * b * S**4
        - 4 * b * b * S**3
        + 4 * b * b * S
        + 2 * b * b
        + 2 * b * S**5
        + 4 * b * S**4
        + 6 * b * S**3
        + 8 * b * S**2
        + 10 * b * S
        + 5 * b
        + 10 * S**6
        + 20 * S**5
        + 18 * S**4
        + 16 * S**3
        + 11 * S**2
        + 6 * S
        + 3
    )
    c2 = div_mod(c2_num, den1**2, modulus)

    n1 = (
        -4 * b * b * S**8
        - 16 * b * b * S**7
        - 16 * b * b * S**6
        + 16 * b * b * S**5
        + 40 * b * b * S**4
        + 16 * b * b * S**3
        - 16 * b * b * S**2
        - 16 * b * b * S
        - 4 * b * b
        + 16 * b * S**9
        + 64 * b * S**8
        + 112 * b * S**7
        + 128 * b * S**6
        + 128 * b * S**5
        + 100 * b * S**4
        + 16 * b * S**3
        - 56 * b * S**2
        - 48 * b * S
        - 12 * b
        + 36 * S**10
        + 144 * S**9
        + 248 * S**8
        + 272 * S**7
        + 224 * S**6
        + 140 * S**5
        + 56 * S**4
        - 12 * S**3
        - 48 * S**2
        - 36 * S
        - 9
    )
    c1 = div_mod(-n1, 4 * den1**4, modulus)

    a0 = 3 * S**4 + 9 * S**3 + 11 * S**2 + 9 * S + 3
    c0 = (div_mod(S**7 * a0, den1**3, modulus) + 2 * b * S**7) % modulus
    return c3, c2, c1, c0


def b_from_s_z_mod(s: int, z: int, eps: int, modulus: int) -> int:
    S = (eps * s) % modulus
    den1 = (S + 1) % modulus
    a0 = 3 * S**4 + 9 * S**3 + 11 * S**2 + 9 * S + 3
    c00 = div_mod(S**7 * a0, den1**3, modulus)
    lam = 2 * S**7
    return div_mod(z * z - c00, lam, modulus)


def surface_residue(s: int, u: int, z: int, eps: int, modulus: int) -> int:
    """Return the first-halving surface equation residue modulo modulus."""
    b = b_from_s_z_mod(s, z, eps, modulus)
    c3, c2, c1, c0 = coeffs_from_s_b_mod(s, b, eps, modulus)
    if (c0 - z * z) % modulus != 0:
        raise AssertionError("b formula failed to force c0=z^2")
    v = div_mod(u * u - c3, 2, modulus)
    w = div_mod(v * v + 2 * z - c2, 2 * u, modulus)
    return (c1 - (w * w - 2 * v * z)) % modulus


def lifts(base: int, modulus: int) -> Iterable[int]:
    step_count = modulus // 5
    for t in range(step_count):
        yield (base + 5 * t) % modulus


def run(max_power: int, sample_limit: int) -> None:
    print("contact-7 first-halving local p=5 diagnostic")
    print("open F_5 surface points:", BASE_POINTS)
    print("H4 divisible by 2 over J(F_5): false at all four points")
    print("Therefore second-halving lifts in these good open neighborhoods: 0")

    for power in range(1, max_power + 1):
        modulus = 5**power
        print(f"\nmodulus 5^{power} = {modulus}")
        total_all = 0
        for s0, u0, z0, eps in BASE_POINTS:
            total = 0
            samples: list[tuple[int, int, int]] = []
            for s in lifts(s0, modulus):
                for u in lifts(u0, modulus):
                    for z in lifts(z0, modulus):
                        if surface_residue(s, u, z, eps, modulus) == 0:
                            total += 1
                            if len(samples) < sample_limit:
                                samples.append((s, u, z))
            total_all += total
            print(
                "base",
                (s0, u0, z0, eps),
                "first_halving_lifts",
                total,
                "second_halving_lifts",
                0,
                "samples",
                samples,
            )
        print("total_first_halving_lifts", total_all, "total_second_halving_lifts", 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-power", type=int, default=3)
    parser.add_argument("--sample-limit", type=int, default=8)
    args = parser.parse_args()
    run(args.max_power, args.sample_limit)


if __name__ == "__main__":
    main()
