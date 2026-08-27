#!/usr/bin/env python3
"""Combine the three certified orbit-12 Hensel seeds by CRT.

The result is a residue class for the five projective CK coordinates with
the same marked pair {r1^2,r2^2}.  It satisfies both CK equations modulo the
full product; it is a lattice seed, not an exact rational point.
"""

from math import prod


MODULI = (11**3, 19**2, 23**2)
SEEDS = (
    (1, 242, 959, 9, 120),
    (1, 248, 98, 3, 11),
    (1, 392, 118, 8, 10),
)


def crt(residues, moduli):
    modulus = prod(moduli)
    value = 0
    for residue, q in zip(residues, moduli):
        cofactor = modulus // q
        value += residue * cofactor * pow(cofactor, -1, q)
    return value % modulus


def centered(value, modulus):
    value %= modulus
    return value if value <= modulus // 2 else value - modulus


def radicands(r, q):
    a = [(x * x) % q for x in r]
    return (
        (-(a[0] - a[2]) * (a[0] - a[3]) * (a[0] - a[4])) % q,
        ((a[2] - a[1]) * (a[0] - a[3]) * (a[0] - a[4])) % q,
        ((a[3] - a[1]) * (a[0] - a[2]) * (a[0] - a[4])) % q,
        ((a[4] - a[1]) * (a[0] - a[2]) * (a[0] - a[3])) % q,
    )


def main():
    modulus = prod(MODULI)
    combined = tuple(
        centered(crt([seed[i] for seed in SEEDS], MODULI), modulus)
        for i in range(5)
    )

    for q, seed in zip(MODULI, SEEDS):
        assert tuple(x % q for x in combined) == seed
        squares = {x * x % q for x in range(q)}
        assert all(g in squares for g in radicands(combined, q))

    linear = sum(combined)
    cubic = sum(x**3 for x in combined)
    assert linear % modulus == 0
    assert cubic % modulus == 0

    print("ELKIES22210_ORBIT12_CRT_SEED")
    print("moduli", MODULI)
    print("modulus_product", modulus)
    print("marked_pair", (1, 2))
    print("centered_residue", combined)
    print("linear_sum", linear, "quotient", linear // modulus)
    print("cubic_sum", cubic, "quotient", cubic // modulus)
    for q in MODULI:
        print("radicands_mod", q, radicands(combined, q))
    print("EXACT_POINT", False)
    print("DONE")


if __name__ == "__main__":
    main()
