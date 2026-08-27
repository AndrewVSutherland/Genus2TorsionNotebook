#!/usr/bin/env python3
"""Affine good-reduction residue check for the q^2=1 diagonal.

This is a lightweight diagnostic, not a complete p-adic proof.  It checks
whether the simplified square classes

    F1 ~ -8*(d^2-1)/((d-n+1)*(d+n-1)),
    F4 ~  8*G(d,n)

can both be nonzero squares modulo p on the affine chart where the displayed
denominators and the D-surface denominator are nonzero.
"""

from __future__ import annotations


PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]


def inv_mod(a: int, p: int) -> int:
    return pow(a % p, -1, p)


def g_form(d: int, n: int, p: int) -> int:
    return (
        d**4
        + d * d * n * n
        + 2 * d * d * n
        - 2 * d * d
        + n * n
        - 2 * n
        + 1
    ) % p


def first_affine_solution(p: int):
    squares = {(a * a) % p for a in range(1, p)}
    for d in range(p):
        for n in range(p):
            den = (d * d + n * n - 1) % p
            left = (d - n + 1) % p
            right = (d + n - 1) % p
            if (
                d % p == 0
                or n % p == 0
                or (n - 1) % p == 0
                or (d - 1) % p == 0
                or (d + 1) % p == 0
                or den == 0
                or left == 0
                or right == 0
            ):
                continue
            g = g_form(d, n, p)
            if g == 0:
                continue
            f1 = (-8 * (d * d - 1) * inv_mod(left * right, p)) % p
            f4 = (8 * g) % p
            if f1 in squares and f4 in squares:
                return d, n, f1, f4
    return None


def main() -> None:
    print("M2248 q^2=1 diagonal affine local residue check")
    print("columns: p has_good_affine_residue witness")
    for p in PRIMES:
        sol = first_affine_solution(p)
        print(p, sol is not None, sol if sol is not None else "-")


if __name__ == "__main__":
    main()
