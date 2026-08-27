#!/usr/bin/env python3
"""Bounded Hensel search in the nonconstant-B root parameter.

For tau=v/u, evaluation of the norm identity at the root of B shows that

    (1+tau*r)/(1+tau) = t^2.

Writing t=1+tau*s and clearing u^2 gives the additional equation

    u^2*(r-1-2*s) - u*v*(s^2+2*s) - v^2*s^2 = 0.

Fixing rational s cuts the one-dimensional incidence curve to a finite
scheme.  Match all open F_5 and F_11 branches, Hensel lift the resulting
five-equation system, reconstruct bounded rational coordinates, and test
the original equations exactly.  The constant-B loci u*v=0 are excluded.
"""

import argparse
from fractions import Fraction
from math import gcd
import sys

sys.path.insert(0, "code")
import z49_structural_contact_iterate as contact
from z49_structural_3adic import (
    evaluate_all_mod, evaluate_exact, mat_vec, matrix_inverse_mod,
    rational_reconstruction_bounded, term_data,
)
from z49_structural_padic_complement import BRANCHES


# Complete open smooth incidence fiber, independently enumerated by the
# defining equations.  The u=0 and v=0 rows are retained here and excluded
# uniformly by the root-chart open conditions below.
BRANCHES[13] = [
    (1,2,1,4,10), (1,2,12,9,10),
    (1,4,1,0,0), (1,4,12,0,0),
    (8,5,0,3,3), (8,5,0,10,3),
    (9,10,1,7,7), (9,10,12,6,7),
    (9,11,6,12,4), (9,11,7,1,4),
    (9,12,1,5,2), (9,12,12,8,2),
]


def rational_parameters(height):
    return [
        (numerator, denominator)
        for denominator in range(1, height + 1)
        for numerator in range(-height, height + 1)
        if gcd(abs(numerator), denominator) == 1
    ]


def constraint(values, s):
    _a, _b, u, v, r = values
    return u*u*(r-1-2*s) - u*v*(s*s+2*s) - v*v*s*s


def constraint_jacobian(seed, s, p):
    _a, _b, u, v, r = seed
    return [
        0,
        0,
        2*u*(r-1-2*s)-v*(s*s+2*s),
        -u*(s*s+2*s)-2*v*s*s,
        u*u,
    ]


def lift(polys, maxima, inverse, initial, sfrac, p, precision):
    current = [value % p for value in initial]
    modulus = p
    numerator, denominator = sfrac
    for _ in range(1, precision):
        next_modulus = modulus*p
        sval = numerator*pow(denominator, -1, next_modulus) % next_modulus
        fvalues = evaluate_all_mod(polys, current, next_modulus, maxima)
        fvalues.append(constraint(current, sval) % next_modulus)
        assert all(value % modulus == 0 for value in fvalues)
        rhs = [-(value//modulus) % p for value in fvalues]
        correction = mat_vec(inverse, rhs, p)
        current = [
            (current[i]+modulus*correction[i]) % next_modulus
            for i in range(5)
        ]
        modulus = next_modulus
    return current, modulus


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parameter-height", type=int, default=50)
    parser.add_argument("--coordinate-bound", type=int, default=5000)
    args = parser.parse_args()
    precisions = {5: 12, 11: 8, 13: 7}
    assert all(
        p**precisions[p] > 2*args.coordinate_bound**2 for p in precisions
    )

    _x, variables, _h, _f, _solved, equations = contact.derive_system()
    polys = [term_data(poly) for poly in equations]
    maxima = [
        max(monomial[j] for terms in polys for monomial, _ in terms)
        for j in range(5)
    ]
    parameters = rational_parameters(args.parameter_height)
    attempted = singular = reconstructed = exact_hits = 0
    by_prime = {}

    for p, branches in BRANCHES.items():
        p_attempted = p_reconstructed = 0
        for sfrac in parameters:
            numerator, denominator = sfrac
            if denominator % p == 0:
                continue
            sres = numerator*pow(denominator, -1, p) % p
            for branch in branches:
                if branch[2] == 0 or branch[3] == 0:
                    continue
                if constraint(branch, sres) % p:
                    continue
                seed = {
                    variable: value for variable, value in zip(variables, branch)
                }
                jacobian = [
                    [int(poly.diff(variable).subs(seed)) % p
                     for variable in variables]
                    for poly in equations
                ]
                jacobian.append(constraint_jacobian(branch, sres, p))
                inverse = matrix_inverse_mod(jacobian, p)
                if inverse is None:
                    singular += 1
                    continue
                attempted += 1
                p_attempted += 1
                residues, modulus = lift(
                    polys, maxima, inverse, branch, sfrac, p, precisions[p]
                )
                rec = [
                    rational_reconstruction_bounded(
                        value, modulus, args.coordinate_bound
                    )
                    for value in residues
                ]
                if any(value is None for value in rec):
                    continue
                reconstructed += 1
                p_reconstructed += 1
                point = [Fraction(*value) for value in rec]
                sq = Fraction(*sfrac)
                if constraint(point, sq) != 0:
                    continue
                if not all(evaluate_exact(poly, point) == 0 for poly in polys):
                    continue
                a, b, u, v, r = point
                if (u == 0 or v == 0 or r == 1 or
                        a+b == Fraction(5, 2) or
                        u+v == 0 or u+v*r == 0):
                    continue
                exact_hits += 1
                print("EXACT_OPEN_HIT", "p", p, "s", sq, "point", point)
        by_prime[p] = (p_attempted, p_reconstructed)

    print("Z49_STRUCTURAL_ROOT_PARAMETER_HENSEL")
    print("parameter_height", args.parameter_height,
          "parameters", len(parameters),
          "coordinate_bound", args.coordinate_bound)
    for p in sorted(by_prime):
        print("prime", p, "branch_lifts", by_prime[p][0],
              "full_reconstructions", by_prime[p][1])
    print("attempted_branch_lifts", attempted,
          "singular_matches", singular,
          "all_five_reconstructed", reconstructed,
          "exact_open_hits", exact_hits)
    print("Z49_STRUCTURAL_ROOT_PARAMETER_HENSEL_DONE")


if __name__ == "__main__":
    main()
