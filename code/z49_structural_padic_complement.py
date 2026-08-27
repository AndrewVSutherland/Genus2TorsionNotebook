#!/usr/bin/env python3
"""Complementary p=5 and p=11 reconstruction on the order-49 curve.

The reduced iterated-contact incidence has four open smooth points over F_5
and ten over F_11.  At each point the (a,b,u,v) Jacobian block is invertible,
so r is a local parameter.  This enumerates rational r of bounded height
integral at p, Hensel-lifts (a,b,u,v), reconstructs at the same bound, and
exact-tests the four defining equations.
"""

from fractions import Fraction
from math import gcd
import argparse
import sys

sys.path.insert(0, "code")
import z49_structural_contact_iterate as contact
from z49_structural_3adic import (
    evaluate_all_mod, evaluate_exact, mat_vec, matrix_inverse_mod,
    rational_reconstruction_bounded, term_data,
)


BRANCHES = {
    5: [
        (2,1,1,3,0), (2,1,2,2,2),
        (2,1,3,3,2), (2,1,4,2,0),
    ],
    11: [
        (5,7,3,0,7), (5,7,8,0,7),
        (6,7,4,5,6), (6,7,7,6,6),
        (7,6,2,2,9), (7,6,5,4,8),
        (7,6,6,7,8), (7,6,9,9,9),
        (10,1,3,7,3), (10,1,8,4,3),
    ],
}


def rational_parameters(height, p):
    by_residue = {i: [] for i in range(p)}
    for denominator in range(1, height+1):
        if denominator % p == 0:
            continue
        for numerator in range(-height, height+1):
            if gcd(abs(numerator), denominator) != 1:
                continue
            residue = numerator*pow(denominator, -1, p) % p
            by_residue[residue].append((numerator, denominator))
    return by_residue


def lift(polys, maxima, inverse, initial, rfrac, p, precision):
    current = list(initial)
    modulus = p
    numerator, denominator = rfrac
    for _ in range(1, precision):
        next_modulus = modulus*p
        rvalue = numerator*pow(denominator, -1, next_modulus) % next_modulus
        values = current + [rvalue]
        fvalues = evaluate_all_mod(polys, values, next_modulus, maxima)
        assert all(value % modulus == 0 for value in fvalues)
        rhs = [-(value//modulus) % p for value in fvalues]
        correction = mat_vec(inverse, rhs, p)
        current = [(current[i]+modulus*correction[i]) % next_modulus
                   for i in range(4)]
        modulus = next_modulus
    return current, modulus


def run_prime(p, height, equations, variables, polys, maxima):
    precision = 1
    while p**precision <= 2*height**2:
        precision += 1
    params = rational_parameters(height, p)
    lifted = reconstructed = exact_hits = 0
    branch_data = []
    for branch in BRANCHES[p]:
        values = {variable: value for variable, value in zip(variables, branch)}
        jacobian = [[int(poly.diff(variables[j]).subs(values)) % p
                     for j in range(4)] for poly in equations]
        inverse = matrix_inverse_mod(jacobian, p)
        assert inverse is not None
        branch_data.append((branch, inverse))

    for branch, inverse in branch_data:
        for rfrac in params[branch[4]]:
            residues, modulus = lift(polys, maxima, inverse, branch[:4],
                                     rfrac, p, precision)
            lifted += 1
            rec = [rational_reconstruction_bounded(value, modulus, height)
                   for value in residues]
            if any(value is None for value in rec):
                continue
            reconstructed += 1
            aq, bq, uq, vq = [Fraction(*value) for value in rec]
            rq = Fraction(*rfrac)
            point = [aq,bq,uq,vq,rq]
            if not all(evaluate_exact(poly, point) == 0 for poly in polys):
                continue
            if rq == 1 or aq+bq == Fraction(5,2) or uq+vq == 0 or uq+vq*rq == 0:
                continue
            exact_hits += 1
            print("EXACT_OPEN_HIT",p,point)
    param_count = sum(len(values) for values in params.values())
    print("prime",p,"height",height,"precision",precision,
          "modulus",p**precision,"rational_r",param_count,
          "branches",len(branch_data),"lifted",lifted,
          "all_four_reconstructed",reconstructed,
          "exact_open_hits",exact_hits)
    return exact_hits


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=200)
    parser.add_argument("--primes", default="5,11")
    args = parser.parse_args()
    _x, variables, _h, _f, _solved, equations = contact.derive_system()
    polys = [term_data(poly) for poly in equations]
    maxima = [max(monomial[j] for terms in polys for monomial, _ in terms)
              for j in range(5)]
    print("Z49_STRUCTURAL_PADIC_COMPLEMENT")
    total_hits = 0
    for p in [int(value) for value in args.primes.split(",")]:
        total_hits += run_prime(p,args.height,equations,variables,polys,maxima)
    print("total_exact_open_hits",total_hits)
    print("coverage p-integral listed smooth branches coordinate_height<=",
          args.height)
    print("Z49_STRUCTURAL_PADIC_COMPLEMENT_DONE")


if __name__ == "__main__":
    main()
