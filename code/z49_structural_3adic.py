#!/usr/bin/env python3
"""Bounded rational search on the two open 3-adic order-49 branches.

Modulo 3 the reduced iterated-contact system has exactly two open points,

    (a,b,u,v,r) = (1,1,+1,0,0), (1,1,-1,0,0).

The four Jacobian columns (a,u,v,r) are invertible at both, so b is a
3-adic local parameter.  Enumerate rational b of bounded naive height that
are integral and congruent to 1 mod 3, Hensel-lift the other four coordinates,
then rationally reconstruct them with the same height bound and exact-test
the original four equations.  Nonintegral 3-adic boundary charts are not
covered.
"""

from fractions import Fraction
from math import gcd
import argparse
import sys

sys.path.insert(0, "code")
import z49_structural_contact_iterate as contact


def rational_parameters(height):
    out = []
    for denominator in range(1, height+1):
        if denominator % 3 == 0:
            continue
        for numerator in range(-height, height+1):
            if gcd(abs(numerator), denominator) != 1:
                continue
            if numerator*pow(denominator, -1, 3) % 3 != 1:
                continue
            out.append((numerator, denominator))
    return out


def matrix_inverse_mod(matrix, p):
    n = len(matrix)
    aug = [[entry % p for entry in matrix[i]] +
           [1 if i == j else 0 for j in range(n)] for i in range(n)]
    for column in range(n):
        pivot = next((i for i in range(column, n)
                      if aug[i][column]), None)
        if pivot is None:
            return None
        aug[column], aug[pivot] = aug[pivot], aug[column]
        inv = pow(aug[column][column], -1, p)
        aug[column] = [(entry*inv) % p for entry in aug[column]]
        for i in range(n):
            if i == column:
                continue
            z = aug[i][column]
            aug[i] = [(aug[i][j]-z*aug[column][j]) % p
                      for j in range(2*n)]
    return [row[n:] for row in aug]


def mat_vec(matrix, vector, p):
    return [sum(matrix[i][j]*vector[j] for j in range(len(vector))) % p
            for i in range(len(matrix))]


def term_data(poly):
    return [(monomial, int(coefficient))
            for monomial, coefficient in poly.terms()]


def evaluate_all_mod(polys, values, modulus, maxima):
    powers = []
    for value, maximum in zip(values, maxima):
        row = [1]
        for _ in range(maximum):
            row.append(row[-1]*value % modulus)
        powers.append(row)
    answers = []
    for terms in polys:
        total = 0
        for monomial, coefficient in terms:
            term = coefficient
            for j, exponent in enumerate(monomial):
                if exponent:
                    term = term*powers[j][exponent] % modulus
            total += term
        answers.append(total % modulus)
    return answers


def evaluate_exact(terms, values):
    total = Fraction(0)
    for monomial, coefficient in terms:
        term = Fraction(coefficient)
        for value, exponent in zip(values, monomial):
            term *= value**exponent
        total += term
    return total


def rational_reconstruction_bounded(residue, modulus, bound):
    r0, r1 = modulus, residue % modulus
    s0, s1 = 0, 1
    while r1 and abs(r1) > bound:
        q = r0//r1
        r0, r1 = r1, r0-q*r1
        s0, s1 = s1, s0-q*s1
    if r1 == 0 or s1 == 0:
        return None
    numerator, denominator = r1, s1
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    common = gcd(abs(numerator), denominator)
    numerator //= common
    denominator //= common
    if abs(numerator) > bound or denominator > bound:
        return None
    if (numerator-residue*denominator) % modulus:
        return None
    return numerator, denominator


def lift(polys, maxima, inverse, branch_u, bfrac, precision):
    p = 3
    current = [1, branch_u % p, 0, 0]  # a,u,v,r
    modulus = p
    numerator, denominator = bfrac
    for _ in range(1, precision):
        next_modulus = modulus*p
        bvalue = numerator*pow(denominator, -1, next_modulus) % next_modulus
        values = [current[0], bvalue, current[1], current[2], current[3]]
        fvalues = evaluate_all_mod(polys, values, next_modulus, maxima)
        assert all(value % modulus == 0 for value in fvalues)
        rhs = [-(value//modulus) % p for value in fvalues]
        correction = mat_vec(inverse, rhs, p)
        current = [(current[i]+modulus*correction[i]) % next_modulus
                   for i in range(4)]
        modulus = next_modulus
    return current, modulus


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=100)
    parser.add_argument("--precision", type=int, default=12)
    args = parser.parse_args()
    assert 3**args.precision > 2*args.height**2

    _x, variables, _h, _f, _solved, equations = contact.derive_system()
    polys = [term_data(poly) for poly in equations]
    maxima = [max(monomial[j] for terms in polys for monomial, _ in terms)
              for j in range(5)]
    unknown_indices = [0,2,3,4]
    branches = []
    for branch_u in (1, -1):
        values = {variables[0]: 1, variables[1]: 1,
                  variables[2]: branch_u % 3,
                  variables[3]: 0, variables[4]: 0}
        jacobian = [[int(poly.diff(variables[j]).subs(values)) % 3
                     for j in unknown_indices] for poly in equations]
        inverse = matrix_inverse_mod(jacobian, 3)
        assert inverse is not None
        branches.append((branch_u, inverse))

    bvalues = rational_parameters(args.height)
    lifted = reconstructed = exact_hits = 0
    hits = []
    for bfrac in bvalues:
        bq = Fraction(*bfrac)
        for branch_u, inverse in branches:
            residues, modulus = lift(polys, maxima, inverse, branch_u,
                                     bfrac, args.precision)
            lifted += 1
            rec = [rational_reconstruction_bounded(value, modulus,
                                                   args.height)
                   for value in residues]
            if any(value is None for value in rec):
                continue
            reconstructed += 1
            aq, uq, vq, rq = [Fraction(*value) for value in rec]
            values = [aq, bq, uq, vq, rq]
            if not all(evaluate_exact(poly, values) == 0 for poly in polys):
                continue
            # Exact open conditions in the reduced chart.
            if rq == 1 or aq+bq == Fraction(5,2) or uq+vq == 0 or uq+vq*rq == 0:
                continue
            exact_hits += 1
            hits.append(values)
            print("EXACT_OPEN_HIT",values)

    print("Z49_STRUCTURAL_3ADIC")
    print("height",args.height,"precision",args.precision,
          "modulus",3**args.precision)
    print("rational_b",len(bvalues),"branches",2,"lifted",lifted)
    print("all_four_reconstructed",reconstructed,
          "exact_open_hits",exact_hits)
    print("coverage integral_at_3 b=1_mod3 coordinate_height<=",args.height)
    print("Z49_STRUCTURAL_3ADIC_DONE")


if __name__ == "__main__":
    main()
