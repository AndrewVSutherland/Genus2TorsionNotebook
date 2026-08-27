#!/usr/bin/env python3
"""Small Hensel/rational-reconstruction probe on order-49 contact slices.

The reduced iterated-contact system has four equations in (a,b,u,v,r).
At the recorded F_5 and F_11 points the first four Jacobian columns are
invertible, so fixing a small rational lift of r gives a unique p-adic lift.
This script lifts digit by digit and tests rational reconstruction exactly.
It is a diagnostic, not an exhaustive rational-point search.
"""

from math import gcd, isqrt
import sys

sys.path.insert(0, "code")
import z49_structural_contact_iterate as contact


def eval_terms(terms, values, modulus=None):
    total = 0
    for monomial, coefficient in terms:
        term = coefficient
        for value, exponent in zip(values, monomial):
            term *= value**exponent
        total += term
    return total if modulus is None else total % modulus


def solve_mod(matrix, rhs, p):
    n = len(matrix)
    aug = [[entry % p for entry in matrix[i]] + [rhs[i] % p]
           for i in range(n)]
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
                      for j in range(n+1)]
    return [aug[i][-1] for i in range(n)]


def rational_reconstruction(residue, modulus):
    bound = isqrt(modulus//2)
    r0, r1 = modulus, residue % modulus
    s0, s1 = 0, 1
    while abs(r1) > bound and r1:
        q = r0//r1
        r0, r1 = r1, r0-q*r1
        s0, s1 = s1, s0-q*s1
    if r1 == 0 or s1 == 0 or abs(s1) > bound:
        return None
    numerator, denominator = r1, s1
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    g = gcd(abs(numerator), denominator)
    numerator //= g
    denominator //= g
    if abs(numerator) > bound or denominator > bound:
        return None
    if (numerator-residue*denominator) % modulus:
        return None
    return numerator, denominator


def exact_zero(polys, fractions):
    from fractions import Fraction
    values = [Fraction(n, d) for n, d in fractions]
    return all(eval_terms(poly, values) == 0 for poly in polys)


def lift_branch(p, initial, fixed_r, precision, polys, jacobian):
    current = list(initial)
    modulus = p
    for exponent in range(1, precision):
        values = current + [fixed_r]
        fvalues = [eval_terms(poly, values) for poly in polys]
        assert all(value % modulus == 0 for value in fvalues)
        matrix = [[eval_terms(jacobian[i][j], values, p)
                   for j in range(4)] for i in range(4)]
        rhs = [-(value//modulus) for value in fvalues]
        correction = solve_mod(matrix, rhs, p)
        if correction is None:
            return None
        current = [current[i]+modulus*correction[i] for i in range(4)]
        modulus *= p
    recon = [rational_reconstruction(value, modulus) for value in current]
    exact = (all(item is not None for item in recon) and
             exact_zero(polys, recon + [(fixed_r, 1)]))
    return modulus, current, recon, exact


def main():
    _x, variables, _h, _f, _solved, equations = contact.derive_system()
    polys = [contact.term_data(poly) for poly in equations]
    jacobian = [[contact.term_data(poly.diff(variable))
                 for variable in variables[:4]] for poly in equations]
    branches = [
        (5, (2,1,1,3), 0),
        (5, (2,1,2,2), 2),
        (5, (2,1,3,3), 2),
        (5, (2,1,4,2), 0),
        (11, (5,7,3,0), -4),
        (11, (6,7,4,5), -5),
    ]
    print("Z49_STRUCTURAL_HENSEL")
    for p, initial, fixed_r in branches:
        result = lift_branch(p, initial, fixed_r, 30,
                             polys, jacobian)
        if result is None:
            print("branch",p,initial,"r",fixed_r,"SINGULAR")
            continue
        modulus, _current, recon, exact = result
        print("branch",p,initial,"r",fixed_r,
              "digits",30,"modulus",modulus,
              "reconstruction",recon,"exact",exact)
    print("Z49_STRUCTURAL_HENSEL_DONE")


if __name__ == "__main__":
    main()
