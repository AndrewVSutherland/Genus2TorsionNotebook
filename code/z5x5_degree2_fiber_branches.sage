#!/usr/bin/env sage
"""Enumerate open roots of fixed finite [5,5] contact fibers.

This bounded scan covers every (u,v,r0,r1,r2) over the requested finite
field.  The monic square root, when it exists, uniquely determines H.
"""

import argparse
from itertools import product

from sage.all import GF, PolynomialRing, matrix


BASE_RESIDUES = {11: (0, 0, 4), 19: (1, 3, 16)}


def monic_square_root(poly, x):
    if poly.degree() != 10 or poly[10] != 1:
        return None
    field = poly.base_ring()
    root = x**5
    for degree in range(9, 4, -1):
        root += ((poly[degree] - (root**2)[degree]) / field(2))*x**(degree - 5)
    return root if root**2 == poly else None


def fiber_jacobian_determinant(prime, values):
    field = GF(prime)
    names = (
        "a", "b", "k", "u", "v", "r0", "r1", "r2",
        "s0", "s1", "s2", "s3", "s4",
    )
    ring = PolynomialRing(field, names, order="degrevlex")
    (a, b, k, u, v, r0, r1, r2,
     s0, s1, s2, s3, s4) = ring.gens()
    poly_ring = PolynomialRing(ring, "X")
    X = poly_ring.gen()
    f = (1 + a*X + b*X**2)**2 - k*X**5
    q = X**2 + u*X + v
    R = r0 + r1*X + r2*X**2
    H = X**5 + s4*X**4 + s3*X**3 + s2*X**2 + s1*X + s0
    identity = H**2 - f*R**2 - q**5
    equations = [ring(identity[i]) for i in range(10)]
    variables = ring.gens()[3:]
    jacobian = matrix(
        ring,
        [[equation.derivative(variable) for variable in variables]
         for equation in equations],
    )
    point = tuple(field(value) for value in BASE_RESIDUES[prime] + values)
    evaluated = jacobian.apply_map(lambda entry: field(entry(*point)))
    return evaluated.det()


def enumerate_branches(prime):
    field = GF(prime)
    poly_ring = PolynomialRing(field, "x")
    x = poly_ring.gen()
    a, b, k = map(field, BASE_RESIDUES[prime])
    f = (1 + a*x + b*x**2)**2 - k*x**5
    branches = []
    q_count = 0
    candidate_count = 0
    for u, v in product(field, repeat=2):
        q = x**2 + u*x + v
        if q.discriminant() == 0 or q.gcd(f).degree() > 0:
            continue
        q_count += 1
        for r0, r1, r2 in product(field, repeat=3):
            R = r0 + r1*x + r2*x**2
            if q.gcd(R).degree() > 0:
                continue
            candidate_count += 1
            H = monic_square_root(q**5 + f*R**2, x)
            if H is None:
                continue
            values = (
                u, v, r0, r1, r2,
                H[0], H[1], H[2], H[3], H[4],
            )
            determinant = fiber_jacobian_determinant(prime, values)
            branches.append((tuple(int(value) for value in values), int(determinant)))
    return f, q_count, candidate_count, branches


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="11,19")
    args = parser.parse_args()
    for prime in (int(value) for value in args.primes.split(",")):
        if prime not in BASE_RESIDUES:
            raise ValueError(f"no recorded base residue for p={prime}")
        f, q_count, candidate_count, branches = enumerate_branches(prime)
        smooth = [branch for branch in branches if branch[1] % prime]
        print(f"\n## p={prime}")
        print(f"base={BASE_RESIDUES[prime]}")
        print(f"f={f}")
        print(f"open_q_count={q_count}")
        print(f"candidate_q_R_count={candidate_count}")
        print(f"branch_count={len(branches)}")
        print(f"smooth_branch_count={len(smooth)}")
        for index, (values, determinant) in enumerate(smooth):
            print(f"branch_{index}={values}; fiber_det={determinant}")


if __name__ == "__main__":
    main()
