#!/usr/bin/env sage
"""Bounded two-prime Hensel/reconstruction search on the [5,5] contact cover.

The ten fiber variables are uniquely Hensel lifted above a fixed rational
base (a,b,k), starting from the smooth points recorded modulo 11 and 19 in
agent_z5x5_degree2_contact_probe.sage.  Bases are enumerated by a declared
height bound, and every apparent reconstruction is checked in QQ.
"""

import argparse
from itertools import product

from sage.all import (
    CRT,
    GF,
    HyperellipticCurve,
    Integer,
    PolynomialRing,
    QQ,
    ZZ,
    gcd,
    matrix,
    vector,
    Zmod,
)


NAMES = (
    "a", "b", "k", "u", "v", "r0", "r1", "r2",
    "s0", "s1", "s2", "s3", "s4",
)
FIBER_NAMES = NAMES[3:]
PRIMES = (11, 19)
BASE_RESIDUES = {
    11: (0, 0, 4),
    19: (1, 3, 16),
}
FIBER_RESIDUES = {
    11: (1, 7, 7, 0, 5, 2, 1, 4, 4, 2),
    19: (1, 3, 9, 0, 1, 18, 11, 10, 9, 4),
}


def contact_equations():
    ring = PolynomialRing(ZZ, NAMES, order="degrevlex")
    (a, b, k, u, v, r0, r1, r2,
     s0, s1, s2, s3, s4) = ring.gens()
    polynomial_ring = PolynomialRing(ring, "X")
    X = polynomial_ring.gen()
    f = (1 + a*X + b*X**2)**2 - k*X**5
    q = X**2 + u*X + v
    R = r0 + r1*X + r2*X**2
    H = X**5 + s4*X**4 + s3*X**3 + s2*X**2 + s1*X + s0
    identity = H**2 - f*R**2 - q**5
    equations = tuple(ring(identity[i]) for i in range(10))
    fiber_jacobian = matrix(
        ring,
        [[equation.derivative(variable) for variable in ring.gens()[3:]]
         for equation in equations],
    )
    return ring, equations, fiber_jacobian


RING, EQUATIONS, FIBER_JACOBIAN = contact_equations()


def evaluate(poly, values):
    return poly(*values)


def rational_mod(value, modulus_ring):
    value = QQ(value)
    return modulus_ring(value.numerator()) / modulus_ring(value.denominator())


def base_mod(base, modulus_ring):
    return tuple(rational_mod(value, modulus_ring) for value in base)


def fiber_inverse(prime, fiber_residue=None):
    field = GF(prime)
    if fiber_residue is None:
        fiber_residue = FIBER_RESIDUES[prime]
    values = tuple(field(value) for value in BASE_RESIDUES[prime] + fiber_residue)
    evaluated = FIBER_JACOBIAN.apply_map(lambda entry: field(evaluate(entry, values)))
    if not evaluated.is_invertible():
        raise ArithmeticError(f"fiber Jacobian is singular modulo {prime}")
    return evaluated.inverse()


FIBER_INVERSES = {prime: fiber_inverse(prime) for prime in PRIMES}


def monic_square_root(poly, x):
    if poly.degree() != 10 or poly[10] != 1:
        return None
    field = poly.base_ring()
    root = x**5
    for degree in range(9, 4, -1):
        root += ((poly[degree] - (root**2)[degree]) / field(2))*x**(degree - 5)
    return root if root**2 == poly else None


def enumerate_open_branches(prime):
    """Exhaust the fixed fiber by scanning all q and R over GF(prime)."""
    field = GF(prime)
    polynomial_ring = PolynomialRing(field, "x")
    x = polynomial_ring.gen()
    a, b, k = map(field, BASE_RESIDUES[prime])
    f = (1 + a*x + b*x**2)**2 - k*x**5
    branches = []
    for u, v in product(field, repeat=2):
        q = x**2 + u*x + v
        if q.discriminant() == 0 or q.gcd(f).degree() > 0:
            continue
        for r0, r1, r2 in product(field, repeat=3):
            R = r0 + r1*x + r2*x**2
            if q.gcd(R).degree() > 0:
                continue
            H = monic_square_root(q**5 + f*R**2, x)
            if H is None:
                continue
            residue = tuple(int(value) for value in (
                u, v, r0, r1, r2,
                H[0], H[1], H[2], H[3], H[4],
            ))
            try:
                inverse = fiber_inverse(prime, residue)
            except ArithmeticError:
                continue
            branches.append((residue, inverse))
    return branches


def hensel_lift(prime, precision, base, fiber_residue=None, inverse=None):
    """Lift the recorded fiber root modulo prime**precision digit by digit."""
    if fiber_residue is None:
        fiber_residue = FIBER_RESIDUES[prime]
    if inverse is None:
        inverse = fiber_inverse(prime, fiber_residue)
    fiber = [Integer(value) for value in fiber_residue]
    field = GF(prime)
    modulus = Integer(prime)

    for exponent in range(1, precision):
        next_modulus = modulus * prime
        residue_ring = Zmod(next_modulus)
        values = base_mod(base, residue_ring) + tuple(residue_ring(value) for value in fiber)
        quotients = []
        for equation in EQUATIONS:
            residual = Integer(evaluate(equation, values))
            if residual % modulus:
                raise ArithmeticError("input is not a root at the current precision")
            quotients.append(field(-(residual // modulus)))
        correction = inverse * vector(field, quotients)
        fiber = [
            Integer((value + modulus*Integer(delta)) % next_modulus)
            for value, delta in zip(fiber, correction)
        ]
        modulus = next_modulus

    residue_ring = Zmod(modulus)
    values = base_mod(base, residue_ring) + tuple(residue_ring(value) for value in fiber)
    if any(evaluate(equation, values) for equation in EQUATIONS):
        raise ArithmeticError("final Hensel residual is nonzero")
    return modulus, tuple(fiber)


def reconstruct_pair(left, right, left_modulus, right_modulus):
    modulus = left_modulus * right_modulus
    residue = Integer(CRT(left, right, left_modulus, right_modulus))
    try:
        return QQ(residue.rational_reconstruction(modulus))
    except (ArithmeticError, ValueError):
        return None


def reconstruct_fiber(lift11, lift19):
    modulus11, fiber11 = lift11
    modulus19, fiber19 = lift19
    reconstructed = tuple(
        reconstruct_pair(x11, x19, modulus11, modulus19)
        for x11, x19 in zip(fiber11, fiber19)
    )
    return reconstructed


def exact_residuals(base, fiber):
    if any(value is None for value in fiber):
        return None
    values = tuple(QQ(value) for value in base + fiber)
    return tuple(QQ(evaluate(equation, values)) for equation in EQUATIONS)


def rational_representatives(residue11, residue19, height_bound):
    modulus = PRIMES[0] * PRIMES[1]
    representatives = set()
    for denominator in range(1, height_bound + 1):
        if gcd(denominator, modulus) != 1:
            continue
        for numerator in range(-height_bound, height_bound + 1):
            if gcd(numerator, denominator) != 1:
                continue
            if (numerator - residue11*denominator) % PRIMES[0]:
                continue
            if (numerator - residue19*denominator) % PRIMES[1]:
                continue
            representatives.add(QQ(numerator) / denominator)
    return sorted(
        representatives,
        key=lambda value: (
            max(abs(value.numerator()), value.denominator()),
            abs(value.numerator()) + value.denominator(),
            value,
        ),
    )


def verify_hit(base, fiber):
    values = tuple(QQ(value) for value in base + fiber)
    (a, b, k, u, v, r0, r1, r2,
     s0, s1, s2, s3, s4) = values
    polynomial_ring = PolynomialRing(QQ, "x")
    x = polynomial_ring.gen()
    f = (1 + a*x + b*x**2)**2 - k*x**5
    q = x**2 + u*x + v
    R = r0 + r1*x + r2*x**2
    H = x**5 + s4*x**4 + s3*x**3 + s2*x**2 + s1*x + s0
    identity = H**2 - f*R**2 - q**5
    return {
        "f": f,
        "q": q,
        "R": R,
        "H": H,
        "identity": identity,
        "f_discriminant": f.discriminant(),
        "q_discriminant": q.discriminant(),
        "gcd_q_f": q.gcd(f),
        "gcd_q_R": q.gcd(R),
        "curve": HyperellipticCurve(f) if f.degree() == 5 and f.discriminant() else None,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height-bound", type=int, default=30)
    parser.add_argument("--precision", type=int, default=7)
    parser.add_argument("--max-bases", type=int, default=0,
                        help="0 means all bases within the height bound")
    parser.add_argument("--all-branches", action="store_true",
                        help="exhaust every smooth open branch in both fixed fibers")
    args = parser.parse_args()

    base_lists = []
    for index, name in enumerate(NAMES[:3]):
        values = rational_representatives(
            BASE_RESIDUES[11][index], BASE_RESIDUES[19][index], args.height_bound
        )
        base_lists.append(values)
        print(f"{name}_candidates={values}")

    bases = list(product(*base_lists))
    bases.sort(key=lambda base: (
        max(max(abs(x.numerator()), x.denominator()) for x in base),
        sum(abs(x.numerator()) + x.denominator() for x in base),
        base,
    ))
    if args.max_bases:
        bases = bases[:args.max_bases]
    print(f"height_bound={args.height_bound}")
    print(f"precision_bound={args.precision}")
    print(f"base_count={len(bases)}")
    if args.all_branches:
        branch_data = {prime: enumerate_open_branches(prime) for prime in PRIMES}
    else:
        branch_data = {
            prime: [(FIBER_RESIDUES[prime], FIBER_INVERSES[prime])]
            for prime in PRIMES
        }
    print("branch_counts=", {prime: len(branch_data[prime]) for prime in PRIMES})
    print(f"branch_pair_count={len(branch_data[11])*len(branch_data[19])}")

    exact_hits = []
    fully_reconstructed = 0
    tested_branch_pairs = 0
    for base_index, base in enumerate(bases, start=1):
        lifts = {
            prime: [
                hensel_lift(prime, args.precision, base, residue, inverse)
                for residue, inverse in branch_data[prime]
            ]
            for prime in PRIMES
        }
        for branch11_index, lift11 in enumerate(lifts[11]):
            for branch19_index, lift19 in enumerate(lifts[19]):
                tested_branch_pairs += 1
                fiber = reconstruct_fiber(lift11, lift19)
                residuals = exact_residuals(base, fiber)
                complete = residuals is not None
                if complete:
                    fully_reconstructed += 1
                if complete and not any(residuals):
                    verification = verify_hit(base, fiber)
                    exact_hits.append((base, fiber, verification))
                    print(f"EXACT_HIT_BASE_INDEX={base_index}")
                    print(f"branch_indices=({branch11_index},{branch19_index})")
                    print(f"base={base}")
                    print(f"fiber={fiber}")
                    for key, value in verification.items():
                        if key != "curve":
                            print(f"{key}={value}")

    print(f"tested_branch_pairs={tested_branch_pairs}")
    print(f"fully_reconstructed_at_final_precision={fully_reconstructed}")
    print(f"exact_hit_count={len(exact_hits)}")
    if not exact_hits:
        scope = "ALL_SMOOTH_OPEN_BRANCH_PAIRS" if args.all_branches else "RECORDED_BRANCH_PAIR"
        print(f"BOUNDED_RESULT=NO_EXACT_RATIONAL_FIBER_ON_{scope}")


if __name__ == "__main__":
    main()
