#!/usr/bin/env python3
"""Deep Hensel probes on the simplest open order-49 slices.

Slices tested:
  * v=0 at its two open F_11 points;
  * r=0 at its two open F_5 points.

In each case the remaining 4x4 Jacobian is invertible.  Lift to p^30,
rationally reconstruct under the uniqueness bound sqrt(p^30/2), and exact
test the four reduced incidence equations.
"""

from fractions import Fraction
from math import isqrt
import sys

sys.path.insert(0, "code")
import z49_structural_contact_iterate as contact
from z49_structural_3adic import (
    evaluate_all_mod, evaluate_exact, mat_vec, matrix_inverse_mod,
    rational_reconstruction_bounded, term_data,
)


def lift(polys, maxima, equations, variables, p, initial, fixed_index,
         fixed_value, precision):
    unknown = [i for i in range(5) if i != fixed_index]
    seed = {variable: value for variable, value in zip(variables, initial)}
    jacobian = [[int(poly.diff(variables[j]).subs(seed)) % p
                 for j in unknown] for poly in equations]
    inverse = matrix_inverse_mod(jacobian, p)
    assert inverse is not None
    current = [initial[j] % p for j in unknown]
    modulus = p
    for _ in range(1, precision):
        next_modulus = modulus*p
        values = []
        cursor = 0
        for j in range(5):
            if j == fixed_index:
                values.append(fixed_value % next_modulus)
            else:
                values.append(current[cursor])
                cursor += 1
        fvalues = evaluate_all_mod(polys, values, next_modulus, maxima)
        assert all(value % modulus == 0 for value in fvalues)
        rhs = [-(value//modulus) % p for value in fvalues]
        correction = mat_vec(inverse, rhs, p)
        current = [(current[i]+modulus*correction[i]) % next_modulus
                   for i in range(4)]
        modulus = next_modulus
    bound = isqrt(modulus//2)
    recon = [rational_reconstruction_bounded(value, modulus, bound)
             for value in current]
    exact = False
    if all(value is not None for value in recon):
        values = []
        cursor = 0
        for j in range(5):
            if j == fixed_index:
                values.append(Fraction(fixed_value))
            else:
                values.append(Fraction(*recon[cursor]))
                cursor += 1
        exact = all(evaluate_exact(poly, values) == 0 for poly in polys)
    return modulus, bound, recon, exact


def main():
    _x, variables, _h, _f, _solved, equations = contact.derive_system()
    polys = [term_data(poly) for poly in equations]
    maxima = [max(monomial[j] for terms in polys for monomial, _ in terms)
              for j in range(5)]
    cases = [
        ("v=0",11,(5,7,3,0,7),3,0),
        ("v=0",11,(5,7,8,0,7),3,0),
        ("r=0",5,(2,1,1,3,0),4,0),
        ("r=0",5,(2,1,4,2,0),4,0),
    ]
    print("Z49_STRUCTURAL_SLICE_HENSEL")
    for label,p,initial,index,value in cases:
        modulus,bound,recon,exact = lift(
            polys,maxima,equations,variables,p,initial,index,value,30)
        print("slice",label,"p",p,"seed",initial,"modulus",modulus,
              "reconstruction_bound",bound,"reconstruction",recon,
              "exact",exact)
    print("Z49_STRUCTURAL_SLICE_HENSEL_DONE")


if __name__ == "__main__":
    main()
