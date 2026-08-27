#!/usr/bin/env python3
"""Reconstruct and certify a closed form for the linked HLP formal branch.

This is deliberately a Padé reconstruction followed by exact substitution,
not a Groebner elimination.  Thirty formal coefficients are obtained from
19-by-19 rational linear solves.  Candidate rational functions are accepted
only if every coefficient is matched and all 18 defining polynomial
identities simplify identically to zero over Q(s).
"""

from fractions import Fraction as Q
import sys

import sympy as sp

sys.path.insert(0, "code")
import hlp_z60_linked_slice_formal as formal
from hlp_z60_marked_tangent import rank, rref


N = 30
NAMES = ["u0", "u1", "l1", "q50", "q51",
         "a0", "a1", "a2", "a3", "a4", "a5", "k5",
         "q30", "q31", "h0", "h1", "h2", "h3", "k3"]


def solve_linear(a, b):
    aug = [list(a[i]) + [b[i]] for i in range(len(a))]
    if rank(a) != len(a[0]) or rank(aug) != len(a[0]):
        return None
    _rr, pivots, red = rref(aug)
    out = [Q(0)] * len(a[0])
    for i, p in enumerate(pivots):
        if p < len(out):
            out[p] = red[i][-1]
    return out


def pade(c, max_degree=4, holdout=4):
    """Return minimal (numerator,denominator), denominator constant one."""
    for total in range(0, 2 * max_degree + 1):
        for d in range(0, min(max_degree, total) + 1):
            m = total - d
            if m > max_degree or m + d + 1 > len(c) - holdout:
                continue
            if d == 0:
                q = [Q(1)]
            else:
                rows, rhs = [], []
                for n in range(m + 1, m + d + 1):
                    rows.append([c[n - j] if n - j >= 0 else Q(0)
                                 for j in range(1, d + 1)])
                    rhs.append(-c[n])
                sol = solve_linear(rows, rhs)
                if sol is None:
                    continue
                q = [Q(1)] + sol
            p = []
            for n in range(m + 1):
                p.append(sum(q[j] * c[n - j]
                             for j in range(min(d, n) + 1)))
            good = True
            for n in range(m + 1, len(c)):
                if sum(q[j] * c[n - j]
                       for j in range(min(d, n) + 1)) != 0:
                    good = False
                    break
            if good:
                while len(p) > 1 and p[-1] == 0:
                    p.pop()
                while len(q) > 1 and q[-1] == 0:
                    q.pop()
                return p, q
    return None


def qsym(v):
    return sp.Rational(v.numerator, v.denominator)


def poly_sym(c, z):
    return sum(qsym(v) * z**i for i, v in enumerate(c))


def reconstruct(series, s):
    z = s**2
    ans = []
    shapes = []
    for name, c in zip(NAMES, series):
        odd = all(c[i] == 0 for i in range(0, len(c), 2))
        even = all(c[i] == 0 for i in range(1, len(c), 2))
        if not (odd or even):
            raise RuntimeError(f"mixed parity for {name}")
        seq = c[1::2] if odd else c[0::2]
        pq = pade(seq)
        if pq is None:
            ans.append(None)
            shapes.append(None)
            continue
        p, q = pq
        value = poly_sym(p, z) / poly_sym(q, z)
        if odd:
            value *= s
        ans.append(sp.factor(value))
        shapes.append((len(p) - 1, len(q) - 1, "odd" if odd else "even"))
    return ans, shapes


def coeffs(poly, x, degree):
    p = sp.Poly(sp.together(poly), x)
    return [p.coeff_monomial(x**i) for i in range(degree + 1)]


def verify(v, s):
    x = sp.symbols("x")
    (u0, u1, l1, q50, q51, a0, a1, a2, a3, a4, a5, k5,
     q30, q31, h0, h1, h2, h3, k3) = v
    q0 = x**2 - sp.Rational(1728, 125)
    u = x**2 + u1*x + u0
    line = -2775 + l1*x
    f = sp.expand((q0*line)**2 - 46250000*q0*u**2)
    q5 = x**2 + q51*x + q50
    aa = a5*x**5 + a4*x**4 + a3*x**3 + a2*x**2 + a1*x + a0
    q3 = x**2 + q31*x + q30
    hh = h3*x**3 + h2*x**2 + h1*x + h0
    e5 = sp.cancel(aa**2 - q0**2*f - k5*q5**5)
    e3 = sp.cancel(hh**2 - f - k3*q3**3)
    n5 = sp.Poly(sp.together(e5).as_numer_denom()[0], x)
    n3 = sp.Poly(sp.together(e3).as_numer_denom()[0], x)
    return n5.is_zero and n3.is_zero, sp.factor(f)


def main():
    free, series, _ff = formal.formal_marked_branch(N)
    assert free == 17
    s = sp.symbols("s")
    values, shapes = reconstruct(series, s)
    print("HLP_Z60_LINKED_SLICE_RECONSTRUCT")
    print("formal_order", N, "max_linear_system", "19x19")
    for name, shape in zip(NAMES, shapes):
        print(name, "shape", shape if shape is not None else "NONE_LE_4")
    if any(v is None for v in values):
        print("LOW_DEGREE_RATIONAL_PARAMETRIZATION", False)
        print("EXACT_18_IDENTITIES", "SKIPPED_NO_CANDIDATE")
    else:
        ok, _f = verify(values, s)
        print("LOW_DEGREE_RATIONAL_PARAMETRIZATION", ok)
        print("EXACT_18_IDENTITIES", ok)
    print("HLP_Z60_LINKED_SLICE_RECONSTRUCT_DONE")


if __name__ == "__main__":
    main()
