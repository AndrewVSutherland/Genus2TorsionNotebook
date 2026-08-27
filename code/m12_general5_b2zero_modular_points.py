#!/usr/bin/env python3
"""Finite-field point triage for the b2=0 M(12)+5 norm cover.

No Groebner basis is formed.  The script exhausts p^5 tuples on either
the five-variable root-of-B quotient (generic linear B) or the separate
constant-B sign quotient.  It records raw, cheap-open, and fully open
counts, formal Jacobian ranks at fully open points, signed-double-cover
lifts, and the projection to compact-M(12) base pairs (b,w).

Examples:
  python3 code/m12_general5_b2zero_modular_points.py generic 7 11
  python3 code/m12_general5_b2zero_modular_points.py constant 7 11 13
"""

from argparse import ArgumentParser
from collections import Counter
from itertools import product
from math import comb

from m12_general5_b2zero_rootquotient_local import (
    compact_f,
    coprime,
    shift_scale,
)


def coeff(a, i):
    return a[i] if i < len(a) else 0


def p_add(a, b, p):
    n = max(len(a), len(b))
    return [((coeff(a, i) + coeff(b, i)) % p) for i in range(n)]


def p_scale(a, scalar, p):
    return [(scalar * x) % p for x in a]


def p_mul(a, b, p):
    ans = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            ans[i + j] = (ans[i + j] + x * y) % p
    return ans


def p_pow(a, n, p):
    ans = [1]
    base = a
    while n:
        if n & 1:
            ans = p_mul(ans, base, p)
        base = p_mul(base, base, p)
        n //= 2
    return ans


def p_sub(a, b, p):
    return p_add(a, p_scale(b, -1, p), p)


class Jet:
    """Value and first formal derivatives over F_p."""

    __slots__ = ("v", "g", "p")

    def __init__(self, value, gradient, p):
        self.p = p
        self.v = value % p
        self.g = tuple(x % p for x in gradient)

    def _coerce(self, other):
        if isinstance(other, Jet):
            return other
        return Jet(other, (0,) * len(self.g), self.p)

    def __add__(self, other):
        other = self._coerce(other)
        return Jet(self.v + other.v,
                   (a + b for a, b in zip(self.g, other.g)), self.p)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, (-a for a in self.g), self.p)

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, other):
        other = self._coerce(other)
        return Jet(self.v * other.v,
                   (self.v * b + other.v * a
                    for a, b in zip(self.g, other.g)), self.p)

    __rmul__ = __mul__

    def __pow__(self, n):
        if n == 0:
            return self._coerce(1)
        if n == 1:
            return self
        return Jet(pow(self.v, n, self.p),
                   (n * pow(self.v, n - 1, self.p) * a for a in self.g),
                   self.p)


def jvar(value, index, n, p):
    gradient = [0] * n
    gradient[index] = 1
    return Jet(value, gradient, p)


def j_padd(a, b):
    zero = a[0] * 0 if a else b[0] * 0
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else zero)
            + (b[i] if i < len(b) else zero) for i in range(n)]


def j_pscale(a, scalar):
    return [scalar * x for x in a]


def j_pmul(a, b):
    zero = a[0] * 0
    ans = [zero for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            ans[i + j] = ans[i + j] + x * y
    return ans


def j_ppow(a, n):
    one = a[0] * 0 + 1
    ans = [one]
    base = a
    while n:
        if n & 1:
            ans = j_pmul(ans, base)
        base = j_pmul(base, base)
        n //= 2
    return ans


def j_compact_f(b, w):
    one = b * 0 + 1
    L = [b, 2 * b - 1]
    H = [w, 1 + w * b]
    one_x_sq = [one, 2 * one, one]
    x_sq = [0 * one, 0 * one, one]
    inner = j_padd(
        j_pmul(L, j_pmul(H, H)),
        j_pscale(j_pmul([b], j_pmul(
            one_x_sq,
            j_padd(j_pmul([w], L), j_pscale(x_sq, -1)))), 4))
    return j_pmul(L, inner)


def j_shift_scale(f, c, d):
    zero = f[0] * 0
    g = [zero for _ in range(6)]
    for j, fj in enumerate(f):
        for i in range(j + 1):
            g[i] = g[i] + fj * comb(j, i) * (d ** i) * ((-c) ** (j - i))
    return g


def matrix_rank_mod(rows, p):
    a = [[x % p for x in row] for row in rows]
    rank = 0
    ncols = len(a[0]) if a else 0
    for col in range(ncols):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, p)
        a[rank] = [(inv * x) % p for x in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][col]:
                scalar = a[i][col]
                a[i] = [(x - scalar * y) % p
                        for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def generic_equations(b, w, c, d, e, p):
    f = compact_f(b, w, p)
    g = shift_scale(f, c, d, p)
    g0, g1, g2, g3, g4, g5 = g
    inv2 = pow(2, -1, p)
    inv8 = pow(8, -1, p)
    ell3 = ((5 * e * inv2) * g0 - g1) % p
    ell4 = ((5 * inv2 + 15 * e * e * inv8) * g0 - g2) % p
    ell5 = ((5 * inv2 + 15 * e * e * inv8) * g0 - g3) % p
    ell6 = ((5 * e * inv2) * g0 - g4) % p
    ell7 = (g0 - g5) % p
    delta = pow(2 - e, 3, p)
    equations = (
        ell3 - ell7,
        (32 * e * e - 33 * e + 58) * ell3 - 20 * ell5,
        (19 * e - 6) * ell3 - 8 * ell6,
        4 * (19 * e - 6) * ell3 * ell3
        - 32 * ell4 * ell3 + 5 * delta * g0 * g0,
    )
    return tuple(x % p for x in equations), f, g, ell3


def generic_jacobian(point, p):
    b, w, c, d, e = [jvar(x, i, 5, p) for i, x in enumerate(point)]
    f = j_compact_f(b, w)
    g0, g1, g2, g3, g4, g5 = j_shift_scale(f, c, d)
    inv2 = pow(2, -1, p)
    inv8 = pow(8, -1, p)
    ell3 = (5 * inv2) * e * g0 - g1
    ell4 = (5 * inv2 + (15 * inv8) * e**2) * g0 - g2
    ell5 = (5 * inv2 + (15 * inv8) * e**2) * g0 - g3
    ell6 = (5 * inv2) * e * g0 - g4
    ell7 = g0 - g5
    delta = (2 - e)**3
    eqs = [
        ell3 - ell7,
        (32 * e**2 - 33 * e + 58) * ell3 - 20 * ell5,
        (19 * e - 6) * ell3 - 8 * ell6,
        4 * (19 * e - 6) * ell3**2 - 32 * ell4 * ell3
        + 5 * delta * g0**2,
    ]
    assert all(eq.v == 0 for eq in eqs)
    return [list(eq.g) for eq in eqs]


def constant_residuals(b, w, u, v, s, p):
    f = compact_f(b, w, p)
    q = [v, u, 1]
    a = [0] * 6
    a[5] = 1
    inv2 = pow(2, -1, p)
    for degree in range(9, 4, -1):
        e_poly = p_sub(p_sub(p_mul(a, a, p), p_scale(f, s, p), p),
                       p_pow(q, 5, p), p)
        a[degree - 5] = (a[degree - 5]
                         - coeff(e_poly, degree) * inv2) % p
    e_poly = p_sub(p_sub(p_mul(a, a, p), p_scale(f, s, p), p),
                   p_pow(q, 5, p), p)
    return tuple(coeff(e_poly, i) % p for i in range(5)), f


def constant_jacobian(point, p):
    b, w, u, v, s = [jvar(x, i, 5, p) for i, x in enumerate(point)]
    f = j_compact_f(b, w)
    q = [v, u, b * 0 + 1]
    zero = b * 0
    a = [zero for _ in range(6)]
    a[5] = zero + 1
    inv2 = pow(2, -1, p)
    for degree in range(9, 4, -1):
        e_poly = j_padd(j_padd(j_pmul(a, a), j_pscale(f, -s)),
                        j_pscale(j_ppow(q, 5), -1))
        a[degree - 5] = a[degree - 5] - inv2 * e_poly[degree]
    e_poly = j_padd(j_padd(j_pmul(a, a), j_pscale(f, -s)),
                    j_pscale(j_ppow(q, 5), -1))
    eqs = e_poly[:5]
    assert all(eq.v == 0 for eq in eqs)
    return [list(eq.g) for eq in eqs]


def nonsingular_quintic(f, p):
    if len(f) != 6 or f[5] % p == 0:
        return False
    derivative = [(i * f[i]) % p for i in range(1, 6)]
    return coprime(f, derivative, p)


def report(label, p, raw, cheap, full, signed, rank_counts,
           base_fibers, samples):
    histogram = Counter(base_fibers.values())
    print(f"{label} p={p} raw={raw} cheap_open={cheap} full_open={full} "
          f"signed_lifts={signed} rank_counts={sorted(rank_counts.items())} "
          f"base_count={len(base_fibers)} "
          f"fiber_size_histogram={sorted(histogram.items())}")
    print(f"{label}_BASE_FIBERS p={p} {sorted(base_fibers.items())}")
    print(f"{label}_SAMPLES p={p} {samples}")


def enumerate_generic(p, sample_limit):
    squares = {x * x % p for x in range(1, p)}
    raw = cheap = full = signed = 0
    rank_counts = Counter()
    base_fibers = Counter()
    samples = []
    f_cache = {(b, w): compact_f(b, w, p)
               for b, w in product(range(p), repeat=2)}
    for b, w, c, d in product(range(p), repeat=4):
        f = f_cache[(b, w)]
        g = shift_scale(f, c, d, p)
        for e in range(p):
            eqs, _, _, ell3 = generic_equations(b, w, c, d, e, p)
            if any(eqs):
                continue
            raw += 1
            if (b == 0 or w == 0 or b == 1 or (2 * b - 1) % p == 0
                    or d == 0 or (e * e - 4) % p == 0 or ell3 == 0):
                continue
            cheap += 1
            if not nonsingular_quintic(f, p):
                continue
            if not coprime([1, e, 1], g, p):
                continue
            full += 1
            tau = (-5 * pow(2 - e, 3, p)
                   * pow(8 * ell3, -1, p)) % p
            is_signed = tau in squares
            if is_signed:
                signed += 2
            point = (b, w, c, d, e)
            rank = matrix_rank_mod(generic_jacobian(point, p), p)
            rank_counts[rank] += 1
            base_fibers[(b, w)] += 1
            if len(samples) < sample_limit:
                samples.append((point, rank, tau, is_signed))
    report("GENERIC_QUOTIENT", p, raw, cheap, full, signed, rank_counts,
           base_fibers, samples)


def enumerate_constant(p, sample_limit):
    squares = {x * x % p for x in range(1, p)}
    raw = cheap = full = signed = 0
    rank_counts = Counter()
    base_fibers = Counter()
    samples = []
    f_cache = {(b, w): compact_f(b, w, p)
               for b, w in product(range(p), repeat=2)}
    for b, w, u, v, s in product(range(p), repeat=5):
        eqs, f = constant_residuals(b, w, u, v, s, p)
        if any(eqs):
            continue
        raw += 1
        if b == 0 or w == 0 or b == 1 or (2 * b - 1) % p == 0 or s == 0:
            continue
        cheap += 1
        q = [v, u, 1]
        if not nonsingular_quintic(f, p):
            continue
        if (u * u - 4 * v) % p == 0 or not coprime(q, f, p):
            continue
        full += 1
        is_signed = s in squares
        if is_signed:
            signed += 2
        point = (b, w, u, v, s)
        rank = matrix_rank_mod(constant_jacobian(point, p), p)
        rank_counts[rank] += 1
        base_fibers[(b, w)] += 1
        if len(samples) < sample_limit:
            samples.append((point, rank, is_signed))
    report("CONSTANT_QUOTIENT", p, raw, cheap, full, signed, rank_counts,
           base_fibers, samples)


def main():
    parser = ArgumentParser()
    parser.add_argument("locus", choices=("generic", "constant"))
    parser.add_argument("primes", nargs="+", type=int)
    parser.add_argument("--sample-limit", type=int, default=12)
    args = parser.parse_args()
    for p in args.primes:
        if p in (2, 5):
            raise ValueError("primes must be odd and different from 5")
        if args.locus == "generic":
            enumerate_generic(p, args.sample_limit)
        else:
            enumerate_constant(p, args.sample_limit)


if __name__ == "__main__":
    main()
